/*
 * mhaptools convert: convert bisulfite sequencing alignments into mHap haplotypes.
 *
 * Rewritten in C on the read-level rules of mhapasm (bsread.c, from
 * https://github.com/JiantaoShi/mHapASM): reads are laid out on the reference with their
 * CIGAR (soft clips, insertions, deletions), the bisulfite strand comes from
 * the aligner's tags (Bismark XG, bwa-meth/BISCUIT YD, BSMAP ZS) or the SAM
 * flags, mates are merged into fragments, and CpGs are called with the
 * wgbs_tools rules (CpG context checked on the read).
 *
 * Output: mHap records (chr, start, end, haplotype, count, strand), where
 * start/end are the 1-based positions of the first/last CpG C and the
 * haplotype has one 0/1 per consecutive CpG. Records are sorted, identical
 * records collapsed, bgzipped and tabix-indexed (columns 1-3, 1-based).
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>
#include <time.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include "bsread.h"

#include "convert.h"
#include "version.h"

KHASH_MAP_INIT_STR(cpgchr, int)

typedef struct {
    const char *in_fn, *cpg_fn, *out_fn, *region, *bed_fn, *ref_fn;
    int min_mapq, min_bq, excl_flags, max_frag, max_ch, threads;
    double max_unconv;
    int non_dir, taps, split, qname, index;
} opt_t;

typedef struct {
    hts_pos_t start, end;   /* 0-based positions of the first and last CpG C */
    char *hap;              /* one '0'/'1' per CpG */
    char *qname;            /* --qname only */
    char strand;
} rec_t;

typedef struct {
    opt_t *o;
    samFile *fp;
    sam_hdr_t *hdr;
    hts_itr_t *itr;          /* multi-region iterator for -r/-b */
    bam1_t *b;
    /* CpGs: from a tabix-indexed file per contig, or all loaded at start */
    htsFile *cpg_fp;
    tbx_t *tbx;
    khash_t(cpgchr) *cpg_names;
    cpgs_t *cpg_all;
    int n_cpg_all;
    cpgs_t cpg;              /* CpGs of the current contig */
    int tid;                 /* current contig */
    bs_pairer_t *pairer;
    bs_buf_t buf;
    rec_t *rec;              /* records not written yet */
    int64_t nrec, mrec, next_flush;
    BGZF *out;
    kstring_t ks, line;
    kstring_t no_cpg_chr;    /* contigs absent from the CpG file */
    int n_no_cpg_chr;
    /* run summary */
    int64_t n_in, n_filtered, n_ch, n_reads, n_frags, n_paired, n_nocpg, n_dropped, n_split;
    int64_t n_records, n_lines;
} cvt_t;

/***************************************************************
 *                        CpG positions                        *
 ***************************************************************/

/* Without a tabix index, the whole CpG file is loaded once. */
static int load_all_cpgs(cvt_t *c, const char *fn)
{
    htsFile *fp = hts_open(fn, "r");
    if (!fp) { fprintf(stderr, "[mhaptools convert] cannot open %s\n", fn); return -1; }
    c->cpg_names = kh_init(cpgchr);
    int cur = -1, ret;
    while (hts_getline(fp, '\n', &c->ks) >= 0) {
        char *s = c->ks.s, *tab = strchr(s, '\t'), *q;
        if (!tab || s[0] == '#') continue;
        *tab = '\0';
        long long pos = strtoll(tab + 1, &q, 10);
        if (q == tab + 1 || pos < 1) continue;
        if (cur < 0 || strcmp(kh_key(c->cpg_names, (khint_t) cur), s)) {
            khint_t k = kh_get(cpgchr, c->cpg_names, s);
            if (k == kh_end(c->cpg_names)) {
                k = kh_put(cpgchr, c->cpg_names, bs_strdup(s), &ret);
                kh_val(c->cpg_names, k) = c->n_cpg_all;
                c->cpg_all = bs_realloc(c->cpg_all, (c->n_cpg_all + 1) * sizeof(cpgs_t));
                memset(&c->cpg_all[c->n_cpg_all++], 0, sizeof(cpgs_t));
            }
            cur = (int) k;
        }
        cpgs_t *g = &c->cpg_all[kh_val(c->cpg_names, (khint_t) cur)];
        if (g->n == g->m) {
            g->m = g->m ? g->m * 2 : 4096;
            g->pos = bs_realloc(g->pos, g->m * sizeof(hts_pos_t));
        }
        g->pos[g->n++] = pos - 1;
    }
    hts_close(fp);
    for (int i = 0; i < c->n_cpg_all; i++) {   /* sort, drop duplicates, clear the mask */
        cpgs_t *g = &c->cpg_all[i];
        int sorted = 1, j = 0;
        for (int k = 1; k < g->n && sorted; k++) sorted = g->pos[k] >= g->pos[k - 1];
        if (!sorted) {
            for (int k = 1; k < g->n; k++) {   /* insertion sort, CpG files are nearly sorted */
                hts_pos_t v = g->pos[k];
                int l = k - 1;
                while (l >= 0 && g->pos[l] > v) { g->pos[l + 1] = g->pos[l]; l--; }
                g->pos[l + 1] = v;
            }
        }
        for (int k = 0; k < g->n; k++)
            if (!j || g->pos[k] != g->pos[j - 1]) g->pos[j++] = g->pos[k];
        g->n = j;
        g->mask = bs_malloc(g->n ? g->n : 1);
        memset(g->mask, 0, g->n);
    }
    return 0;
}

static int cpg_lookup(void *h, const char *name)
{
    khint_t k = kh_get(cpgchr, (khash_t(cpgchr) *) h, name);
    return k == kh_end((khash_t(cpgchr) *) h) ? -1 : kh_val((khash_t(cpgchr) *) h, k);
}

/* Make the CpGs of BAM contig `tid` current. */
static void load_contig_cpgs(cvt_t *c, int tid)
{
    const char *chr = sam_hdr_tid2name(c->hdr, tid);
    int n;
    if (c->tbx) {
        n = bs_cpgs_load(c->cpg_fp, c->tbx, chr, 0, HTS_POS_MAX, &c->cpg, &c->ks);
    } else {
        /* borrow the preloaded array; same "chr" prefix fallback as bsread */
        int id = cpg_lookup(c->cpg_names, chr);
        if (id < 0) {
            kstring_t alt = {0, 0, NULL};
            if (!strncmp(chr, "chr", 3)) kputs(chr + 3, &alt); else ksprintf(&alt, "chr%s", chr);
            id = cpg_lookup(c->cpg_names, alt.s);
            if (id < 0 && (!strcmp(chr, "MT") || !strcmp(chr, "chrM"))) id = cpg_lookup(c->cpg_names, !strcmp(chr, "MT") ? "chrM" : "MT");
            free(alt.s);
        }
        if (id >= 0) c->cpg = c->cpg_all[id];
        else memset(&c->cpg, 0, sizeof(cpgs_t));
        n = id < 0 ? -1 : c->cpg.n;
    }
    if (n < 0) {   /* reported once at the end */
        if (c->n_no_cpg_chr++ < 5) ksprintf(&c->no_cpg_chr, "%s%s", c->no_cpg_chr.l ? ", " : "", chr);
    }
}

/***************************************************************
 *                        records                              *
 ***************************************************************/

static int rec_cmp(const void *a, const void *b)
{
    const rec_t *x = a, *y = b;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    if (x->end != y->end) return x->end < y->end ? -1 : 1;
    int d = strcmp(x->hap, y->hap);
    if (d) return d;
    if (x->strand != y->strand) return x->strand < y->strand ? -1 : 1;
    if (x->qname && y->qname) return strcmp(x->qname, y->qname);
    return 0;
}

static int write_line(cvt_t *c)
{
    if (bgzf_write(c->out, c->line.s, c->line.l) < 0) {
        fprintf(stderr, "[mhaptools convert] error writing %s\n", c->o->out_fn);
        exit(1);
    }
    c->n_lines++;
    return 0;
}

/* Write, in sorted order and with identical records collapsed, every record
 * that starts before `below`; no record created later can start before it. */
static void flush_records(cvt_t *c, hts_pos_t below)
{
    if (!c->nrec) return;
    qsort(c->rec, c->nrec, sizeof(rec_t), rec_cmp);
    const char *chr = sam_hdr_tid2name(c->hdr, c->tid);
    int64_t i = 0;
    while (i < c->nrec && c->rec[i].start < below) {
        int64_t j = i + 1;
        if (!c->o->qname)
            while (j < c->nrec && rec_cmp(&c->rec[i], &c->rec[j]) == 0) j++;
        rec_t *r = &c->rec[i];
        c->line.l = 0;
        ksprintf(&c->line, "%s\t%" PRIhts_pos "\t%" PRIhts_pos "\t%s\t%" PRId64 "\t%c",
                 chr, r->start + 1, r->end + 1, r->hap, j - i, r->strand);
        if (c->o->qname) ksprintf(&c->line, "\t%s", r->qname);
        kputc('\n', &c->line);
        write_line(c);
        for (int64_t k = i; k < j; k++) { free(c->rec[k].hap); free(c->rec[k].qname); }
        i = j;
    }
    memmove(c->rec, c->rec + i, (c->nrec - i) * sizeof(rec_t));
    c->nrec -= i;
}

static void add_record(cvt_t *c, const frag_t *f, int lo, int hi, char strand)
{
    if (c->nrec == c->mrec) {
        c->mrec = c->mrec ? c->mrec * 2 : 65536;
        c->rec = bs_realloc(c->rec, c->mrec * sizeof(rec_t));
    }
    rec_t *r = &c->rec[c->nrec++];
    r->start = c->cpg.pos[f->cpg[lo]];
    r->end = c->cpg.pos[f->cpg[hi - 1]];
    r->strand = strand;
    r->hap = bs_malloc(hi - lo + 1);
    for (int k = lo; k < hi; k++) r->hap[k - lo] = (f->meth[k] ^ c->o->taps) ? '1' : '0';
    r->hap[hi - lo] = '\0';
    r->qname = c->o->qname ? bs_strdup(f->qname) : NULL;
    c->n_records++;
}

/* 1 if a CpG in [ia, ib] lies inside the aligned part of a mate, i.e. it was
 * sequenced but got no call (error, failed context, deletion, mate conflict). */
static int gap_sequenced(const cvt_t *c, const frag_t *f, int32_t ia, int32_t ib)
{
    for (int32_t i = ia; i <= ib; i++)
        for (int k = 0; k < f->nmate; k++)
            if (c->cpg.pos[i] >= f->m[k].beg && c->cpg.pos[i] + 1 < f->m[k].end) return 1;
    return 0;
}

static int frag_bottom(const frag_t *f)
{
    for (int k = 0; k < f->nmate; k++) if (f->m[k].read1) return f->m[k].bottom;
    return f->m[0].bottom;
}

/* A fragment is complete: call its CpGs and turn the calls into mHap records.
 * An mHap haplotype covers consecutive CpGs, so the calls are cut where CpGs
 * are missing. CpGs between non-overlapping mates were not sequenced and the
 * fragment is cut there (as in mHapTools); a sequenced CpG without a call
 * drops the fragment, or cuts it with --split. */
static void frag_done(frag_t *f, void *data)
{
    cvt_t *c = data;
    c->n_frags++;
    if (f->nmate == 2) c->n_paired++;
    bs_frag_call(f, &c->cpg, c->o->min_bq, &c->buf);
    if (f->ncall == 0) { c->n_nocpg++; bs_frag_free(f); return; }
    char strand = c->o->non_dir ? '*' : frag_bottom(f) ? '-' : '+';
    if (!c->o->split) {
        for (int k = 1; k < f->ncall; k++) {
            if (f->cpg[k] > f->cpg[k - 1] + 1 && gap_sequenced(c, f, f->cpg[k - 1] + 1, f->cpg[k] - 1)) {
                c->n_dropped++;
                bs_frag_free(f);
                return;
            }
        }
    }
    int lo = 0, pieces = 0;
    for (int k = 1; k <= f->ncall; k++) {
        if (k < f->ncall && f->cpg[k] == f->cpg[k - 1] + 1) continue;
        add_record(c, f, lo, k, strand);
        pieces++;
        lo = k;
    }
    if (pieces > 1) c->n_split++;
    bs_frag_free(f);
}

/***************************************************************
 *                        driver                               *
 ***************************************************************/

/* All reads of the current contig are in: write everything out. */
static void finish_contig(cvt_t *c)
{
    if (c->tid < 0) return;
    bs_pairer_flush(c->pairer, HTS_POS_MAX);
    flush_records(c, HTS_POS_MAX);
    if (c->tbx) bs_cpgs_free(&c->cpg);
    else memset(&c->cpg, 0, sizeof(cpgs_t));
    c->tid = -1;
}

static int next_read(cvt_t *c)
{
    return c->itr ? sam_itr_next(c->fp, c->itr, c->b) : sam_read1(c->fp, c->hdr, c->b);
}

static int run(cvt_t *c)
{
    const opt_t *o = c->o;
    hts_pos_t last_pos = -1;
    int ret;
    c->tid = -1;
    c->next_flush = 65536;
    while ((ret = next_read(c)) >= 0) {
        bam1_t *b = c->b;
        if (b->core.tid < 0) break;          /* unmapped reads come last */
        c->n_in++;
        if (b->core.tid != c->tid) {
            if (c->tid >= 0 && b->core.tid < c->tid) {
                fprintf(stderr, "[mhaptools convert] %s is not sorted by coordinate\n", o->in_fn);
                return -1;
            }
            finish_contig(c);
            c->tid = b->core.tid;
            last_pos = -1;
            load_contig_cpgs(c, c->tid);
        }
        hts_pos_t pos = b->core.pos;
        if (pos < last_pos) {
            fprintf(stderr, "[mhaptools convert] %s is not sorted by coordinate\n", o->in_fn);
            return -1;
        }
        last_pos = pos;
        bs_pairer_flush(c->pairer, pos);     /* mates that never showed up */
        if (c->nrec >= c->next_flush) {
            flush_records(c, pos - o->max_frag);
            c->next_flush = c->nrec + 65536;
        }
        if ((b->core.flag & o->excl_flags) || b->core.qual < o->min_mapq || c->cpg.n == 0) {
            c->n_filtered++;
            continue;
        }
        mate_t m;
        if (bs_mate_init(&m, b) < 0) { c->n_filtered++; continue; }
        if (!o->taps && (o->max_ch >= 0 || o->max_unconv >= 0)) {
            int nc, nt;
            bs_mate_ch(&m, &c->cpg, &nc, &nt);
            if ((o->max_ch >= 0 && nc > o->max_ch) ||
                (o->max_unconv >= 0 && nc >= 3 && nc > o->max_unconv * (nc + nt))) {
                c->n_ch++;
                bs_mate_free(&m);
                continue;
            }
        }
        c->n_reads++;
        bs_pairer_add(c->pairer, b, &m);
    }
    if (ret < -1) { fprintf(stderr, "[mhaptools convert] error reading %s\n", o->in_fn); return -1; }
    finish_contig(c);
    return 0;
}

/* Region strings for the multi-region iterator, from -r or a BED file. */
static char **read_regions(const opt_t *o, sam_hdr_t *hdr, int *n_out)
{
    char **regs = NULL;
    int n = 0, m = 0;
    if (o->region) {
        regs = bs_malloc(sizeof(char *));
        regs[n++] = bs_strdup(o->region);
        *n_out = n;
        return regs;
    }
    htsFile *fp = hts_open(o->bed_fn, "r");
    if (!fp) { fprintf(stderr, "[mhaptools convert] cannot open %s\n", o->bed_fn); return NULL; }
    kstring_t ks = {0, 0, NULL};
    int64_t n_bad = 0;
    while (hts_getline(fp, '\n', &ks) >= 0) {
        if (!ks.l || ks.s[0] == '#' || !strncmp(ks.s, "track", 5) || !strncmp(ks.s, "browser", 7)) continue;
        char *chr = strtok(ks.s, " \t"), *sb = strtok(NULL, " \t"), *se = strtok(NULL, " \t\r"), *e1, *e2;
        if (!chr || !sb || !se) { n_bad++; continue; }
        long long beg = strtoll(sb, &e1, 10), end = strtoll(se, &e2, 10);
        int tid = bs_bam_tid(hdr, chr);
        if (*e1 || *e2 || beg < 0 || end <= beg || tid < 0) { n_bad++; continue; }
        if (n == m) { m = m ? m * 2 : 1024; regs = bs_realloc(regs, m * sizeof(char *)); }
        kstring_t r = {0, 0, NULL};
        ksprintf(&r, "{%s}:%lld-%lld", sam_hdr_tid2name(hdr, tid), beg + 1, end);
        regs[n++] = r.s;
    }
    free(ks.s);
    hts_close(fp);
    if (n_bad) fprintf(stderr, "[mhaptools convert] skipped %" PRId64 " BED line(s) (malformed or contig absent from the BAM)\n", n_bad);
    *n_out = n;
    return regs;
}

static void usage(FILE *fp)
{
    fprintf(fp,
"mhaptools convert (mHapTools " HAP_VERSION_TEXT "): convert bisulfite sequencing alignments into mHap haplotypes\n"
"\n"
"Usage: mhaptools convert -i <in.bam> -c <CpG.gz> [-r chr:beg-end | -b regions.bed] [-o out.mhap.gz] [options]\n"
"\n"
"Input:\n"
"  -i, --input FILE        coordinate-sorted SAM/BAM/CRAM ('-' for stdin); -r and -b need an index\n"
"  -c, --cpg FILE          CpG positions (chr, 1-based C position, ...); plain or (b)gzipped,\n"
"                          read per contig when tabix-indexed, otherwise loaded at once\n"
"  -r, --region STR        only reads overlapping this region (chr:beg-end)\n"
"  -b, --bed FILE          only reads overlapping these BED regions\n"
"  -T, --reference FILE    reference FASTA (CRAM input only)\n"
"Output:\n"
"  -o, --output FILE       mHap file [out.mhap.gz]; bgzipped and tabix-indexed if it ends in .gz\n"
"  -n, --non-directional   report strand '*' (haplotypes not grouped by strand)\n"
"  -m, --mode STR          BS (bisulfite, EM-seq: C = methylated) or TAPS (T = methylated) [BS]\n"
"      --qname             one line per record with the read name in column 7 (no collapsing)\n"
"      --no-index          do not write the tabix index\n"
"Filters:\n"
"  -q, --min-mapq INT      minimum mapping quality [10]\n"
"  -F, --exclude-flags INT skip reads with any of these flags [0xF04]\n"
"  -B, --min-bq INT        minimum base quality of a CpG call [0]\n"
"  -L, --max-frag INT      maximum mate distance of a fragment, in bp [1000]\n"
"      --max-unconv FLOAT  drop incompletely converted reads: at least 3 of their\n"
"                          cytosines outside CpGs are unconverted and they make up more\n"
"                          than FLOAT of the C+T (G+A on the bottom strand) there;\n"
"                          -1 = no filter; ignored with -m TAPS [0.2]\n"
"      --max-ch INT        drop reads with more than INT unconverted cytosines outside\n"
"                          CpGs (counts library-prep tails and normal CH methylation\n"
"                          too); -1 = no filter [-1]\n"
"      --split             cut a fragment at a sequenced CpG without a call instead of\n"
"                          dropping the fragment\n"
"Other:\n"
"  -@, --threads INT       extra threads for decompression and compression [0]\n"
"  -h, --help              show this help\n"
"  -v, --version           show the version\n");
}

int main_convert(int argc, char *argv[])
{
    bs_prog = "mhaptools convert";
    opt_t o;
    memset(&o, 0, sizeof(o));
    o.out_fn = "out.mhap.gz";
    o.min_mapq = 10;
    o.excl_flags = 0xF04;
    o.max_frag = 1000;
    o.max_ch = -1;
    o.max_unconv = 0.2;
    o.index = 1;

    static const struct option lopts[] = {
        {"input", required_argument, NULL, 'i'},
        {"cpg", required_argument, NULL, 'c'},
        {"region", required_argument, NULL, 'r'},
        {"bed", required_argument, NULL, 'b'},
        {"reference", required_argument, NULL, 'T'},
        {"output", required_argument, NULL, 'o'},
        {"non-directional", no_argument, NULL, 'n'},
        {"mode", required_argument, NULL, 'm'},
        {"qname", no_argument, NULL, 1},
        {"no-index", no_argument, NULL, 2},
        {"max-ch", required_argument, NULL, 3},
        {"max-unconv", required_argument, NULL, 5},
        {"split", no_argument, NULL, 4},
        {"min-mapq", required_argument, NULL, 'q'},
        {"exclude-flags", required_argument, NULL, 'F'},
        {"min-bq", required_argument, NULL, 'B'},
        {"max-frag", required_argument, NULL, 'L'},
        {"threads", required_argument, NULL, '@'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}
    };
    int ch;
    while ((ch = getopt_long(argc, argv, "i:c:r:b:T:o:nm:q:F:B:L:@:hv", lopts, NULL)) >= 0) {
        switch (ch) {
        case 'i': o.in_fn = optarg; break;
        case 'c': o.cpg_fn = optarg; break;
        case 'r': o.region = optarg; break;
        case 'b': o.bed_fn = optarg; break;
        case 'T': o.ref_fn = optarg; break;
        case 'o': o.out_fn = optarg; break;
        case 'n': o.non_dir = 1; break;
        case 'm':
            if (!strcasecmp(optarg, "TAPS")) o.taps = 1;
            else if (strcasecmp(optarg, "BS")) { fprintf(stderr, "[mhaptools convert] -m must be BS or TAPS\n"); return 1; }
            break;
        case 1: o.qname = 1; break;
        case 2: o.index = 0; break;
        case 3: o.max_ch = atoi(optarg); break;
        case 5: o.max_unconv = atof(optarg); break;
        case 4: o.split = 1; break;
        case 'q': o.min_mapq = atoi(optarg); break;
        case 'F': o.excl_flags = (int) strtol(optarg, NULL, 0); break;
        case 'B': o.min_bq = atoi(optarg); break;
        case 'L': o.max_frag = atoi(optarg); break;
        case '@': o.threads = atoi(optarg); break;
        case 'h': usage(stdout); return 0;
        case 'v': puts(HAP_VERSION_TEXT); return 0;
        default: usage(stderr); return 1;
        }
    }
    if (!o.in_fn || !o.cpg_fn || optind != argc) { usage(stderr); return 1; }
    if (o.region && o.bed_fn) { fprintf(stderr, "[mhaptools convert] use either -r or -b\n"); return 1; }
    if (o.max_frag < 1) { fprintf(stderr, "[mhaptools convert] -L must be >= 1\n"); return 1; }

    clock_t t0 = clock();
    cvt_t c;
    memset(&c, 0, sizeof(c));
    c.o = &o;
    if (!(c.fp = sam_open(o.in_fn, "r"))) { fprintf(stderr, "[mhaptools convert] cannot open %s\n", o.in_fn); return 1; }
    if (o.ref_fn && hts_set_fai_filename(c.fp, o.ref_fn) < 0) { fprintf(stderr, "[mhaptools convert] cannot use reference %s\n", o.ref_fn); return 1; }
    if (o.threads > 0) hts_set_threads(c.fp, o.threads);
    if (!(c.hdr = sam_hdr_read(c.fp))) { fprintf(stderr, "[mhaptools convert] cannot read the header of %s\n", o.in_fn); return 1; }

    hts_idx_t *idx = NULL;
    char **regs = NULL;
    int nregs = 0;
    if (o.region || o.bed_fn) {
        if (!(idx = sam_index_load(c.fp, o.in_fn))) { fprintf(stderr, "[mhaptools convert] -r/-b need an index of %s (samtools index)\n", o.in_fn); return 1; }
        if (!(regs = read_regions(&o, c.hdr, &nregs))) return 1;
        if (nregs && !(c.itr = sam_itr_regarray(idx, c.hdr, regs, nregs))) { fprintf(stderr, "[mhaptools convert] invalid region(s)\n"); return 1; }
    }

    if ((c.tbx = tbx_index_load3(o.cpg_fn, NULL, HTS_IDX_SILENT_FAIL))) {
        if (!(c.cpg_fp = hts_open(o.cpg_fn, "r"))) { fprintf(stderr, "[mhaptools convert] cannot open %s\n", o.cpg_fn); return 1; }
    } else if (load_all_cpgs(&c, o.cpg_fn) < 0) {
        return 1;
    }

    int gz = strlen(o.out_fn) > 3 && !strcmp(o.out_fn + strlen(o.out_fn) - 3, ".gz");
    if (!(c.out = bgzf_open(o.out_fn, gz ? "w" : "wu"))) { fprintf(stderr, "[mhaptools convert] cannot write %s\n", o.out_fn); return 1; }
    if (gz && o.threads > 0) bgzf_mt(c.out, o.threads, 256);

    c.b = bam_init1();
    c.pairer = bs_pairer_init(o.max_frag, frag_done, &c);
    int rc = (regs && !nregs) ? 0 : run(&c);

    if (bgzf_close(c.out) < 0) { fprintf(stderr, "[mhaptools convert] error closing %s\n", o.out_fn); rc = -1; }
    if (rc == 0 && gz && o.index && strcmp(o.out_fn, "-")) {
        if (c.n_lines == 0) {
            fprintf(stderr, "[mhaptools convert] no records; index not written\n");
        } else {
            tbx_conf_t conf = {TBX_GENERIC, 1, 2, 3, '#', 0};
            if (tbx_index_build3(o.out_fn, NULL, 0, o.threads, &conf) < 0) {
                fprintf(stderr, "[mhaptools convert] failed to index %s\n", o.out_fn);
                rc = -1;
            }
        }
    }

    if (c.n_no_cpg_chr)
        fprintf(stderr, "[mhaptools convert] %d contig(s) with reads are absent from the CpG file, their reads were skipped: %s%s\n",
                c.n_no_cpg_chr, c.no_cpg_chr.s, c.n_no_cpg_chr > 5 ? ", ..." : "");
    fprintf(stderr, "[mhaptools convert] %" PRId64 " reads in, %" PRId64 " used (%" PRId64 " filtered, %" PRId64 " incompletely converted), "
            "%" PRId64 " fragments (%" PRId64 " pairs)\n",
            c.n_in, c.n_reads, c.n_filtered, c.n_ch, c.n_frags, c.n_paired);
    fprintf(stderr, "[mhaptools convert] %" PRId64 " fragments without CpG calls, %" PRId64 " dropped for a CpG without a call, "
            "%" PRId64 " cut into several records\n", c.n_nocpg, c.n_dropped, c.n_split);
    fprintf(stderr, "[mhaptools convert] %" PRId64 " records, %" PRId64 " lines written to %s, %.2f s CPU\n",
            c.n_records, c.n_lines, o.out_fn, (double) (clock() - t0) / CLOCKS_PER_SEC);

    bs_pairer_destroy(c.pairer);
    bs_buf_free(&c.buf);
    free(c.rec);
    free(c.ks.s);
    free(c.line.s);
    free(c.no_cpg_chr.s);
    if (c.tbx) { tbx_destroy(c.tbx); hts_close(c.cpg_fp); }
    if (c.cpg_names) {
        for (khint_t k = kh_begin(c.cpg_names); k != kh_end(c.cpg_names); k++)
            if (kh_exist(c.cpg_names, k)) free((char *) kh_key(c.cpg_names, k));
        kh_destroy(cpgchr, c.cpg_names);
        for (int i = 0; i < c.n_cpg_all; i++) bs_cpgs_free(&c.cpg_all[i]);
        free(c.cpg_all);
    }
    for (int i = 0; i < nregs; i++) free(regs[i]);
    free(regs);
    if (c.itr) hts_itr_destroy(c.itr);
    if (idx) hts_idx_destroy(idx);
    bam_destroy1(c.b);
    sam_hdr_destroy(c.hdr);
    sam_close(c.fp);
    return rc ? 1 : 0;
}
