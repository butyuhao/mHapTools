/*
 * bsread: read-level logic shared by mhapasm and mhapconvert.
 *
 * The rules follow wgbs_tools (nloyfer/wgbs_tools):
 *   - read layout:      patter_utils::clean_CIGAR
 *   - strand:           patter_utils::is_bottom, with the aligner's strand tags
 *                       (Bismark XG, bwa-meth/BISCUIT YD, BSMAP ZS) taking precedence
 *   - methylation call: patter::compareSeqToRef (CpG context checked on the read)
 *   - mate merging:     patter_utils::merge_PE (CpGs on which the mates disagree dropped)
 *   - pairing:          match_maker
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/khash.h>
#include "bsread.h"

const char *bs_prog = "bsread";

KHASH_MAP_INIT_STR(bspend, frag_t *)

/***************************************************************
 *                        small helpers                        *
 ***************************************************************/

void *bs_malloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "[%s] out of memory\n", bs_prog); exit(1); }
    return p;
}

void *bs_realloc(void *p, size_t n)
{
    p = realloc(p, n ? n : 1);
    if (!p) { fprintf(stderr, "[%s] out of memory\n", bs_prog); exit(1); }
    return p;
}

char *bs_strdup(const char *s)
{
    size_t n = strlen(s) + 1;
    return memcpy(bs_malloc(n), s, n);
}

/* first index i with v[i] >= x */
int bs_lower_bound(const hts_pos_t *v, int n, hts_pos_t x)
{
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (v[mid] < x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static int name2id_fallback(int (*lookup)(void *, const char *), void *h, const char *name)
{
    static const char *mito[] = {"chrM", "MT", "M", "chrMT"};
    int id = lookup(h, name);
    if (id >= 0) return id;
    for (int i = 0; i < 4; i++) {
        if (strcmp(name, mito[i])) continue;
        for (int j = 0; j < 4; j++)
            if (j != i && (id = lookup(h, mito[j])) >= 0) return id;
        return -1;
    }
    if (!strncmp(name, "chr", 3)) return lookup(h, name + 3);
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "chr%s", name);
    id = lookup(h, ks.s);
    free(ks.s);
    return id;
}

static int bam_lookup(void *h, const char *name) { return sam_hdr_name2tid((sam_hdr_t *) h, name); }
static int tbx_lookup(void *h, const char *name) { return tbx_name2id((tbx_t *) h, name); }

int bs_bam_tid(sam_hdr_t *hdr, const char *name) { return name2id_fallback(bam_lookup, hdr, name); }
int bs_tbx_tid(tbx_t *tbx, const char *name) { return name2id_fallback(tbx_lookup, tbx, name); }

/***************************************************************
 *                        CpG positions                        *
 ***************************************************************/

int bs_cpgs_load(htsFile *fp, tbx_t *tbx, const char *chr, hts_pos_t beg, hts_pos_t end,
                 cpgs_t *g, kstring_t *ks)
{
    g->n = 0;
    int tid = bs_tbx_tid(tbx, chr);
    if (tid < 0) return -1;
    hts_itr_t *itr = tbx_itr_queryi(tbx, tid, beg < 0 ? 0 : beg, end);
    if (!itr) return 0;
    while (tbx_itr_next(fp, tbx, itr, ks) >= 0) {
        char *p = strchr(ks->s, '\t'), *q;
        if (!p) continue;
        long long pos = strtoll(p + 1, &q, 10);
        if (q == p + 1 || pos < 1) continue;
        if (g->n && pos - 1 <= g->pos[g->n - 1]) continue;  /* duplicates */
        if (g->n == g->m) {
            g->m = g->m ? g->m * 2 : 4096;
            g->pos = bs_realloc(g->pos, g->m * sizeof(hts_pos_t));
            g->mask = bs_realloc(g->mask, g->m);
        }
        g->pos[g->n] = pos - 1;
        g->mask[g->n++] = 0;
    }
    tbx_itr_destroy(itr);
    return g->n;
}

void bs_cpgs_free(cpgs_t *g)
{
    free(g->pos);
    free(g->mask);
    g->pos = NULL;
    g->mask = NULL;
    g->n = g->m = 0;
}

/***************************************************************
 *                        reads                                *
 ***************************************************************/

/* Bisulfite strand of a read. Strand tags written by the aligner win
 * (Bismark XG, bwa-meth/BISCUIT YD, BSMAP ZS). Otherwise use the wgbs_tools
 * is_bottom rule (OB = read1 reverse / read2 forward), without requiring the
 * proper-pair bit. */
int bs_read_is_bottom(const bam1_t *b)
{
    uint8_t *t;
    if ((t = bam_aux_get(b, "XG")) && *t == 'Z') {
        const char *v = bam_aux2Z(t);
        if (v && !strcmp(v, "GA")) return 1;
        if (v && !strcmp(v, "CT")) return 0;
    }
    if ((t = bam_aux_get(b, "YD"))) {
        char v = 0;
        if (*t == 'A') v = bam_aux2A(t);
        else if (*t == 'Z') { const char *z = bam_aux2Z(t); v = z ? z[0] : 0; }
        if (v == 'r' || v == 'R') return 1;
        if (v == 'f' || v == 'F') return 0;
    }
    if ((t = bam_aux_get(b, "ZS")) && *t == 'Z') {
        const char *v = bam_aux2Z(t);
        if (v && v[0] == '-') return 1;
        if (v && v[0] == '+') return 0;
    }
    uint16_t fl = b->core.flag;
    if (fl & BAM_FPAIRED) {
        if (fl & BAM_FREAD1) return (fl & BAM_FREVERSE) != 0;
        if (fl & BAM_FREAD2) return (fl & BAM_FREVERSE) == 0;
    }
    return (fl & BAM_FREVERSE) != 0;
}

/* Lay the read out on the reference (patter_utils::clean_CIGAR): keep M/=/X,
 * drop I/S/H/P, fill D/N with 'N'. */
int bs_mate_init(mate_t *m, const bam1_t *b)
{
    const uint32_t *cig = bam_get_cigar(b);
    int ncig = b->core.n_cigar;
    hts_pos_t rlen = bam_cigar2rlen(ncig, cig);
    if (rlen <= 0) return -1;
    m->beg = b->core.pos;
    m->end = b->core.pos + rlen;
    m->bottom = bs_read_is_bottom(b);
    m->read1 = !(b->core.flag & BAM_FPAIRED) || (b->core.flag & BAM_FREAD1);
    m->seq = bs_malloc(rlen);
    m->qual = bs_malloc(rlen);
    const uint8_t *s = bam_get_seq(b), *q = bam_get_qual(b);
    hts_pos_t r = 0;
    int32_t qi = 0;
    for (int i = 0; i < ncig; i++) {
        int op = bam_cigar_op(cig[i]);
        int32_t len = bam_cigar_oplen(cig[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int32_t j = 0; j < len; j++, r++, qi++) {
                m->seq[r] = seq_nt16_str[bam_seqi(s, qi)];
                m->qual[r] = q[qi];
            }
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            memset(m->seq + r, 'N', len);
            memset(m->qual + r, 0, len);
            r += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qi += len;
        }
    }
    return 0;
}

void bs_mate_free(mate_t *m)
{
    free(m->seq);
    free(m->qual);
    m->seq = NULL;
    m->qual = NULL;
}

/* CpG calls of one read (patter::compareSeqToRef). The whole CpG must lie in
 * the aligned span and the read must show the CpG context itself: C/T + G on
 * the top strand, C + G/A on the bottom strand. Only informative calls are
 * returned. */
int bs_mate_calls(const mate_t *m, const cpgs_t *g, int min_bq, int32_t *idx, uint8_t *st)
{
    int n = 0;
    for (int i = bs_lower_bound(g->pos, g->n, m->beg); i < g->n && g->pos[i] + 1 < m->end; i++) {
        if (g->mask[i]) continue;
        hts_pos_t c = g->pos[i] - m->beg;
        char base, ctx;
        uint8_t q;
        if (!m->bottom) { base = m->seq[c]; ctx = m->seq[c + 1]; q = m->qual[c]; }
        else            { base = m->seq[c + 1]; ctx = m->seq[c]; q = m->qual[c + 1]; }
        if (q < min_bq) continue;
        int s = -1;
        if (!m->bottom) { if (ctx == 'G') s = base == 'C' ? 1 : base == 'T' ? 0 : -1; }
        else            { if (ctx == 'C') s = base == 'G' ? 1 : base == 'A' ? 0 : -1; }
        if (s < 0) continue;
        idx[n] = i;
        st[n++] = (uint8_t) s;
    }
    return n;
}

/* Cytosines of a read outside CpGs (reference-free: positions are classified
 * by the CpG list only). *n_c counts unconverted ones (C on the top strand, G
 * on the bottom strand); *n_t counts T (A), i.e. converted cytosines plus the
 * reference's own T (A). An incompletely converted read has a high
 * n_c / (n_c + n_t); a converted read has a few n_c from sequencing errors,
 * non-CpG methylation or bases added in library preparation. */
void bs_mate_ch(const mate_t *m, const cpgs_t *g, int *n_c, int *n_t)
{
    int i = bs_lower_bound(g->pos, g->n, m->beg - 1);
    char cyt = m->bottom ? 'G' : 'C', conv = m->bottom ? 'A' : 'T';
    *n_c = *n_t = 0;
    for (hts_pos_t r = 0; r < m->end - m->beg; r++) {
        char b = m->seq[r];
        if (b != cyt && b != conv) continue;
        hts_pos_t c = m->beg + r - (m->bottom ? 1 : 0);   /* C of the CpG this base would belong to */
        while (i < g->n && g->pos[i] < c) i++;
        if (i < g->n && g->pos[i] == c) continue;
        if (b == cyt) (*n_c)++; else (*n_t)++;
    }
}

/***************************************************************
 *                        fragments                            *
 ***************************************************************/

frag_t *bs_frag_new(const bam1_t *b, const mate_t *m)
{
    frag_t *f = bs_malloc(sizeof(frag_t));
    memset(f, 0, sizeof(frag_t));
    f->qname = bs_strdup(bam_get_qname(b));
    f->m[0] = *m;
    f->nmate = 1;
    f->heap_idx = -1;
    return f;
}

void bs_frag_free(frag_t *f)
{
    for (int i = 0; i < f->nmate; i++) bs_mate_free(&f->m[i]);
    free(f->cpg);
    free(f->meth);
    free(f->qname);
    free(f);
}

static void grow_buffers(bs_buf_t *buf, int need)
{
    if (need <= buf->m) return;
    buf->m = need * 2;
    for (int k = 0; k < 2; k++) {
        buf->idx[k] = bs_realloc(buf->idx[k], buf->m * sizeof(int32_t));
        buf->st[k] = bs_realloc(buf->st[k], buf->m);
    }
}

void bs_buf_free(bs_buf_t *buf)
{
    for (int k = 0; k < 2; k++) { free(buf->idx[k]); free(buf->st[k]); }
    memset(buf, 0, sizeof(*buf));
}

/* Span and CpG calls of a complete fragment. The calls of the mates are
 * merged as in merge_PE: a CpG seen by one mate is kept, a CpG on which the
 * mates disagree is dropped. */
void bs_frag_call(frag_t *f, const cpgs_t *g, int min_bq, bs_buf_t *buf)
{
    f->beg = f->m[0].beg;
    f->end = f->m[0].end;
    int n[2] = {0, 0};
    for (int k = 0; k < f->nmate; k++) {
        if (f->m[k].beg < f->beg) f->beg = f->m[k].beg;
        if (f->m[k].end > f->end) f->end = f->m[k].end;
        grow_buffers(buf, bs_lower_bound(g->pos, g->n, f->m[k].end) - bs_lower_bound(g->pos, g->n, f->m[k].beg) + 1);
        n[k] = bs_mate_calls(&f->m[k], g, min_bq, buf->idx[k], buf->st[k]);
    }
    free(f->cpg);
    free(f->meth);
    f->cpg = bs_malloc((n[0] + n[1]) * sizeof(int32_t));
    f->meth = bs_malloc(n[0] + n[1]);
    int32_t *i0 = buf->idx[0], *i1 = buf->idx[1];
    uint8_t *s0 = buf->st[0], *s1 = buf->st[1];
    int i = 0, j = 0, k = 0;
    while (i < n[0] || j < n[1]) {
        if (j >= n[1] || (i < n[0] && i0[i] < i1[j])) {
            f->cpg[k] = i0[i]; f->meth[k++] = s0[i++];
        } else if (i >= n[0] || i1[j] < i0[i]) {
            f->cpg[k] = i1[j]; f->meth[k++] = s1[j++];
        } else {
            if (s0[i] == s1[j]) { f->cpg[k] = i0[i]; f->meth[k++] = s0[i]; }
            i++; j++;
        }
    }
    f->ncall = k;
}

/***************************************************************
 *                        mate pairing                         *
 ***************************************************************/

struct bs_pairer_s {
    khash_t(bspend) *pend;  /* reads waiting for their mate, by QNAME */
    frag_t **heap;          /* the same reads as a min-heap on the expected mate position */
    int nheap, mheap;
    int max_frag;
    bs_frag_done_f done;
    void *data;
};

static void heap_swap(bs_pairer_t *p, int i, int j)
{
    frag_t *t = p->heap[i];
    p->heap[i] = p->heap[j];
    p->heap[j] = t;
    p->heap[i]->heap_idx = i;
    p->heap[j]->heap_idx = j;
}

static void heap_up(bs_pairer_t *p, int i)
{
    while (i > 0 && p->heap[(i - 1) / 2]->mpos > p->heap[i]->mpos) {
        heap_swap(p, i, (i - 1) / 2);
        i = (i - 1) / 2;
    }
}

static void heap_down(bs_pairer_t *p, int i)
{
    for (;;) {
        int l = 2 * i + 1, r = l + 1, m = i;
        if (l < p->nheap && p->heap[l]->mpos < p->heap[m]->mpos) m = l;
        if (r < p->nheap && p->heap[r]->mpos < p->heap[m]->mpos) m = r;
        if (m == i) return;
        heap_swap(p, i, m);
        i = m;
    }
}

static void heap_push(bs_pairer_t *p, frag_t *f)
{
    if (p->nheap == p->mheap) {
        p->mheap = p->mheap ? p->mheap * 2 : 256;
        p->heap = bs_realloc(p->heap, p->mheap * sizeof(frag_t *));
    }
    p->heap[p->nheap] = f;
    f->heap_idx = p->nheap++;
    heap_up(p, f->heap_idx);
}

static void heap_remove(bs_pairer_t *p, frag_t *f)
{
    int i = f->heap_idx;
    if (i != --p->nheap) {
        p->heap[i] = p->heap[p->nheap];
        p->heap[i]->heap_idx = i;
        heap_down(p, i);
        heap_up(p, i);
    }
    f->heap_idx = -1;
}

bs_pairer_t *bs_pairer_init(int max_frag, bs_frag_done_f done, void *data)
{
    bs_pairer_t *p = bs_malloc(sizeof(bs_pairer_t));
    memset(p, 0, sizeof(*p));
    p->pend = kh_init(bspend);
    p->max_frag = max_frag;
    p->done = done;
    p->data = data;
    return p;
}

void bs_pairer_add(bs_pairer_t *p, const bam1_t *b, const mate_t *m)
{
    const bam1_core_t *co = &b->core;
    int expect = (co->flag & BAM_FPAIRED) && !(co->flag & BAM_FMUNMAP) && co->mtid == co->tid
                 && llabs((long long) (co->mpos - co->pos)) <= p->max_frag;
    if (expect) {
        khint_t k = kh_get(bspend, p->pend, bam_get_qname(b));
        if (k != kh_end(p->pend)) {
            frag_t *f = kh_val(p->pend, k);
            kh_del(bspend, p->pend, k);
            heap_remove(p, f);
            f->m[f->nmate++] = *m;
            p->done(f, p->data);
            return;
        }
        if (co->mpos >= co->pos) {   /* mate still to come */
            frag_t *f = bs_frag_new(b, m);
            f->mpos = co->mpos;
            int ret;
            k = kh_put(bspend, p->pend, f->qname, &ret);
            if (ret == 0) {          /* same QNAME already waiting: keep both as singles */
                frag_t *old = kh_val(p->pend, k);
                kh_key(p->pend, k) = f->qname;
                heap_remove(p, old);
                p->done(old, p->data);
            }
            kh_val(p->pend, k) = f;
            heap_push(p, f);
            return;
        }
    }
    p->done(bs_frag_new(b, m), p->data);
}

void bs_pairer_flush(bs_pairer_t *p, hts_pos_t pos)
{
    while (p->nheap && p->heap[0]->mpos < pos) {
        frag_t *f = p->heap[0];
        heap_remove(p, f);
        khint_t k = kh_get(bspend, p->pend, f->qname);
        if (k != kh_end(p->pend)) kh_del(bspend, p->pend, k);
        p->done(f, p->data);
    }
}

/* Frees reads still waiting for their mate without handing them out; flush
 * first to keep them. */
void bs_pairer_destroy(bs_pairer_t *p)
{
    if (!p) return;
    for (int i = 0; i < p->nheap; i++) bs_frag_free(p->heap[i]);
    kh_destroy(bspend, p->pend);
    free(p->heap);
    free(p);
}
