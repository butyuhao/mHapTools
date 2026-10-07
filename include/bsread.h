/*
 * bsread: read-level logic shared by mhapasm and mhapconvert.
 *
 * Reads are laid out on the reference with their CIGAR, their bisulfite strand
 * is taken from the aligner's tags or the SAM flags, mates are paired into
 * fragments, and CpG methylation is called with the wgbs_tools rules
 * (patter::compareSeqToRef, patter_utils::merge_PE).
 */

#ifndef BSREAD_H
#define BSREAD_H

#include <stdint.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

extern const char *bs_prog;   /* program name used in messages */

typedef struct {
    hts_pos_t beg, end; /* 0-based reference span [beg, end) */
    int bottom;         /* 1 if the read comes from the bottom (OB/CTOB) strand */
    int read1;          /* 1 for read 1 of a pair and for single-end reads */
    char *seq;          /* bases laid out on the reference, 'N' in deletions */
    uint8_t *qual;      /* base qualities laid out on the reference, 0 in deletions */
} mate_t;

typedef struct {
    char *qname;
    hts_pos_t beg, end; /* union of the mate spans */
    hts_pos_t mpos;     /* expected mate position while waiting for the mate */
    int heap_idx;       /* position in the pending-mate heap, -1 if not pending */
    int nmate;
    mate_t m[2];
    int ncall;
    int32_t *cpg;       /* called CpGs, as ascending indexes into the CpG array */
    uint8_t *meth;      /* 1 = methylated, 0 = unmethylated */
} frag_t;

typedef struct {
    hts_pos_t *pos;     /* 0-based positions of the CpG C, ascending */
    uint8_t *mask;      /* 1 = ignore this CpG */
    int n, m;
} cpgs_t;

typedef struct {        /* scratch space for the per-mate CpG calls */
    int32_t *idx[2];
    uint8_t *st[2];
    int m;
} bs_buf_t;

void *bs_malloc(size_t n);
void *bs_realloc(void *p, size_t n);
char *bs_strdup(const char *s);
int bs_lower_bound(const hts_pos_t *v, int n, hts_pos_t x);

/* Contig lookups that accept a missing or extra "chr" prefix and the usual
 * names of the mitochondrial genome. */
int bs_bam_tid(sam_hdr_t *hdr, const char *name);
int bs_tbx_tid(tbx_t *tbx, const char *name);

/* CpGs of [beg, end) from a tabix-indexed file (chr, 1-based C position, ...).
 * Returns the number of CpGs, or -1 if the contig is absent from the file. */
int bs_cpgs_load(htsFile *fp, tbx_t *tbx, const char *chr, hts_pos_t beg, hts_pos_t end,
                 cpgs_t *g, kstring_t *ks);
void bs_cpgs_free(cpgs_t *g);

int bs_read_is_bottom(const bam1_t *b);
int bs_mate_init(mate_t *m, const bam1_t *b);
void bs_mate_free(mate_t *m);
int bs_mate_calls(const mate_t *m, const cpgs_t *g, int min_bq, int32_t *idx, uint8_t *st);
void bs_mate_ch(const mate_t *m, const cpgs_t *g, int *n_c, int *n_t);

frag_t *bs_frag_new(const bam1_t *b, const mate_t *m);
void bs_frag_free(frag_t *f);
void bs_frag_call(frag_t *f, const cpgs_t *g, int min_bq, bs_buf_t *buf);
void bs_buf_free(bs_buf_t *buf);

/* Pairs mates by QNAME, like wgbs_tools match_maker. A read whose mate is
 * unmapped, on another contig or farther than max_frag becomes a fragment on
 * its own; so does a read whose mate has not shown up once the input is past
 * the mate position (see bs_pairer_flush). Finished fragments are handed to
 * `done`, which takes ownership. */
typedef void (*bs_frag_done_f)(frag_t *f, void *data);
typedef struct bs_pairer_s bs_pairer_t;

bs_pairer_t *bs_pairer_init(int max_frag, bs_frag_done_f done, void *data);
void bs_pairer_add(bs_pairer_t *p, const bam1_t *b, const mate_t *m);
void bs_pairer_flush(bs_pairer_t *p, hts_pos_t pos);
void bs_pairer_destroy(bs_pairer_t *p);

#endif
