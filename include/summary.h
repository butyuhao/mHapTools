//
// Created by Yuhao Dan on 2020/8/2.
//

#ifndef HAPTOOLS__SUMMARY_H_
#define HAPTOOLS__SUMMARY_H_

#include <string>
#include <map>
#include "../htslib-1.10.2/htslib/hts.h"
#include "mhap.h"

namespace std {

typedef struct region_t {
  string chr;
  mhap_pos_t beg;
  mhap_pos_t end;
}region_t;

typedef struct summary_t {
  int is_empty() {
    if (n_reads == 0 && m_base == 0 && t_base == 0 && n_reads_k4 == 0 && n_dr == 0 &&
        n_reads_r == 0 && m_base_r == 0 && t_base_r == 0 && n_reads_k4_r == 0 && n_dr_r == 0) {
      return 1;
    }
    return 0;
  }
  int is_direction_plus_empty() {
    if (n_reads == 0 && m_base == 0 && t_base == 0 && n_reads_k4 == 0 && n_dr == 0) {
      return 1;
    }
    return 0;
  }
  int is_direction_minus_empty() {
    if (n_reads_r == 0 && m_base_r == 0 && t_base_r == 0 && n_reads_k4_r == 0 && n_dr_r == 0) {
      return 1;
    }
    return 0;
  }
  // unstranded or direction +
  mhap_pos_t n_reads;
  mhap_pos_t m_base;
  mhap_pos_t t_base;
  mhap_pos_t n_reads_k4;
  mhap_pos_t n_dr;
  mhap_pos_t n_mr;
  // when direction is -
  mhap_pos_t n_reads_r;
  mhap_pos_t m_base_r;
  mhap_pos_t t_base_r;
  mhap_pos_t n_reads_k4_r;
  mhap_pos_t n_dr_r;
  mhap_pos_t n_mr_r;

} summary_t;

class ContextSummary {
 public:
  ContextSummary(): fn_hap(NULL), fn_cpg(NULL), fn_bed(NULL), fn_out(NULL),
                    fp_bed(NULL), hap_idx(NULL), stranded(false), region(NULL),
                    fp_hap_gz(NULL), fp_cpg(NULL), genome_wide(false), region_chr_match(false),
                    region_beg_end_match(false){};
  ~ContextSummary();

  char *fn_hap;
  char *fn_cpg;
  char *fn_bed;
  char *fn_out;
  char *region;
  BGZF* fp_hap_gz;
  htsFile *fp_cpg;
  FILE *fp_bed;
  mhap_idx_t *hap_idx;
  bool stranded;
  bool genome_wide;
  bool region_chr_match;
  bool region_beg_end_match;

  vector<string> summary_result;
  map<string, map<mhap_pos_t, summary_t> > genome_wide_map;
};
int main_summary(int argc, char *argv[]);
}

#endif //HAPTOOLS__SUMMARY_H_
