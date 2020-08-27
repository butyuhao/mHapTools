//
// Created by Yuhao Dan on 2020/8/2.
//

#ifndef HAPTOOLS__SUMMARY_H_
#define HAPTOOLS__SUMMARY_H_

#include <string>
#include <map>
#include "../htslib-1.10.2/htslib/hts.h"
#include "hap.h"

namespace std {

typedef struct region_t {
  string chr;
  hap_pos_t beg;
  hap_pos_t end;
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
  hap_pos_t n_reads;
  hap_pos_t m_base;
  hap_pos_t t_base;
  hap_pos_t n_reads_k4;
  hap_pos_t n_dr;
  // when direction is -
  hap_pos_t n_reads_r;
  hap_pos_t m_base_r;
  hap_pos_t t_base_r;
  hap_pos_t n_reads_k4_r;
  hap_pos_t n_dr_r;

} summary_t;

class ContextSummary {
 public:
  ContextSummary():fn_hap(NULL), fn_bed(NULL), fn_out(NULL), fp_hap(NULL),
                    fp_bed(NULL), stranded(false), region(NULL){};
  ~ContextSummary();

  char *fn_hap;
  char *fn_cpg;
  char *fn_bed;
  char *fn_out;
  char *region;
  hapFile *fp_hap;
  htsFile *fp_cpg;
  FILE *fp_bed;
  bool stranded;
  bool genome_wide = false;
  bool region_chr_match = false;
  bool region_beg_end_match = false;

  vector<string> summary_result;
  map<string, map<hap_pos_t, summary_t> > genome_wide_map;
};
int main_summary(int argc, char *argv[]);
}

#endif //HAPTOOLS__SUMMARY_H_
