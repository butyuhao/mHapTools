//
// Created by Yuhao Dan on 2020/8/2.
//

#ifndef HAPTOOLS__SUMMARY_H_
#define HAPTOOLS__SUMMARY_H_

#include "hap.h"
#include <string>

namespace std {

typedef struct region_t {
  string chr;
  hap_pos_t beg;
  hap_pos_t end;
}region_t;

typedef struct summary_t {
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
  char *fn_bed;
  char *fn_out;
  hapFile *fp_hap;
  FILE *fp_bed;
  bool stranded;
  char *region;

  vector<string> summary_result;

};
int main_summary(int argc, char *argv[]);
}

#endif //HAPTOOLS__SUMMARY_H_