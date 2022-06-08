//
// Created by Yuhao Dan on 2020/7/26.
//

#ifndef HAPTOOLS__BETA_H_
#define HAPTOOLS__BETA_H_

#include <map>
#include <vector>
#include <unordered_map>
#include "mhap.h"
#include "../htslib-1.10.2/htslib/hts.h"

namespace std {

typedef struct beta_t {
  // unstranded or direction +
  int methy_reads;
  int total_reads;
  // when direction is -
  int methy_reads_r;
  int total_reads_r;

  bool is_in_bed;
} beta_t;

int main_beta(int argc, char *argv[]);

class ContextBeta {
 public:
  ContextBeta() :fn_hap(NULL), fn_cpg(NULL), fn_bed(NULL), fn_out(NULL),
  fp_hap_gz(NULL), fp_cpg(NULL), fp_bed(NULL), stranded(false), is_process_beta_checked(false) {}
  ~ContextBeta();

  char *fn_hap;
  char *fn_cpg;
  char *fn_bed;
  char *fn_out;
  BGZF *fp_hap_gz;
  htsFile *fp_cpg;
  FILE *fp_bed;
  bool stranded;
  bool is_process_beta_checked;

  map<string, map<mhap_pos_t, beta_t>, less<string>> beta_map;
  unordered_map<string, vector<hts_pos_t> > cpg_pos_map;
  map<string, beta_t, less<string> > beta_with_bed_results;

};

}// namespace std

#endif //HAPTOOLS__BETA_H_
