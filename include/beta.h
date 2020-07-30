//
// Created by Yuhao Dan on 2020/7/26.
//

#ifndef HAPTOOLS__BETA_H_
#define HAPTOOLS__BETA_H_

#include <map>
#include <vector>
#include <unordered_map>
#include "hap.h"
#include "../htslib-1.10.2/htslib/hts.h"

namespace std {

typedef struct beta_t {
  // unstranded or direction +
  int methy_reads;
  int total_reads;
  // when direction is -
  int methy_reads_r;
  int total_reads_r;
} beta_t;

int main_beta(int argc, char *argv[]);

class ContextBeta {
 public:
  ContextBeta() :fn_hap(NULL), fn_cpg(NULL), fn_bed(NULL), fn_out(NULL),
  fp_hap(NULL), fp_cpg(NULL), fp_bed(NULL), stranded(false) {}
  ~ContextBeta();

  char *fn_hap;
  char *fn_cpg;
  char *fn_bed;
  char *fn_out;
  hapFile *fp_hap;
  htsFile *fp_cpg;
  FILE *fp_bed;
  bool stranded;

  map<string, map<hap_pos_t, beta_t> > beta_map;
  unordered_map<string, vector<hts_pos_t> > cpg_pos_map;


};

}// namespace std

#endif //HAPTOOLS__BETA_H_
