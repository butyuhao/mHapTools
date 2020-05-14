#include <iostream>
#include <string.h>
#include "convert.cpp"
#include <htslib/regidx.h>
using namespace std;

int main(int argc, char *argv[]) {

//  string cpg_idx_fn = string("/Users/butyuhao/Documents/GitHub/haptools/data/hg19_CpG.gz.tbi");
//  string cpg_path = string("/Users/butyuhao/Documents/GitHub/haptools/data/hg19_CpG.gz");
//  htsFile *fp = hts_open(cpg_path.c_str(), "r");
//  //get the tbi index
//  tbx_t *idx_cpg = tbx_index_load(cpg_idx_fn.c_str());
//  int tid = tbx_name2id(idx_cpg, "1");
//  hts_itr_t *cpg_itr = tbx_itr_queryi(idx_cpg, tid, 2000, 20000);
//  kstring_t ksbuf = {0, 0, NULL};
//  while(tbx_itr_next(fp, idx_cpg, cpg_itr, &ksbuf) >= 0) {
//    cout << ksbuf.s << endl;
//  }
//
//  while(1);

  if (argc < 2) {
    cout << "See usage" << endl;
  }
  int ret = 0;
  if (strcmp(argv[1], "convert") == 0)    ret = main_convert(argc - 1, argv + 1);

  return ret;


}