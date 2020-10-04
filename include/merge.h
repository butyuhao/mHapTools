#ifndef HAPTOOLS__MERGE_H_
#define HAPTOOLS__MERGE_H_

#include <string>
#include <vector>
#include <unordered_map>
#include "mhap.h"
#include "../htslib-1.10.2/htslib/sam.h"


namespace std {

class ContextMerge {
 public:
  ContextMerge () : fn_cpg(NULL), fn_hap1(NULL), fn_hap2(NULL), fn_out(NULL),
                    fp_hap1_gz(NULL), fp_hap2_gz(NULL), fp_cpg(NULL){};
  ~ContextMerge();

  char *fn_hap1;
  char *fn_hap2;
  char *fn_cpg;
  char *fn_out;
  BGZF *fp_hap1_gz;
  BGZF *fp_hap2_gz;

  htsFile *fp_cpg;
  unordered_map<string, vector<hts_pos_t> > cpg_pos_map;

};

int main_merge(int argc, char *argv[]);

} //namespace std 



#endif //HAPTOOLS__MERGE_H_