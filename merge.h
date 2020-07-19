#ifndef HAPTOOLS__MERGE_H_
#define HAPTOOLS__MERGE_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <htslib/sam.h>


namespace std {

class ContextMerge {
 public:
  ContextMerge () : fn_cpg(NULL), fn_hap1(NULL), fn_hap2(NULL), fn_out(NULL),
                    fp_cpg(NULL){};
  ~ContextMerge() {};

  char *fn_hap1;
  char *fn_hap2;
  char *fn_cpg;
  char *fn_out;

  htsFile *fp_cpg;
  unordered_map<string, vector<hts_pos_t> > cpg_pos_map;

};

} //namespace std 



#endif //HAPTOOLS__MERGE_H_