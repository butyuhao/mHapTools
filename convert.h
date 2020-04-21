//
// Created by Yuhao Dan on 2020/4/13.
//

#ifndef BAM2HAP__CONVERT_H_
#define BAM2HAP__CONVERT_H_
#include <vector>
#include <string_view>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/kseq.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <vector>

using namespace std;

class context {
 public:
  context () {};
  ~context();
  void init_ctx();
  bool parse_region();
  void print_region();

  htsFile *fp_bam;
  htsFile *fp_cpg;
  bam_hdr_t *hdr_bam;/* -i bam文件头部的指针 */

  bam1_t *aln;

  char *fn_bam;
  char *fn_cpg;
  char *bam_path;   /* -i option */
  char *output_path;  /* -o option */
  char *aligner;      /* -a option */
  char *bed_file;     /* -b option */
  char *cpg_path;      /* -c option */
  char *region;       /* -r option */

  vector<int> cpg_pos;

  //region
  string i_chr;
  int i_start = -1;
  int i_end = -1;


};
#endif //BAM2HAP__CONVERT_H_
