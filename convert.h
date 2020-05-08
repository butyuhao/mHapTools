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
#include <map>
#include<fstream>

#define test_mode false

using namespace std;


uint8_t _base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};

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
  uint8_t *bam_aux_p;//bam文件辅助信息的指针

  char *fn_bam;
  char *fn_cpg;
  char *bam_path;   /* -i option */
  char *output_path;  /* -o option */
  char *aligner;      /* -a option */
  char *bed_file;     /* -b option */
  char *cpg_path;      /* -c option */
  char *region;       /* -r option */

  vector<uint32_t> cpg_pos;

  //region
  string i_chr;
  uint32_t i_start = -1;
  uint32_t i_end = -1;

  map<string, u_int32_t > res_map_sort;
  vector<pair<string, u_int32_t>> vt;

};

struct HT_s {
  HT_s() {}
  HT_s (char *_h_chr, uint32_t _h_start, uint32_t _h_end, string &_hap_met, int _count, int8_t _WC)
      : h_chr(_h_chr), h_start(_h_start), h_end(_h_end), hap_met(_hap_met), count(_count), WC(_WC){}
  string to_str() {
    return string(h_chr) + '\t' + to_string(h_start) + '\t' + to_string(h_end) + '\t' + hap_met + to_string(WC);
  }
  u_int32_t get_h_start() {
    return h_start;
  }
  char *h_chr;
  uint32_t h_start;
  uint32_t h_end;
  string hap_met;
  int count;
  int8_t WC;
};

class sam_read {
 public:
  sam_read() {}
  bool init(context &ctx);
  bool haplo_type();
  bool _get_bismark_std();
  bool _get_XM_tag(context &ctx);
  bool _get_ZS_tag(context &ctx);
  bool _get_bismark_QC(context &ctxainain);

  char *XM_tag = NULL;
  char *ZS_tag = NULL;
  int read_WC = 0;
  bool QC = true;

  char *read_name = NULL;
  uint16_t flag = 0; //
  int read_map_quality = 0; //mapping quality
  uint32_t *read_cigar = NULL;

  vector<char> seq;//the sequence of the reads

  char *read_chr = NULL; //contig name (chromosome)
  uint32_t read_start = 0; //left most position of alignment in zero based coordianate (+1)
  uint32_t read_end = 0;
  uint32_t read_len = 0; //length of the read.
  uint8_t *read_qual = NULL; //quality string

  context *ctx = NULL;

  string _hap_seq = "";
  vector<int8_t> _hap_qual;
  vector<uint32_t> _cpg;
  string _hap_met = "";

  HT_s HT = HT_s();
  HT_s merged_HT = HT_s();

  enum direction {
    DIRECTION_PLUS = 0,
    DIRECTION_MINUS
  };
};
#endif //BAM2HAP__CONVERT_H_
