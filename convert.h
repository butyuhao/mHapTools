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
#define test_mode true

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

  vector<uint32_t> cpg_pos;

  //region
  string i_chr;
  int i_start = -1;
  int i_end = -1;


};

struct HT_s {
  HT_s() {}
  HT_s (char *_h_chr, uint32_t _h_start, uint32_t _h_end, string &_hap_met, int _count, int8_t _WC)
      : h_chr(_h_chr), h_start(_h_start), h_end(_h_end), hap_met(_hap_met), count(_count), WC(_WC){}
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
  void _get_bismark_std();
  bool _get_XM_tag();
  bool _get_ZS_tag();
  void _get_bismark_QC();
  uint8_t* _my_bam_get_seq(context &ctx);
  uint8_t *XM_tag;
  uint8_t *ZS_tag;
  int8_t WC = -1;
  bool QC = false;

  char *qname;
  uint16_t flag; //
  char *rname;
  uint32_t pos;
  uint8_t mapq; //mapping quality
  uint32_t *cigar;
  char *rnext;
  int pnext;
  int32_t tlen;
  vector<char> seq;//the sequence of the reads

  char *chr; //contig name (chromosome)
  int32_t start; //left most position of alignment in zero based coordianate (+1)
  int32_t end;
  uint32_t len; //length of the read.
  uint8_t *qual; //quality string

  context *ctx;

  string _hap_seq;
  vector<int8_t> _hap_qual;
  vector<uint32_t> _cpg;
  string _hap_met;

  HT_s HT;

  enum direction {
    DIRECTION_PLUS = 0,
    DIRECTION_MINUS
  };
};
#endif //BAM2HAP__CONVERT_H_
