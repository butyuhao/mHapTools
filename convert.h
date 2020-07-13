//
// Created by Yuhao Dan on 2020/4/13.
//

#ifndef BAM2HAP__CONVERT_H_
#define BAM2HAP__CONVERT_H_

#include <map>
#include <unordered_map>
#include <vector>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/sam.h>

#define test_mode true

class Context;
struct HT_s;

namespace std {

const uint8_t kbase[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};

class Context {
 public:
  Context () :fp_bam(NULL), fp_cpg(NULL), idx_cpg(NULL), has_idx_cpg(false),
              idx_bam(NULL), cpg_itr(NULL), sam_itr(NULL), hdr_bam(NULL),
              aln(NULL), bam_aux_p(NULL), fn_bam(NULL), output_path(NULL),
              aligner(NULL), fn_bed(NULL), fn_cpg(NULL), region(NULL){};
  ~Context();

  bool parse_region();


  htsFile *fp_bam;
  htsFile *fp_cpg;
  tbx_t *idx_cpg;
  bool has_idx_cpg;
  hts_idx_t *idx_bam;
  hts_itr_t *cpg_itr;
  hts_itr_t *sam_itr;
  bam_hdr_t *hdr_bam;
  bam1_t *aln;
  uint8_t *bam_aux_p;

  // options
  char *fn_bam;   /* -i option */
  char *output_path;  /* -o option */
  char *aligner;      /* -a option */
  char *fn_bed;     /* -b option */
  char *fn_cpg;      /* -c option */
  char *region;       /* -r option */

  vector<hts_pos_t> cpg_pos;
  unordered_map<int, vector<hts_pos_t>> cpg_pos_map;

  //region
  string i_chr;
  int i_tid;
  hts_pos_t i_beg;
  hts_pos_t i_end;

  map<string, int> res_map;

  int region_to_parse;

};

struct HT_s {
  HT_s() {}
  HT_s (char *_h_chr, hts_pos_t _h_start, hts_pos_t _h_end, string &_hap_met, int _count, int8_t _WC)
      : h_chr(_h_chr), h_start(_h_start), h_end(_h_end), hap_met(_hap_met), count(_count), WC(_WC), ht_count(0){}

  void get_WC_symbol();
  string to_str() {
    return string(h_chr) + '\t' + to_string(h_start) + '\t' + to_string(h_end) + '\t' + hap_met + to_string(WC);
  }
  int ht_count;
  char *h_chr;
  hts_pos_t h_start;
  hts_pos_t h_end;
  string hap_met;
  int count;
  char WC_symbol;
  int8_t WC;

};

enum Direction {
  DIRECTION_PLUS = 0,
  DIRECTION_MINUS,
  DIRECTION_UNKNOWN,
};

enum RegionToParse {
  SINGLE_REGION = 0,
  MULTI_REGION,
  WHOLE_FILE,
};

class SamRead {
 public:
  SamRead() {}
  ~SamRead();
  bool init(Context &ctx);
  bool haplo_type();
  bool _get_bismark_std();
  bool _get_XM(Context &ctx);
  bool _get_ZS(Context &ctx);
  bool _get_bismark_QC(Context &ctxainain);

  char *XM_string = NULL;
  char *ZS_string = NULL;
  int read_WC = 0;
  bool QC = true;

  char *read_name = NULL;
  uint16_t flag = 0; //
  int read_map_quality = 0; //mapping quality
  uint32_t *read_cigar = NULL;

  char *seq = NULL;//  vector<char> seq;//the sequence of the reads

  char *read_chr = NULL; //contig name (chromosome)
  hts_pos_t read_start = 0; //left most position of alignment in zero based coordianate (+1)
  hts_pos_t read_end = 0;
  hts_pos_t read_len = 0; //length of the read.
  uint8_t *read_qual = NULL; //quality string

  Context *ctx = NULL;

  string _hap_seq = "";
  vector<int8_t> _hap_qual;
  vector<hts_pos_t> _cpg;
  string _hap_met = "";

  HT_s HT = HT_s();
  HT_s merged_HT = HT_s();

};
} // namespace std

#endif //BAM2HAP__CONVERT_H_
