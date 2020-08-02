//
// Created by Yuhao Dan on 2020/4/13.
//

#ifndef BAM2HAP__CONVERT_H_
#define BAM2HAP__CONVERT_H_

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include "../htslib-1.10.2/htslib/hts.h"
#include "../htslib-1.10.2/htslib/tbx.h"
#include "../htslib-1.10.2/htslib/sam.h"

class Context;
struct HT_s;

namespace std {

const uint8_t kbase[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};

class Context {
 public:
  Context () :fp_bam(NULL), fp_cpg(NULL), idx_cpg(NULL), has_idx_cpg(false),
              idx_bam(NULL), cpg_itr(NULL), sam_itr(NULL), hdr_bam(NULL),
              aln(NULL), bam_aux_p(NULL), fn_bam(NULL), output_path(NULL),
              fn_bed(NULL), fn_cpg(NULL), region(NULL),stranded(false),
              non_directional(false){};
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
  char *fn_bed;     /* -b option */
  char *fn_cpg;      /* -c option */
  char *region;       /* -r option */
  bool stranded;
  bool non_directional;

  vector<hts_pos_t> cpg_pos;
  unordered_map<int, vector<hts_pos_t> > cpg_pos_map;

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
  DIRECTION_UNKNOWN = 0,
  DIRECTION_PLUS,
  DIRECTION_MINUS,
};

enum RegionToParse {
  SINGLE_REGION = 0,
  MULTI_REGION,
  WHOLE_FILE,
};

class SamRead {
 public:
  SamRead() : XM_string(NULL), ZS_string(NULL), read_WC(0), QC (true), read_name(NULL),
              flag(0), read_map_quality(0), read_cigar(NULL), seq(NULL), read_chr(NULL),
              read_start(0), read_end(0), read_len(0), read_qual(NULL), ctx(NULL),
              _hap_seq(""), _hap_met(""){}
  ~SamRead();
  bool init(Context &ctx);
  bool haplo_type();
  bool _get_XM(Context &ctx);
  bool _get_bismark_QC(Context &ctxainain);

  char *XM_string;
  char *ZS_string;
  int read_WC;
  bool QC;

  char *read_name;
  uint16_t flag; //
  int read_map_quality; //mapping quality
  uint32_t *read_cigar;

  char *seq;//  vector<char> seq;//the sequence of the reads

  char *read_chr; //contig name (chromosome)
  hts_pos_t read_start; //left most position of alignment in zero based coordianate (+1)
  hts_pos_t read_end;
  hts_pos_t read_len; //length of the read.
  uint8_t *read_qual; //quality string

  Context *ctx;

  string _hap_seq;
  vector<int8_t> _hap_qual;
  vector<hts_pos_t> _cpg;
  string _hap_met;

  HT_s HT = HT_s();
  HT_s merged_HT = HT_s();

};

int main_convert(int argc, char *argv[]);

} // namespace std

#endif //BAM2HAP__CONVERT_H_
