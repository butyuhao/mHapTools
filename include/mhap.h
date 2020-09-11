//
// Created by Yuhao Dan on 2020/7/9.
//

#ifndef HAPTOOLS__HAP_H_
#define HAPTOOLS__HAP_H_

#include "../htslib-1.10.2/htslib/tbx.h"
#include "../htslib-1.10.2/htslib/hts.h"
#include "../htslib-1.10.2/htslib/bgzf.h"
#include <vector>
#include <string>
#include <iostream>
#define hapFile FILE
#define HAP_BUF_SIZE 2048
#define HAP_LINE_DELIM_TAB '\t'
#define HAP_DEFAULT_DIRECTION '-'
#define HAP_NULL_STRING ""
#define hap_itr_queryi tbx_itr_queryi
#define hap_idx_t tbx_t
#define hap_tid tbx_tid
#define hap_idx_destroy tbx_destroy
#define hap_itr_t hts_itr_t
#define hap_itr_destroy hts_itr_destroy




typedef int64_t hap_pos_t;

namespace std {

typedef struct hap_t{
  string chr;
  hap_pos_t chr_beg;
  hap_pos_t chr_end;
  string hap_str;
  int hap_count;
  char hap_direction;
  vector<hap_pos_t> cpg_pos;//存储cpg位置;

  void print() {
    cout << chr << '\t' << chr_beg << '\t' << chr_end << '\t' << hap_str
         << '\t' << hap_count << '\t' << hap_direction << endl;
  }

  bool is_valid() {
    if (hap_direction == '+' || hap_direction == '-') {
      return true;
    } else {
      return false;
    }
  }

} hap_t;

hapFile* hap_open(const char *filename, const char *mode);
// Returns 0 on success,
//        -1 on EOF,
int hap_read(hapFile *const fp, hap_t *h_line_t);
int hap_close(hapFile *fp);

hap_idx_t *hap_index_load(const char *fn);

int hap_name2id(hap_idx_t *hap_idx_t, const char *ss);

/// Create a .tbi index for hap file
/** @param hap_fn     Name of the hap file
    @return 0 (success), -1 (general failure) or -2 (compression not BGZF)

    the index name is automatically derived from @p fn.

    suppose the hap file is named XXX.hap.gz,
    the index file is named XXX.hap.gz.tbi
*/
int hap_index_build(const char *hap_fn);

} //namespace std
#endif //HAPTOOLS__HAP_H_

