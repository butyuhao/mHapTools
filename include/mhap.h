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
#define mHapFile FILE
#define MHAP_BUF_SIZE 2048
#define MHAP_LINE_DELIM_TAB '\t'
#define MHAP_DEFAULT_DIRECTION '-'
#define MHAP_NULL_STRING ""
#define mhap_itr_queryi tbx_itr_queryi
#define mhap_idx_t tbx_t
#define mhap_tid tbx_tid
#define mhap_idx_destroy tbx_destroy
#define mhap_itr_t hts_itr_t
#define mhap_itr_destroy hts_itr_destroy




typedef int64_t mhap_pos_t;

namespace std {

typedef struct mhap_t{
  string chr;
  mhap_pos_t chr_beg;
  mhap_pos_t chr_end;
  string mhap_str;
  int mhap_count;
  char mhap_direction;
  vector<mhap_pos_t> cpg_pos;//存储cpg位置;

  void print() {
    cout << chr << '\t' << chr_beg << '\t' << chr_end << '\t' << mhap_str
         << '\t' << mhap_count << '\t' << mhap_direction << endl;
  }

  bool is_valid() {
    if (mhap_direction == '+' || mhap_direction == '-') {
      return true;
    } else {
      return false;
    }
  }

} mhap_t;

mHapFile* mhap_open(const char *filename, const char *mode);

void parse_mhap_line(char *line, int buf_len, mhap_t *h_line_t);

// Returns 0 on success,
//        -1 on EOF,
int mhap_read(mHapFile *const fp, mhap_t *h_line_t);
int mhap_close(mHapFile *fp);

mhap_idx_t *mhap_index_load(const char *fn);

int mhap_name2id(mhap_idx_t *mhap_idx_t, const char *ss);

/// Create a .tbi index for hap file
/** @param hap_fn     Name of the hap file
    @return 0 (success), -1 (general failure) or -2 (compression not BGZF)

    the index name is automatically derived from @p fn.

    suppose the hap file is named XXX.hap.gz,
    the index file is named XXX.hap.gz.tbi
*/
int mhap_index_build(const char *hap_fn);

} //namespace std
#endif //HAPTOOLS__HAP_H_

