//
// Created by Yuhao Dan on 2020/7/9.
//

#ifndef HAPTOOLS__HAP_H_
#define HAPTOOLS__HAP_H_

#define hapFile FILE
#define HAP_BUF_SIZE 2048
#define HAP_LINE_DELIM_TAB '\t'
#define HAP_DEFAULT_DIRECTION '-'
#define HAP_NULL_STRING ""
#include <vector>
#include <string>
#include <iostream>
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
}


#endif //HAPTOOLS__HAP_H_

