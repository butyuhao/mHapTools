#include <unordered_map>
#include <cstring>
#include "mhap.h"
namespace std {

//todo: 将hap的读取与解析改成c++风格而不是c风格
mHapFile* mhap_open(const char *filename, const char *mode) {
  return fopen(filename, mode);
}

void parse_mhap_line(char *line, int buf_len, mhap_t *h_line_t) {
  char *p, *q;
  p = q = line;

  while(*q != MHAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->chr = string(p);
  (*q) = MHAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != MHAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->chr_beg = atoll(p);
  (*q) = MHAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != MHAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->chr_end = atoll(p);
  (*q) = MHAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != MHAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->mhap_str = string(p);
  (*q) = MHAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != MHAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->mhap_count = atoi(p);
  (*q) = '\n';
  p = q + 1;

  h_line_t->mhap_direction = *p;

};

int mhap_read(mHapFile *const fp, mhap_t *h_line_t) {
  // Returns 0 on success,
  //        -1 on EOF,
  static char buf [MHAP_BUF_SIZE];
  fgets(buf, MHAP_BUF_SIZE, fp);
  if (feof(fp)) {
    return -1;
  }
  int buf_len = strlen(buf);

  if (buf_len == 1) {
    return -1;
  }

  parse_mhap_line(buf, buf_len, h_line_t);

  return 0;
}

int mhap_close(mHapFile *fp) {
  if (fp) {
    fclose(fp);
  }
  return 0;
}

mhap_idx_t *mhap_index_load(const char *fn) {
  return tbx_index_load(fn);
}

int mhap_name2id(mhap_idx_t *mhap_idx_t, const char *ss) {
  return tbx_name2id(mhap_idx_t, ss);
}


int mhap_index_build(const char *hap_fn) {
  /*
 * return: 0 (success), -1 (general failure) or -2 (compression not BGZF)
 */
  tbx_conf_t conf = tbx_conf_t{TBX_GENERIC, 1, 2 ,3};
  const tbx_conf_t* conf_p = &conf;

  int ret = tbx_index_build(hap_fn, 0, conf_p);
  return ret;
}

}// namespace std

