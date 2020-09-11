#include <unordered_map>
#include <cstring>
#include "mhap.h"
namespace std {

//todo: 将hap的读取与解析改成c++风格而不是c风格
hapFile* hap_open(const char *filename, const char *mode) {
  return fopen(filename, mode);
}

void parse_hap_line(char *line, int buf_len, hap_t *h_line_t) {
  char *p, *q;
  p = q = line;

  while(*q != HAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->chr = string(p);
  (*q) = HAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->chr_beg = atoll(p);
  (*q) = HAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->chr_end = atoll(p);
  (*q) = HAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->hap_str = string(p);
  (*q) = HAP_LINE_DELIM_TAB;
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && *q !='\n' && (q - line) < buf_len) {q++;}
  (*q) = '\0';
  h_line_t->hap_count = atoi(p);
  (*q) = '\n';
  p = q + 1;

  h_line_t->hap_direction = *p;

};

int hap_read(hapFile *const fp, hap_t *h_line_t) {
  // Returns 0 on success,
  //        -1 on EOF,
  static char buf [HAP_BUF_SIZE];
  fgets(buf, HAP_BUF_SIZE, fp);
  if (feof(fp)) {
    return -1;
  }
  int buf_len = strlen(buf);

  if (buf_len == 1) {
    return -1;
  }

  parse_hap_line(buf, buf_len, h_line_t);

  return 0;
}

int hap_close(hapFile *fp) {
  if (fp) {
    fclose(fp);
  }
  return 0;
}

hap_idx_t *hap_index_load(const char *fn) {
  return tbx_index_load(fn);
}

int hap_name2id(hap_idx_t *hap_idx_t, const char *ss) {
  return tbx_name2id(hap_idx_t, ss);
}


int hap_index_build(const char *hap_fn) {
  /*
 * return: 0 (success), -1 (general failure) or -2 (compression not BGZF)
 */
  tbx_conf_t conf = tbx_conf_t{TBX_GENERIC, 1, 2 ,3};
  const tbx_conf_t* conf_p = &conf;

  int ret = tbx_index_build(hap_fn, 0, conf_p);
  return ret;
}

}// namespace std

