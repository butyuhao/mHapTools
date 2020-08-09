#include <unordered_map>
#include <cstring>
#include "./include/hap.h"
namespace std {

hapFile* hap_open(const char *filename, const char *mode) {
  return fopen(filename, mode);
}

void parse_hap_line(char *line, int buf_len, hap_t *h_line_t) {
  char *p, *q;
  p = q = line;

  while(*q != HAP_LINE_DELIM_TAB && (q - line) < buf_len - 1) {q++;}
  (*q) = '\0';
  h_line_t->chr = string(p);
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && (q - line) < buf_len - 1) {q++;}
  (*q) = '\0';
  h_line_t->chr_beg = atoll(p);
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && (q - line) < buf_len - 1) {q++;}
  (*q) = '\0';
  h_line_t->chr_end = atoll(p);
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && (q - line) < buf_len - 1) {q++;}
  (*q) = '\0';
  h_line_t->hap_str = string(p);
  q = p = q + 1;

  while(*q != HAP_LINE_DELIM_TAB && (q - line) < buf_len - 1) {q++;}
  (*q) = '\0';
  h_line_t->hap_count = atoi(p);
  p = q + 1;

  h_line_t->hap_direction = *p;

};

int hap_read(hapFile *const fp, hap_t *h_line_t) {
  // Returns 0 on success,
  //        -1 on EOF,
  static char buf [HAP_BUF_SIZE];
  if (!feof(fp)) {
    fgets(buf, HAP_BUF_SIZE, fp);
  } else {
    return -1;
  }
  int buf_len = strlen(buf);

  if (buf_len == 1) {
    return -1;
  }

  if (buf[buf_len] == '\n') {
    buf[buf_len - 1] = '\0';
  } else {
    //in case the last line does not have \n
    buf[buf_len] = '\0';
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


}// namespace std

