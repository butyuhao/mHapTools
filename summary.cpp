//
// Created by Yuhao Dan on 2020/8/2.
//
#include <getopt.h>
#include "./include/summary.h"
#include "./include/hap.h"
#include "../htslib-1.10.2/htslib/hts.h"

namespace std {

ContextSummary::~ContextSummary() {
  if (fp_hap) {
    hap_close(fp_hap);
  }
  if (fp_bed) {
    fclose(fp_bed);
  }
}

int summary_opt_check(ContextSummary &ctx_sum) {
  if (ctx_sum.fn_bed == NULL && ctx_sum.region == NULL) {
    return 1;
  }
  return 0;
}

int open_hap(ContextSummary &ctx_sum) {
  ctx_sum.fp_hap = hap_open(ctx_sum.fn_hap, "r");
  if (ctx_sum.fp_hap == NULL) {
    return 1;
  }
  return 0;
}

int get_region(char* region_string, region_t &reg_t) {
  if (region_string == NULL) {
    return 1;
  }
  char *p, *q;
  p  = region_string;
  q = region_string;
  while(*q != '\0' && *q != ':') {
    q++;
  }
  if (*q != '\0') {
    *q = '\0';
  } else {
    return 1;
  }
  reg_t.chr = string(p);
  *q = ':';
  p = q + 1;
  if (*p != '\0') {
    while(*q != '\0' && *q != '-') {
      q++;
    }
    if (*q != '\0') {
      *q = '\0';
    }
    reg_t.beg = atoll(p);
    *q = '-';
    p = q + 1;
    if (*p != '\0') {
      reg_t.end = atoll(p);
    } else {
      return 1;
    }
  } else {
    return 1;
  }
  return 0;
}

int get_summary(ContextSummary &ctx_sum) {
  summary_t sum_t = summary_t{0,0,0,0,0,0,0,0,0,0};
  if (ctx_sum.region != NULL) {

  }
  return 0;
}

int get_summary_within_region(ContextSummary &ctx_sum, summary_t &sum_t) {
  //get summary result within ctx_sum.region
  return 0;
}

int main_summary(int argc, char *argv[]) {

  ContextSummary ctx_sum = ContextSummary();

  int long_index;

  static const char *opt_string = "i:o:b:sr:";

  static const struct option long_opts[] = {
      { "fn_in", required_argument, NULL, 'i' },
      { "fn_out", optional_argument, NULL, 'o' },
      { "fn_bed", optional_argument, NULL, 'b' },
      { "stranded", optional_argument, NULL, 's' },
      { "region", optional_argument, NULL, 'r' },
      { NULL, no_argument, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        ctx_sum.fn_hap = optarg;
        break;
      }
      case 'o': {
        ctx_sum.fn_out = optarg;
        break;
      }
      case 'b': {
        ctx_sum.fn_bed = optarg;
        break;
      }
      case 's': {
        ctx_sum.stranded = true;
        break;
      }
      case 'r': {
        ctx_sum.region = optarg;
        break;
      }
      default: {
        break;
      }
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }

  int ret = 0;

  ret = summary_opt_check(ctx_sum);
  if (ret == 1){
    hts_log_error("opt error");
    return 1;
  }

  ret = open_hap(ctx_sum);
  if (ret == 1) {
    hts_log_error("open hap error");
    return 1;
  }

//  region_t reg_t = region_t {"", 0,0};
//
//  get_region(ctx_sum.region, reg_t);

  ret = get_summary(ctx_sum);



  return 0;
}
}//namespace std
