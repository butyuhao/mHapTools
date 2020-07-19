#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include "hap.cpp"
#include "merge.h"


namespace std {

bool load_chr_cpg(ContextMerge &ctx_merge) {
  /*
   * load all the cpg pos to ctx_merge
   */
  ctx_merge.fp_cpg = hts_open(ctx_merge.fn_cpg, "r");
  kstring_t cpg_line = {0,0,NULL};
  unordered_map<string, vector<hts_pos_t>>::iterator cpg_pos_map_itor;
  while (hts_getline(ctx_merge.fp_cpg, KS_SEP_LINE, &cpg_line) > 0) {

    char *p ,*q;
    string chr;
    hts_pos_t cpg_start;
    p = q = cpg_line.s;

    while(*q && *q != '\t') {
      q++;
    }
    *q = '\0';
    chr = string(p);
    *q = '\t';
    p = q + 1;
    q = p;
    while(*q && *q != '\t') {
      q++;
    }
    *q = '\0';
    cpg_start = atoll(p);

    cpg_pos_map_itor = ctx_merge.cpg_pos_map.find(chr);

    if (cpg_pos_map_itor == ctx_merge.cpg_pos_map.end()) {
      vector<hts_pos_t> v;
      v.push_back(cpg_start);
      ctx_merge.cpg_pos_map[chr] = v;
    } else {
      ctx_merge.cpg_pos_map[chr].push_back(cpg_start);
    }
  }

  return true;
}

  bool opt_check(ContextMerge &ctx_merge) {
  if (ctx_merge.fn_hap2 == NULL || ctx_merge.fn_hap1 == NULL || ctx_merge.fn_out == NULL) {
    return false;
  }
  return true;
}

  int main_merge(int argc, char *argv[]) {

  ContextMerge ctx_merge = ContextMerge();

  int long_index;

  static const char *opt_string = "i:o:c:";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "output", optional_argument, NULL, 'o' },
      { "cpg_path", required_argument, NULL, 'c' },
      { NULL, no_argument, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        if (argc < 4) {
          return 1;
        } else {
          ctx_merge.fn_hap1 = optarg;
          ctx_merge.fn_hap2 = argv[optind];
        }
        break;
      }
      case 'c': {
        ctx_merge.fn_cpg = optarg;
      }
      case 'o': {
        ctx_merge.fn_out = optarg;
        break;
      }
      default: {
        break;
      }
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }

  if (!opt_check(ctx_merge)) {
    hts_log_error("opt error");
  }

  load_chr_cpg(ctx_merge);

  hapFile *fp_hap1 = hap_open(ctx_merge.fn_hap1, "rb");
  hapFile *fp_hap2 = hap_open(ctx_merge.fn_hap2, "rb");

  hap_t hap_t_1 = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};

  hap_t hap_t_2 = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};


  while(hap_read(fp_hap1, &hap_t_1) == 0 && hap_read(fp_hap2, &hap_t_2) == 0) {

  }




  return 0;
}
}//namespace std
