//
// Created by Yuhao Dan on 2020/8/2.
//
#include <getopt.h>
#include <fstream>
#include "./include/summary.h"
#include "./include/hap.h"
#include "../htslib-1.10.2/htslib/hts.h"
#include "../htslib-1.10.2/htslib/regidx.h"


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

int get_summary_within_region(ContextSummary &ctx_sum, region_t &reg_t, summary_t &sum_t) {
  //get summary result within reg_t and put to sum_t
  ctx_sum.fp_hap = hap_open(ctx_sum.fn_hap, "r");
  if (ctx_sum.fp_hap == NULL) {
    return 1;
  }
  hap_t hap_line_t = hap_t {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};
  while(hap_read(ctx_sum.fp_hap, &hap_line_t) == 0) {
    if (hap_line_t.chr != reg_t.chr) {
      continue;
    }
    hap_pos_t cur_m_base = 0;
    hap_pos_t cur_t_base = 0;
    if ((hap_line_t.chr_end >= reg_t.beg && hap_line_t.chr_end <= reg_t.end) ||
        (hap_line_t.chr_beg >= reg_t.beg && hap_line_t.chr_beg <= reg_t.end) ||
        (hap_line_t.chr_beg <= reg_t.beg && hap_line_t.chr_end >= reg_t.end)) {
      if (ctx_sum.stranded) {
        if (hap_line_t.hap_direction == '+') {
          ++sum_t.n_reads;
          sum_t.t_base += hap_line_t.hap_str.size();
          cur_t_base = hap_line_t.hap_str.size();
          for (auto b : hap_line_t.hap_str) {
            if (b == '1') {
              ++sum_t.m_base;
              ++cur_m_base;
            }
          }
          if (cur_m_base >= 4) {
            ++sum_t.n_reads_k4;
          }
          if (cur_m_base > 0 && cur_t_base != cur_m_base) {
            ++sum_t.n_dr;
          }
        } else if (hap_line_t.hap_direction == '-') {
          ++sum_t.n_reads_r;
          sum_t.t_base_r += hap_line_t.hap_str.size();
          cur_t_base = hap_line_t.hap_str.size();
          for (auto b : hap_line_t.hap_str) {
            if (b == '1') {
              ++sum_t.m_base_r;
              ++cur_m_base;
            }
          }
          if (cur_m_base >= 4) {
            ++sum_t.n_reads_k4_r;
          }
          if (cur_m_base > 0 && cur_t_base != cur_m_base) {
            ++sum_t.n_dr_r;
          }
        } else {
          hts_log_error("Contain unknown direction info.");
          return 1;
        }
      } else {
        ++sum_t.n_reads;
        sum_t.t_base += hap_line_t.hap_str.size();
        cur_t_base = hap_line_t.hap_str.size();
        for (auto b : hap_line_t.hap_str) {
          if (b == '1') {
            ++sum_t.m_base;
            ++cur_m_base;
          }
        }
        if (cur_m_base >= 4) {
          ++sum_t.n_reads_k4;
        }
        if (cur_m_base > 0 && cur_t_base != cur_m_base) {
          ++sum_t.n_dr;
        }
      }
    }
  }
  hap_close(ctx_sum.fp_hap);
  return 0;
}

void get_summary_str(ContextSummary &ctx_sum, region_t &reg_t, summary_t &sum_t) {
  if (ctx_sum.stranded == true) {
    ctx_sum.summary_result.push_back(reg_t.chr + '\t' + to_string(reg_t.beg) +
    '\t' + to_string(reg_t.end) + '\t' + '+' + '\t' + to_string(sum_t.n_reads)
    + '\t' + to_string(sum_t.m_base) + '\t' + to_string(sum_t.t_base) + '\t' +
    to_string(sum_t.n_reads_k4) + '\t' + to_string(sum_t.n_dr));

    ctx_sum.summary_result.push_back(reg_t.chr + '\t' + to_string(reg_t.beg) +
        '\t' + to_string(reg_t.end) + '\t' + '-' + '\t' + to_string(sum_t.n_reads_r)
        + '\t' + to_string(sum_t.m_base_r) + '\t' + to_string(sum_t.t_base_r) + '\t' +
        to_string(sum_t.n_reads_k4_r) + '\t' + to_string(sum_t.n_dr_r));
  } else {
    ctx_sum.summary_result.push_back(reg_t.chr + '\t' + to_string(reg_t.beg) +
        '\t' + to_string(reg_t.end) + '\t' + '*' + '\t' + to_string(sum_t.n_reads)
        + '\t' + to_string(sum_t.m_base) + '\t' + to_string(sum_t.t_base) + '\t' +
        to_string(sum_t.n_reads_k4) + '\t' + to_string(sum_t.n_dr));
  }
}
int get_summary(ContextSummary &ctx_sum) {
  int ret = 0;
  if (ctx_sum.region != NULL) {
    summary_t sum_t = summary_t{0,0,0,0,0,0,0,0,0,0};
    region_t reg_t = region_t {"", 0,0};
    ret = get_region(ctx_sum.region, reg_t);
    if (ret == 1) {
      return 1;
    }
    ret = get_summary_within_region(ctx_sum, reg_t, sum_t);
    if (ret == 1) {
      return 1;
    }
    get_summary_str(ctx_sum, reg_t, sum_t);
  }
  if (ctx_sum.fn_bed != NULL) {
    regidx_t *idx = regidx_init(ctx_sum.fn_bed,NULL,NULL,0,NULL);
    regitr_t *itr = regitr_init(idx);
    while (regitr_loop(itr)) {
      summary_t sum_t = summary_t{0,0,0,0,0,0,0,0,0,0};
      region_t reg_t = region_t {"", 0,0};
      reg_t.chr = itr->seq;
      reg_t.beg =  itr->beg;
      reg_t.end =  itr->end + 1;
      ret = get_summary_within_region(ctx_sum, reg_t, sum_t);
      if (ret == 1) {
        return 1;
      }
      get_summary_str(ctx_sum, reg_t, sum_t);
    }
  }
  return 0;
}

void saving_summary(ContextSummary &ctx_sum) {
  string out_stream_name;
  if (ctx_sum.fn_out != NULL) {
    out_stream_name = ctx_sum.fn_out;
  } else {
    out_stream_name = string(ctx_sum.fn_hap) + "_summary.txt";
  }
  ofstream out_stream(out_stream_name);
  for (auto s : ctx_sum.summary_result) {
    out_stream << s << endl;
  }
  out_stream.close();
}

int main_summary(int argc, char *argv[]) {

  ContextSummary ctx_sum = ContextSummary();

  int long_index;

  static const char *opt_string = "i:o:b:sr:";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "output", optional_argument, NULL, 'o' },
      { "bed", optional_argument, NULL, 'b' },
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

  cout << "Processing..." << endl;
  ret = get_summary(ctx_sum);

  if (ret == 1) {
    hts_log_error("Get summary error.");
  }
  cout << "Saving..." << endl;
  saving_summary(ctx_sum);

  return 0;
}
}//namespace std
