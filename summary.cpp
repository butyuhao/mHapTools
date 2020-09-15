//
// Created by Yuhao Dan on 2020/8/2.
//
#include <getopt.h>
#include <fstream>
#include "./include/summary.h"
#include "./include/mhap.h"
#include "./include/utils.h"
#include "./htslib-1.10.2/htslib/regidx.h"
#include "./htslib-1.10.2/htslib/kseq.h"

namespace std {

ContextSummary::~ContextSummary() {
  if (fp_hap) {
    mhap_close(fp_hap);
  }
  if (fp_bed) {
    fclose(fp_bed);
  }
  if (hap_idx) {
    mhap_idx_destroy(hap_idx);
  }
  if (fp_hap) {
    mhap_close(fp_hap);
  }
  if (fp_cpg) {
    hts_close(fp_cpg);
  }
}

int summary_opt_check(ContextSummary &ctx_sum) {
  if (ctx_sum.fn_hap == NULL) {
    hts_log_error("Please specify input file");
    return 1;
  }
  if (ctx_sum.genome_wide) {
    if (ctx_sum.fn_cpg == NULL) {
      hts_log_error("If -g is specified, -c is also required");
      return 1;
    }
  } else {
    if (ctx_sum.fn_bed == NULL && ctx_sum.region == NULL) {
      hts_log_error("Please specify -b or -r");
      return 1;
    }
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
  ctx_sum.fp_hap_gz = bgzf_open(ctx_sum.fn_hap, "r");
  if (ctx_sum.fp_hap_gz == NULL) {
    hts_log_error("Fail to open the mhap file.");
    return 1;
  }
  mhap_t hap_line_t = mhap_t {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};

  int hap_tid = mhap_name2id(ctx_sum.hap_idx, reg_t.chr.c_str());
  if (reg_t.beg <= 500) {
    reg_t.beg = 0;
  } else {
    reg_t.beg -= 500;
  }
  reg_t.end += 500;
  hts_itr_t *hap_itr = mhap_itr_queryi(ctx_sum.hap_idx, hap_tid, reg_t.beg, reg_t.end);

  kstring_t ksbuf = {0, 0, NULL};

  while(hts_itr_next(ctx_sum.fp_hap_gz, hap_itr, &ksbuf, ctx_sum.hap_idx) >= 0) {
    parse_mhap_line(ksbuf.s, ksbuf.l, &hap_line_t);
    ctx_sum.region_chr_match = true;
    mhap_pos_t cur_m_base = 0;
    mhap_pos_t cur_t_base = 0;
    if ((hap_line_t.chr_end >= reg_t.beg && hap_line_t.chr_end <= reg_t.end) ||
        (hap_line_t.chr_beg >= reg_t.beg && hap_line_t.chr_beg <= reg_t.end) ||
        (hap_line_t.chr_beg <= reg_t.beg && hap_line_t.chr_end >= reg_t.end)) {
      ctx_sum.region_beg_end_match = true;
      if (ctx_sum.stranded) {
        if (hap_line_t.mhap_direction == '+') {
          ++sum_t.n_reads;
          sum_t.t_base += hap_line_t.mhap_str.size();
          cur_t_base = hap_line_t.mhap_str.size();
          for (auto b : hap_line_t.mhap_str) {
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
        } else if (hap_line_t.mhap_direction == '-') {
          ++sum_t.n_reads_r;
          sum_t.t_base_r += hap_line_t.mhap_str.size();
          cur_t_base = hap_line_t.mhap_str.size();
          for (auto b : hap_line_t.mhap_str) {
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
        sum_t.t_base += hap_line_t.mhap_str.size();
        cur_t_base = hap_line_t.mhap_str.size();
        for (auto b : hap_line_t.mhap_str) {
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

  if (hap_itr) {
    hts_itr_destroy(hap_itr);
  }
  if (ctx_sum.fp_hap_gz) {
    bgzf_close(ctx_sum.fp_hap_gz);
  }
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

  string filename_hap_idx = string(ctx_sum.fn_hap) + ".tbi";
  ctx_sum.hap_idx = mhap_index_load(filename_hap_idx.c_str());
  if (ctx_sum.hap_idx == NULL) {
    hts_log_error("Fail to open the index file of the input mHap file. Please generate .tbi index file for the mHap file.");
    return 1;
  }

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
    if (idx) {
      regidx_destroy(idx);
    }
    if (itr) {
      regitr_destroy(itr);
    }
  }
  return 0;
}

void saving_summary(ContextSummary &ctx_sum) {
  if (!ctx_sum.region_chr_match || !ctx_sum.region_beg_end_match) {
    if (!ctx_sum.region_chr_match) {
      hts_log_warning("Warning: Check the region you specified (especially the chr name), no reads in the mhap file match the region.");
    }
    if (!ctx_sum.region_beg_end_match) {
      hts_log_warning("Warning: Check the region you specified, no reads in the mhap file match the region.");
    }
  }
  string out_stream_name;
  if (ctx_sum.fn_out != NULL) {
    out_stream_name = ctx_sum.fn_out;
  } else {
    out_stream_name = "summary.txt";
  }
  ofstream out_stream(out_stream_name);
  for (auto s : ctx_sum.summary_result) {
    out_stream << s << endl;
  }
  out_stream.close();
}

inline int parse_cpg_line(char *cpg_line, string *chr, hts_pos_t *beg, hts_pos_t *end) {
  if (cpg_line == NULL || chr == NULL || beg == NULL || end == NULL) {
    return 1;
  }
  char *p ,*q;
  p = q = cpg_line;

  while(*q && *q != '\t') {
    q++;
  }
  *q = '\0';
  *chr = string(p);
  *q = '\t';
  p = q + 1;
  q = p;
  while(*q && *q != '\t') {
    q++;
  }
  *q = '\0';
  *beg = atoll(p);
  *q = '\t';
  p = q + 1;
  *end = atoll(p);
  return 0;
}

int load_cpg_init_map(ContextSummary &ctx_sum) {
  int ret = 0;
  map<string, map<mhap_pos_t, summary_t> >::iterator chr_itor;
  map<mhap_pos_t, summary_t>::iterator cpg_itor;
  ctx_sum.fp_cpg = hts_open(ctx_sum.fn_cpg, "r");
  if (ctx_sum.fp_cpg == NULL) {
    hts_log_error("Fail to open CpG file.");
    return 1;
  }
  kstring_t cpg_line = {0,0,NULL};

  while (hts_getline(ctx_sum.fp_cpg, KS_SEP_LINE, &cpg_line) > 0) {
    string chr;
    hts_pos_t beg;
    hts_pos_t end;
    ret = parse_cpg_line(cpg_line.s, &chr, &beg, &end);
    if (ret == 1) {
      hts_log_error("Error: parse_cpg_line()");
      return 1;
    }
    //load genome_summary_map
    chr_itor = ctx_sum.genome_wide_map.find(chr);
    if (chr_itor == ctx_sum.genome_wide_map.end()) {
      map<mhap_pos_t, summary_t> summary_t_map;
      summary_t sum_t  = {0,0,0,0,0,0,0,0,0,0};
      summary_t_map[beg] = sum_t;
      ctx_sum.genome_wide_map[chr] = summary_t_map;
    } else {
      cpg_itor = ctx_sum.genome_wide_map[chr].find(beg);
      if (cpg_itor == ctx_sum.genome_wide_map[chr].end()) {
        summary_t sum_t  = {0,0,0,0,0,0,0,0,0,0};
        ctx_sum.genome_wide_map[chr][beg] = sum_t;
      }
    }
  }
  return 0;
}

int saving_genome_wide(ContextSummary &ctx_sum) {
  string out_stream_name;
  if (ctx_sum.fn_out != NULL) {
    out_stream_name = ctx_sum.fn_out;
  } else {
    out_stream_name = "summary_genome_wide.txt";
  }
  ofstream out_stream(out_stream_name);

  map<string, map<mhap_pos_t, summary_t> >::iterator chr_itor;
  map<mhap_pos_t, summary_t>::iterator cpg_itor;
  for (chr_itor = ctx_sum.genome_wide_map.begin(); chr_itor != ctx_sum.genome_wide_map.end(); chr_itor++) {
    for (cpg_itor = chr_itor->second.begin(); cpg_itor != chr_itor->second.end(); cpg_itor++) {
      if (!cpg_itor->second.is_empty()) {
        if (ctx_sum.stranded) {
          if (!cpg_itor->second.is_direction_plus_empty()) {
            out_stream << chr_itor->first + '\t' + to_string(cpg_itor->first) +
                '\t' + to_string(cpg_itor->first + 1) + '\t' + '+' + '\t' + to_string(cpg_itor->second.n_reads)
                + '\t' + to_string(cpg_itor->second.m_base) + '\t' + to_string(cpg_itor->second.t_base) + '\t' +
                to_string(cpg_itor->second.n_reads_k4) + '\t' + to_string(cpg_itor->second.n_dr) << endl;
          }
          if (!cpg_itor->second.is_direction_minus_empty()) {
            out_stream << chr_itor->first + '\t' + to_string(cpg_itor->first) +
                '\t' + to_string(cpg_itor->first + 1) + '\t' + '-' + '\t' + to_string(cpg_itor->second.n_reads_r)
                + '\t' + to_string(cpg_itor->second.m_base_r) + '\t' + to_string(cpg_itor->second.t_base_r) + '\t' +
                to_string(cpg_itor->second.n_reads_k4_r) + '\t' + to_string(cpg_itor->second.n_dr_r) << endl;
          }
        } else {
          out_stream << chr_itor->first + '\t' + to_string(cpg_itor->first) +
              '\t' + to_string(cpg_itor->first + 1) + '\t' + '*' + '\t' + to_string(cpg_itor->second.n_reads)
              + '\t' + to_string(cpg_itor->second.m_base) + '\t' + to_string(cpg_itor->second.t_base) + '\t' +
              to_string(cpg_itor->second.n_reads_k4) + '\t' + to_string(cpg_itor->second.n_dr) << endl;
        }
      }
    }
  }
  out_stream.close();
  return 0;
}

int process_genome_wide(ContextSummary &ctx_sum) {
  ctx_sum.fp_hap = mhap_open(ctx_sum.fn_hap, "r");
  if (ctx_sum.fp_hap == NULL) {
    hts_log_error("Fail to open mhap file");
    return 1;
  }
  mhap_t hap_line_t = mhap_t {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};
  while(mhap_read(ctx_sum.fp_hap, &hap_line_t) == 0) {
    mhap_pos_t cur_m_base = 0;
    mhap_pos_t cur_t_base = 0;
    map<string, map<mhap_pos_t, summary_t> >::iterator chr_itor;
    map<mhap_pos_t, summary_t>::iterator cpg_itor;
    chr_itor = ctx_sum.genome_wide_map.find(hap_line_t.chr);
    if (chr_itor == ctx_sum.genome_wide_map.end()) {
      hts_log_error("Can not find chr: %s, CpG file and mhap file are not match.", hap_line_t.chr.c_str());
      return 1;
    }
    cpg_itor = ctx_sum.genome_wide_map[hap_line_t.chr].find(hap_line_t.chr_beg);
    if (cpg_itor == ctx_sum.genome_wide_map[hap_line_t.chr].end()) {
      hts_log_error("Can not find CpG begin point %lld, CpG file and mhap file are not match.", hap_line_t.chr_beg);
      return 1;
    }
    summary_t cur_sum_t = summary_t {0,0,0,0,0,0,0,0,0,0};
    if (ctx_sum.stranded) {
      if (hap_line_t.mhap_direction == '+') {
        ++cur_sum_t.n_reads;
        cur_sum_t.t_base += hap_line_t.mhap_str.size();
        cur_t_base = hap_line_t.mhap_str.size();
        for (auto b : hap_line_t.mhap_str) {
          if (b == '1') {
            ++cur_sum_t.m_base;
            ++cur_m_base;
          }
        }
        if (cur_m_base >= 4) {
          ++cur_sum_t.n_reads_k4;
        }
        if (cur_m_base > 0 && cur_t_base != cur_m_base) {
          ++cur_sum_t.n_dr;
        }
      } else if (hap_line_t.mhap_direction == '-') {
        ++cur_sum_t.n_reads_r;
        cur_sum_t.t_base_r += hap_line_t.mhap_str.size();
        cur_t_base = hap_line_t.mhap_str.size();
        for (auto b : hap_line_t.mhap_str) {
          if (b == '1') {
            ++cur_sum_t.m_base_r;
            ++cur_m_base;
          }
        }
        if (cur_m_base >= 4) {
          ++cur_sum_t.n_reads_k4_r;
        }
        if (cur_m_base > 0 && cur_t_base != cur_m_base) {
          ++cur_sum_t.n_dr_r;
        }
      } else {
        hts_log_error("Contain unknown direction info.");
        return 1;
      }
    } else {
      ++cur_sum_t.n_reads;
      cur_sum_t.t_base += hap_line_t.mhap_str.size();
      cur_t_base = hap_line_t.mhap_str.size();
      for (auto b : hap_line_t.mhap_str) {
        if (b == '1') {
          ++cur_sum_t.m_base;
          ++cur_m_base;
        }
      }
      if (cur_m_base >= 4) {
        ++cur_sum_t.n_reads_k4;
      }
      if (cur_m_base > 0 && cur_t_base != cur_m_base) {
        ++cur_sum_t.n_dr;
      }
    }

    for (int i = 0;
         cpg_itor != ctx_sum.genome_wide_map[hap_line_t.chr].end() &&
    i < hap_line_t.mhap_str.size(); cpg_itor++, i++) {

      cpg_itor->second.n_reads += cur_sum_t.n_reads;
      cpg_itor->second.m_base += cur_sum_t.m_base;
      cpg_itor->second.t_base += cur_sum_t.t_base;
      cpg_itor->second.n_reads_k4 += cur_sum_t.n_reads_k4;
      cpg_itor->second.n_dr += cur_sum_t.n_dr;

      if (ctx_sum.stranded) {
        cpg_itor->second.n_reads_r += cur_sum_t.n_reads_r;
        cpg_itor->second.m_base_r += cur_sum_t.m_base_r;
        cpg_itor->second.t_base_r += cur_sum_t.t_base_r;
        cpg_itor->second.n_reads_k4_r += cur_sum_t.n_reads_k4_r;
        cpg_itor->second.n_dr_r += cur_sum_t.n_dr_r;
      }
      //sanity check
      if (i == hap_line_t.mhap_str.size() - 1) {
        if (cpg_itor->first != hap_line_t.chr_end) {
          hts_log_error("CpG file and mhap file do not match.");
          return 1;
        }
      }
    }
  }
  mhap_close(ctx_sum.fp_hap);
  return 0;
}


int get_genome_wide(ContextSummary &ctx_sum) {
  int ret = 0;
  cout << "Loading CpG positions..." << endl;
  ret = load_cpg_init_map(ctx_sum);
  if (ret == 1) {
    return 1;
  }
  cout << "Processing..." << endl;
  ret = process_genome_wide(ctx_sum);
  if (ret == 1) {
    return 1;
  }

  return 0;
}

static void help() {
  cout << "Usage: mhaptools summary -i <in.mhap> -c <CpG.gz> [-r chr:beg-end | -b bed_file.bed ] | [-g] [-s] [-o name.mhap]" << endl;
  cout << "Options:" << endl;
  cout << "  -i  str  input file, mhap format" << endl;
  cout << "  -c  str  CpG file, gz format" << endl;
  cout << "  -r  str  region, e.g. chr1:2000-200000" << endl;
  cout << "  -b  str  bed file, one query region per line" << endl;
  cout << "  -g  flag get genome-wide result" << endl;
  cout << "  -s  flag group results by the direction of mhap reads" << endl;
  cout << "  -o  str  output filename [fn_hap_summary_genome_wide.txt | fn_hap_summary.txt]" << endl;
  cout << "Long options:" << endl;
  cout << "  -i  --input" << endl;
  cout << "  -c  --cpg" << endl;
  cout << "  -r  --region" << endl;
  cout << "  -b  --bed" << endl;
  cout << "  -g  --genome-wide" << endl;
  cout << "  -s  --stranded" << endl;
  cout << "  -o  --output" << endl;
  cout << "Examples:" << endl;
  cout << "- Get summary within a region:" << endl;
  cout << "  mhaptools summary -i in.mhap -c CpG.gz -r chr1:2000-200000" << endl << endl;
  cout << "- Get summary within several regions:" << endl;
  cout << "  mhaptools summary -i in.mhap -c CpG.gz -b bed_file.bed" << endl << endl;
  cout << "- Get genome-wide summary:" << endl;
  cout << "  mhaptools summary -i in.mhap -c CpG.gz -g" << endl << endl;
}

int sum_fn_suffix_check(ContextSummary &ctx_sum) {
  string gz_suffix = ".gz";
  string mhap_suffix = ".mhap";
  string mhap_gz_suffix = ".mhap.gz";
  string bed_suffix = ".bed";
  string txt_suffix = ".txt";
  if (ctx_sum.genome_wide) {
    if (ctx_sum.fn_hap) {
      if (!is_suffix(ctx_sum.fn_hap, mhap_suffix)) {
        hts_log_error("-i opt should be followed by a .mhap file.");
        return 1;
      }
    }
  } else {
    if (ctx_sum.fn_hap) {
      if (!is_suffix(ctx_sum.fn_hap, mhap_gz_suffix)) {
        hts_log_error("-i opt should be followed by a .mhap.gz file.");
        return 1;
      }
    }
  }
  if (ctx_sum.fn_bed) {
    if (!is_suffix(ctx_sum.fn_bed, bed_suffix)) {
      hts_log_error("-b opt should be followed by a .bed file.");
      return 1;
    }
  }
  if (ctx_sum.fn_cpg) {
    if (!is_suffix(ctx_sum.fn_cpg, gz_suffix)) {
      hts_log_error("-c opt should be followed by a .gz file.");
      return 1;
    }
  }
  if (ctx_sum.fn_out) {
    if (!is_suffix(ctx_sum.fn_out, txt_suffix)) {
      hts_log_error("-o opt should be followed by a .txt file.");
      return 1;
    }
  }
  return 0;
}

int main_summary(int argc, char *argv[]) {
  if (argc == optind) {
    help();
    return 0;
  }

  ContextSummary ctx_sum = ContextSummary();

  int long_index;

  static const char *opt_string = "i:o:b:sr:gc:h";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "output", optional_argument, NULL, 'o' },
      { "cpg", optional_argument, NULL, 'c' },
      { "bed", optional_argument, NULL, 'b' },
      { "stranded", optional_argument, NULL, 's' },
      { "region", optional_argument, NULL, 'r' },
      { "genome-wide", optional_argument, NULL, 'g' },
      { "help", optional_argument, NULL, 'h' },
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
      case 'c': {
        ctx_sum.fn_cpg = optarg;
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
      case 'g': {
        ctx_sum.genome_wide = true;
        break;
      }
      case 'h': {
        help();
        return 0;
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
  if (sum_fn_suffix_check(ctx_sum) == 1) {
    hts_log_error("filename suffix error.");
    return 1;
  }

  if (ctx_sum.genome_wide == false) {
    cout << "Processing..." << endl;
    ret = get_summary(ctx_sum);
    if (ret == 1) {
      hts_log_error("Get summary error.");
    }
    cout << "Saving..." << endl;
    saving_summary(ctx_sum);
  } else {
    ret = get_genome_wide(ctx_sum);
    if (ret == 1) {
      hts_log_error("Error: get_genome_wide()");
      return 1;
    }
    cout << "Saving..." << endl;
    ret = saving_genome_wide(ctx_sum);
    if (ret == 1) {
      hts_log_error("Error: saving_genome_wide()");
      return 1;
    }
  }
  return 0;
}
}//namespace std
