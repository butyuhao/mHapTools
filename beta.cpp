//
// Created by Yuhao Dan on 2020/7/26.
//
#include <getopt.h>
#include <map>
#include <fstream>
#include <regidx.h>
#include "./include/beta.h"
#include "./htslib-1.10.2/htslib/kseq.h"
#include "./htslib-1.10.2/htslib/hts.h"
#include "./htslib-1.10.2/htslib/regidx.h"
#include "./include/summary.h"
#include "./include/utils.h"

namespace std {

extern int _lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos);

ContextBeta::~ContextBeta() {
  if (fp_hap_gz) {
    bgzf_close(fp_hap_gz);
  }
  if (fp_cpg) {
    hts_close(fp_cpg);
  }
  if (fp_bed) {
    fclose(fp_bed);
  }
}

bool load_cpg_init_beta_map(ContextBeta &ctx_beta) {
  /*
   * load all the cpg pos to ctx_beta and init beta map
   */
  ctx_beta.fp_cpg = hts_open(ctx_beta.fn_cpg, "r");
  if (ctx_beta.fp_cpg == NULL) {
    hts_log_error("CpG file pointer error.");
    return 1;
  }
  kstring_t cpg_line = {0,0,NULL};
  map<string, map<mhap_pos_t, beta_t> >::iterator chr_itor;
  map<mhap_pos_t, beta_t>::iterator cpg_itor;
  unordered_map<string, vector<hts_pos_t> >::iterator cpg_pos_map_itor;
  while (hts_getline(ctx_beta.fp_cpg, KS_SEP_LINE, &cpg_line) > 0) {

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

    //load beta_map
    chr_itor = ctx_beta.beta_map.find(chr);
    if (chr_itor == ctx_beta.beta_map.end()) {
        beta_t bt_t  = {0, 0,0,0,false};
        map<mhap_pos_t, beta_t, less<mhap_pos_t> > m;
        m[cpg_start] = bt_t;
        ctx_beta.beta_map[chr] = m;
    } else {
        cpg_itor = ctx_beta.beta_map[chr].find(cpg_start);
        if (cpg_itor == ctx_beta.beta_map[chr].end()) {
          beta_t bt_t  = {0, 0,0,0, false};
          ctx_beta.beta_map[chr][cpg_start] = bt_t;
        }
    }

    //load cpg_pos_map
    cpg_pos_map_itor = ctx_beta.cpg_pos_map.find(chr);
    if (cpg_pos_map_itor == ctx_beta.cpg_pos_map.end()) {
      vector<hts_pos_t> v;
      v.push_back(cpg_start);
      ctx_beta.cpg_pos_map[chr] = v;
    } else {
      ctx_beta.cpg_pos_map[chr].push_back(cpg_start);
    }
  }

  return 0;
}

bool saving_beta(ContextBeta &ctx_beta, int mode) {
  //mode==1 whole 输出total reads !=0 的结果
  //mode==2 bed 输出is_in_bed==true的结果
  string out_stream_name;
  if (ctx_beta.fn_out != NULL) {
    out_stream_name = ctx_beta.fn_out;
  } else {
    out_stream_name = "beta.txt";
  }
  ofstream out_stream(out_stream_name);

  map<string, map<mhap_pos_t, beta_t> , less<string> >::iterator chr_itor;
  map<mhap_pos_t, beta_t, less<mhap_pos_t>>::iterator cpg_itor;
  map<string, beta_t, less<string> >::iterator beta_with_bed_result_itor;
  if (mode==1) {
    if (ctx_beta.stranded) {
      for (chr_itor = ctx_beta.beta_map.begin(); chr_itor != ctx_beta.beta_map.end(); chr_itor++) {
        for (cpg_itor = chr_itor->second.begin(); cpg_itor != chr_itor->second.end(); cpg_itor++) {
          if (cpg_itor->second.total_reads != 0) {
            out_stream << chr_itor->first << '\t' << cpg_itor->first << '\t'
                       << cpg_itor->first + 1 << '\t' <<
                       cpg_itor->second.methy_reads * 100 / cpg_itor->second.total_reads <<
                       '\t' << cpg_itor->second.methy_reads << '\t' <<
                       cpg_itor->second.total_reads - cpg_itor->second.methy_reads <<
                       '\t' << '+' << endl;
          }
          if (cpg_itor->second.total_reads_r != 0) {
            out_stream << chr_itor->first << '\t' << cpg_itor->first << '\t'
                       << cpg_itor->first + 1 << '\t' <<
                       cpg_itor->second.methy_reads_r * 100 / cpg_itor->second.total_reads_r <<
                       '\t' << cpg_itor->second.methy_reads_r <<'\t' <<
                       cpg_itor->second.total_reads_r - cpg_itor->second.methy_reads_r<<
                       '\t' << '-' << endl;
          }
        }
      }
    } else if (!ctx_beta.stranded){
      for (chr_itor = ctx_beta.beta_map.begin(); chr_itor != ctx_beta.beta_map.end(); chr_itor++) {
        for (cpg_itor = chr_itor->second.begin(); cpg_itor != chr_itor->second.end(); cpg_itor++) {
          if (cpg_itor->second.total_reads != 0) {
            out_stream << chr_itor->first << '\t' << cpg_itor->first << '\t'
                       << cpg_itor->first + 1 << '\t' <<
                       cpg_itor->second.methy_reads * 100 / cpg_itor->second.total_reads <<
                       '\t' << cpg_itor->second.methy_reads << '\t' <<
                       cpg_itor->second.total_reads - cpg_itor->second.methy_reads <<
                       '\t' << '*' << endl;
          }
        }
      }
    }
  } else if (mode == 2) {
    beta_with_bed_result_itor = ctx_beta.beta_with_bed_results.begin();
    while (beta_with_bed_result_itor != ctx_beta.beta_with_bed_results.end()) {
      if (ctx_beta.stranded) {
        if (beta_with_bed_result_itor->second.total_reads != 0) {
          out_stream << beta_with_bed_result_itor->first << '\t' <<
                     beta_with_bed_result_itor->second.methy_reads * 100 / beta_with_bed_result_itor->second.total_reads
                     <<
                     '\t' << beta_with_bed_result_itor->second.methy_reads << '\t' <<
                     beta_with_bed_result_itor->second.total_reads - beta_with_bed_result_itor->second.methy_reads <<
                     '\t' << '+' << endl;
        } else {
          out_stream << beta_with_bed_result_itor->first << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << '+' << endl;
        }
        if (beta_with_bed_result_itor->second.total_reads_r != 0) {
          out_stream << beta_with_bed_result_itor->first << '\t' <<
                     beta_with_bed_result_itor->second.methy_reads_r * 100
                         / beta_with_bed_result_itor->second.total_reads_r <<
                     '\t' << beta_with_bed_result_itor->second.methy_reads_r << '\t' <<
                     beta_with_bed_result_itor->second.total_reads_r - beta_with_bed_result_itor->second.methy_reads_r
                     <<
                     '\t' << '-' << endl;
        } else {
          out_stream << beta_with_bed_result_itor->first << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << '-' << endl;
        }

      } else if (!ctx_beta.stranded) {
        if (beta_with_bed_result_itor->second.total_reads != 0) {
          out_stream << beta_with_bed_result_itor->first << '\t' <<
                     beta_with_bed_result_itor->second.methy_reads * 100 / beta_with_bed_result_itor->second.total_reads
                     <<
                     '\t' << beta_with_bed_result_itor->second.methy_reads << '\t' <<
                     beta_with_bed_result_itor->second.total_reads - beta_with_bed_result_itor->second.methy_reads <<
                     '\t' << '*' << endl;
        } else {
          out_stream << beta_with_bed_result_itor->first << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << '*' << endl;
        }
      }
      beta_with_bed_result_itor++;
    }

  }

  out_stream.close();
  return 0;
}

int process_beta(ContextBeta &ctx_beta, mhap_t &h_t) {
  map<string, map<mhap_pos_t, beta_t>, less<string> >::iterator chr_itor;

  map<mhap_pos_t, beta_t, less<mhap_pos_t> >::iterator cpg_itor;

  unordered_map<string, vector<hts_pos_t>>::iterator cpg_pos_map_itor;

  if (ctx_beta.is_process_beta_checked == false) {
    //check the existence of the first appear chr in the map
    chr_itor = ctx_beta.beta_map.find(h_t.chr);
    if (chr_itor == ctx_beta.beta_map.end()) {
      hts_log_error("Can't find chr: %s in the input CpG file.", h_t.chr.c_str());
      return 1;
    }
    ctx_beta.is_process_beta_checked = true;
  }

  int pos = _lower_bound(ctx_beta.cpg_pos_map[h_t.chr], h_t.chr_beg);
  cpg_pos_map_itor = ctx_beta.cpg_pos_map.find(h_t.chr);

  int i = 0;
  if (cpg_pos_map_itor != ctx_beta.cpg_pos_map.end()) {
    while (pos < ctx_beta.cpg_pos_map[h_t.chr].size() && cpg_pos_map_itor != ctx_beta.cpg_pos_map.end()) {
      mhap_pos_t cur_cpg_pos = ctx_beta.cpg_pos_map[h_t.chr][pos];
      if (cur_cpg_pos >= h_t.chr_beg && cur_cpg_pos <= h_t.chr_end) {
        if (i < h_t.mhap_str.size()) {
          char cur_hap = h_t.mhap_str[i];
          //check current cpg pos in case it doesn't exists in the map
          cpg_itor = ctx_beta.beta_map[h_t.chr].find(cur_cpg_pos);
          if (cpg_itor == ctx_beta.beta_map[h_t.chr ].end()) {
            hts_log_error("Can't find the CpG position in the input CpG file. Position: %lld ", cur_cpg_pos);
            hts_log_error("mhap read: %s %lld %lld %s %c", h_t.chr.c_str(), h_t.chr_beg, h_t.chr_end, h_t.mhap_str.c_str(), h_t.mhap_direction);
            return 1;
          }
          if (!ctx_beta.stranded) {
            if (cur_hap == '1') {
              ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads += h_t.mhap_count;
              ctx_beta.beta_map[h_t.chr][cur_cpg_pos].methy_reads += h_t.mhap_count;
            } else if (cur_hap == '0') {
              ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads += h_t.mhap_count;
            } else {
              hts_log_error("mhap string value error.");
              return 1;
            }
          } else if (ctx_beta.stranded) {
            switch (h_t.mhap_direction) {
              case '+':
                if (cur_hap == '1') {
                  ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads += h_t.mhap_count;
                  ctx_beta.beta_map[h_t.chr][cur_cpg_pos].methy_reads += h_t.mhap_count;
                } else if (cur_hap == '0') {
                  ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads += h_t.mhap_count;
                } else {
                  hts_log_error("mhap string value error.");
                  return 1;
                }
                break;
              case '-':
                if (cur_hap == '1') {
                  ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads_r += h_t.mhap_count;
                  ctx_beta.beta_map[h_t.chr][cur_cpg_pos].methy_reads_r += h_t.mhap_count;
                } else if (cur_hap == '0') {
                  ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads_r += h_t.mhap_count;
                } else {
                  hts_log_error("mhap string value error.");
                  return 1;
                }
                break;
              case '*':
                hts_log_error("mhap file contains directional reads, do not use --stranded or -s opt");
                return 1;
                break;
              default:
                hts_log_error("mhap string value error.");
                return 1;
                break;
            }
          }
        } else {
          hts_log_warning("length of cpg pos and mhap str doesn't match in mhap read.");
          hts_log_warning("mhap read: %s %lld %lld %s %c", h_t.chr.c_str(), h_t.chr_beg, h_t.chr_end, h_t.mhap_str.c_str(), h_t.mhap_direction);
          return 0;
        }
        i++;
        pos++;
      } else {
        break;
      }
    }
  }
  return 0;
}

int get_beta(ContextBeta &ctx_beta) {
  mhap_t h_t = {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};
  int ret = 0;
  kstring_t str = {0, 0, NULL};
  if (ctx_beta.fn_bed == NULL) {
    while(bgzf_getline(ctx_beta.fp_hap_gz, '\n', &str) >= 0) {
      parse_mhap_line(str.s, str.l, &h_t);
      ret = process_beta(ctx_beta, h_t);
      if (ret == 1) {
        return 1;
      }
    }
  } else {
    while(bgzf_getline(ctx_beta.fp_hap_gz, '\n', &str) >= 0) {
      parse_mhap_line(str.s, str.l, &h_t);
      ret = process_beta(ctx_beta, h_t);
      if (ret == 1) {
        return 1;
      }
    }
    map<string, map<mhap_pos_t, beta_t> , less<string> >::iterator chr_itor;
    map<mhap_pos_t, beta_t, less<mhap_pos_t>>::iterator cpg_itor;

    regidx_t *idx = regidx_init(ctx_beta.fn_bed,NULL,NULL,0,NULL);
    regitr_t *itr = regitr_init(idx);
    while (regitr_loop(itr)) {
      region_t reg_t = region_t {"", 0,0};
      string chr = itr->seq;
      mhap_pos_t beg =  itr->beg;
      mhap_pos_t end =  itr->end + 1;
      int pos = _lower_bound(ctx_beta.cpg_pos_map[chr], beg);

      string key = chr + '\t' + to_string(beg) + '\t' + to_string(end);

      beta_t bt_t  = {0, 0,0,0, false};

      mhap_pos_t  beg_cpg_pos = ctx_beta.cpg_pos_map[chr][pos]; // 找出比指定beg大的第一个位置

      chr_itor = ctx_beta.beta_map.find(chr);
      cpg_itor = chr_itor->second.find(beg_cpg_pos);
      while(cpg_itor->first <= end && cpg_itor!=chr_itor->second.end()) {
        bt_t.methy_reads += cpg_itor->second.methy_reads;
        bt_t.total_reads += cpg_itor->second.total_reads;
        bt_t.methy_reads_r += cpg_itor->second.methy_reads_r;
        bt_t.total_reads_r += cpg_itor->second.total_reads_r;
        cpg_itor++;
      }
      ctx_beta.beta_with_bed_results[key] = bt_t;
    }
  }

  return 0;
}

bool beta_opt_check(ContextBeta &ctx_beta) {
  if (ctx_beta.fn_cpg == NULL || ctx_beta.fn_hap == NULL) {
    hts_log_error("Please specify -i and -c");
    return 1;
  }
  return 0;
}

int beta_fn_suffix_check(ContextBeta &ctx_beta) {
  string mhap_suffix = ".mhap.gz";
  string gz_suffix = ".gz";
  string output_suffix = ".txt";
  string bed_suffix = ".bed";
  if (ctx_beta.fn_hap) {
    if (!is_suffix(ctx_beta.fn_hap, mhap_suffix)) {
      hts_log_error("-i opt should be followed by a .mhap.gz file.");
      return 1;
    }
  }
  if (ctx_beta.fn_cpg) {
    if (!is_suffix(ctx_beta.fn_cpg, gz_suffix)) {
      hts_log_error("-c opt should be followed by a .gz file.");
      return 1;
    }
  }
  if (ctx_beta.fn_out) {
    if (!is_suffix(ctx_beta.fn_out, output_suffix)) {
      hts_log_error("-o opt should be followed by a .txt file.");
      return 1;
    }
  }
  if (ctx_beta.fn_bed) {
    if (!is_suffix(ctx_beta.fn_bed, bed_suffix)) {
      hts_log_error("-b opt should be followed by a .bed file.");
      return 1;
    }
  }
  return 0;
}

static void help() {
  cout << "Usage: mhaptools beta -i <in.mhap.gz> -c <CpG.gz> [-b bed_file.bed] [-s] [-o name.txt]" << endl;
  cout << "Options:" << endl;
  cout << "  -i  str  input file, .mhap.gz format" << endl;
  cout << "  -c  str  CpG file, gz format" << endl;
  cout << "  -b  str  bed file, one query region per line" << endl;
  cout << "  -s  flag group results by the direction of mHap reads" << endl;
  cout << "  -o  str  output filename [beta.txt]" << endl;
  cout << "Long options:" << endl;
  cout << "  -i  --input" << endl;
  cout << "  -c  --cpg" << endl;
  cout << "  -b  --bed" << endl;
  cout << "  -s  --stranded" << endl;
  cout << "  -o  --output" << endl;
  cout << "Examples:" << endl;
  cout << "- Get beta results:" << endl;
  cout << "  mhaptools beta -i in.mhap.gz -c CpG.gz" << endl << endl;
  cout << "- Get beta results, group results by the direction of mHap reads:" << endl;
  cout << "  mhaptools beta -i in.mhap.gz -c CpG.gz -s" << endl << endl;
  cout << "- Get beta results within several regions" << endl;
  cout << "  mhaptools beta -i in.mhap.gz -c CpG.gz -b bed_file.bed" << endl << endl;

}

int main_beta(int argc, char *argv[]) {
  if (argc == optind) {
    help();
    return 0;
  }
  ContextBeta ctx_beta = ContextBeta();

  int long_index;

  static const char *opt_string = "i:o:c:b:sh";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "output", optional_argument, NULL, 'o' },
      { "cpg", required_argument, NULL, 'c' },
      { "bed", optional_argument, NULL, 'b' },
      { "stranded", optional_argument, NULL, 's' },
      { "help", optional_argument, NULL, 'h' },
      { NULL, no_argument, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        ctx_beta.fn_hap = optarg;
        break;
      }
      case 'c': {
        ctx_beta.fn_cpg = optarg;
        break;
      }
      case 'o': {
        ctx_beta.fn_out = optarg;
        break;
      }
      case 'b': {
        ctx_beta.fn_bed = optarg;
        break;
      }
      case 's': {
        ctx_beta.stranded = true;
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

  if (beta_opt_check(ctx_beta) == 1) {
    hts_log_error("opt error");
    return 1;
  }
  if (beta_fn_suffix_check(ctx_beta) == 1) {
    hts_log_error("filename suffix error.");
    return 1;
  }

  int ret = 0;

  ctx_beta.fp_hap_gz = bgzf_open(ctx_beta.fn_hap, "r");

  if (ctx_beta.fp_hap_gz == NULL) {
    hts_log_error("Fail to open input .mhap.gz file.");
    return 0;
  }

  cout << "Loding CpG positions..." << endl;
  ret = load_cpg_init_beta_map(ctx_beta);

  if (ret == 1) {
    hts_log_error("Could not load CpG file.");
    return 1;
  }

  cout << "Processing..." << endl;
  ret = get_beta(ctx_beta);
  if (ret == 1) {
    return 1;
  }

  cout << "Saving..." << endl;
  if (ctx_beta.fn_bed == NULL) {
    ret = saving_beta(ctx_beta, 1); // entire file
  } else {
    ret = saving_beta(ctx_beta, 2); // bed file
  }

  if (ret == 1) {
    return 1;
  }

  return 0;
}

} //namespace std