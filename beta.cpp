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


namespace std {

extern int _lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos);

ContextBeta::~ContextBeta() {
  if (fp_hap) {
    hap_close(fp_hap);
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
  map<string, map<hap_pos_t, beta_t> >::iterator chr_itor;
  map<hap_pos_t, beta_t>::iterator cpg_itor;
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
        beta_t bt_t  = {0, 0,0,0};
        map<hap_pos_t, beta_t> m;
        m[cpg_start] = bt_t;
        ctx_beta.beta_map[chr] = m;
    } else {
        cpg_itor = ctx_beta.beta_map[chr].find(cpg_start);
        if (cpg_itor == ctx_beta.beta_map[chr].end()) {
          beta_t bt_t  = {0, 0,0,0};
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

bool saving_beta(ContextBeta &ctx_beta) {
  string out_stream_name;
  if (ctx_beta.fn_out != NULL) {
    out_stream_name = ctx_beta.fn_out;
  } else {
    out_stream_name = "beta.txt";
  }
  ofstream out_stream(out_stream_name);

  map<string, map<hap_pos_t, beta_t> >::iterator chr_itor;
  map<hap_pos_t, beta_t>::iterator cpg_itor;
  if (ctx_beta.stranded) {
      for (chr_itor = ctx_beta.beta_map.begin(); chr_itor != ctx_beta.beta_map.end(); chr_itor++) {
        for (cpg_itor = chr_itor->second.begin(); cpg_itor != chr_itor->second.end(); cpg_itor++) {
          if (cpg_itor->second.total_reads != 0) {
            out_stream << chr_itor->first << '\t' << cpg_itor->first << '\t'
                       << cpg_itor->first + 1 << '\t' <<cpg_itor->second.methy_reads
                       << '\t' << cpg_itor->second.total_reads << '\t' << '+' << endl;
          }
          if (cpg_itor->second.total_reads_r != 0) {
            out_stream << chr_itor->first << '\t' << cpg_itor->first << '\t'
                       << cpg_itor->first + 1 << '\t' <<cpg_itor->second.methy_reads_r
                       << '\t' << cpg_itor->second.total_reads_r << '\t' << '-' << endl;
          }
        }
      }
    } else if (!ctx_beta.stranded){
      for (chr_itor = ctx_beta.beta_map.begin(); chr_itor != ctx_beta.beta_map.end(); chr_itor++) {
        for (cpg_itor = chr_itor->second.begin(); cpg_itor != chr_itor->second.end(); cpg_itor++) {
          if (cpg_itor->second.total_reads != 0) {
            out_stream << chr_itor->first << '\t' << cpg_itor->first << '\t'
            << cpg_itor->first + 1 << '\t' <<cpg_itor->second.methy_reads
            << '\t' << cpg_itor->second.total_reads << '\t' << '*' << endl;
          }
        }
      }
    }
  out_stream.close();
  return 0;
}

int process_beta(ContextBeta &ctx_beta, hap_t &h_t) {
  map<string, map<hap_pos_t, beta_t> >::iterator chr_itor;

  map<hap_pos_t, beta_t>::iterator cpg_itor;

  unordered_map<string, vector<hts_pos_t>>::iterator cpg_pos_map_itor;

  //check the existence of the chr in the map
  chr_itor = ctx_beta.beta_map.find(h_t.chr);
  if (chr_itor == ctx_beta.beta_map.end()) {
    hts_log_error("Can't find chr: %s in the input CpG file.", h_t.chr.c_str());
    return 1;
  }

  int pos = _lower_bound(ctx_beta.cpg_pos_map[h_t.chr], h_t.chr_beg);
  cpg_pos_map_itor = ctx_beta.cpg_pos_map.find(h_t.chr);

  int i = 0;
  if (cpg_pos_map_itor != ctx_beta.cpg_pos_map.end()) {
    while (pos < ctx_beta.cpg_pos_map[h_t.chr].size() && cpg_pos_map_itor != ctx_beta.cpg_pos_map.end()) {
      hap_pos_t cur_cpg_pos = ctx_beta.cpg_pos_map[h_t.chr][pos];
      if (cur_cpg_pos >= h_t.chr_beg && cur_cpg_pos <= h_t.chr_end) {
        if (i < h_t.hap_str.size()) {
          hap_pos_t cur_cpg_pos = ctx_beta.cpg_pos_map[h_t.chr][pos];
          char cur_hap = h_t.hap_str[i];
          //check current cpg pos in case it doesn't exists in the map
          cpg_itor = ctx_beta.beta_map[h_t.chr].find(cur_cpg_pos);
          if (cpg_itor == ctx_beta.beta_map[h_t.chr ].end()) {
            hts_log_error("Can't find the CpG position in the input CpG file. Position: %lld ", cur_cpg_pos);
            hts_log_error("hap read: %s %lld %lld %s %c", h_t.chr.c_str(), h_t.chr_beg, h_t.chr_end, h_t.hap_str.c_str(), h_t.hap_direction);
            return 1;
          }
          if (!ctx_beta.stranded) {
            if (cur_hap == '1') {
              ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads;
              ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].methy_reads;
            } else if (cur_hap == '0') {
              ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads;
            } else {
              hts_log_error("hap string value error.");
              return 1;
            }
          } else if (ctx_beta.stranded) {
            switch (h_t.hap_direction) {
              case '+':
                if (cur_hap == '1') {
                  ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads;
                  ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].methy_reads;
                } else if (cur_hap == '0') {
                  ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads;
                } else {
                  hts_log_error("hap string value error.");
                  return 1;
                }
                break;
              case '-':
                if (cur_hap == '1') {
                  ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads_r;
                  ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].methy_reads_r;
                } else if (cur_hap == '0') {
                  ++ctx_beta.beta_map[h_t.chr][cur_cpg_pos].total_reads_r;
                } else {
                  hts_log_error("hap string value error.");
                  return 1;
                }
                break;
              case '*':
                hts_log_error("Hap file contains directional reads, do not use --stranded or -s opt");
                return 1;
                break;
              default:
                hts_log_error("hap string value error.");
                return 1;
                break;
            }
          }
        } else {
          hts_log_error("length of cpg pos and hap str doesn't match in hap read.");
          hts_log_error("hap read: %s %lld %lld %s %c", h_t.chr.c_str(), h_t.chr_beg, h_t.chr_end, h_t.hap_str.c_str(), h_t.hap_direction);
          return 1;
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
  hap_t h_t = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};
  int ret = 0;
  if (ctx_beta.fn_bed == NULL) {
    while(hap_read(ctx_beta.fp_hap, &h_t) == 0) {
      ret = process_beta(ctx_beta, h_t);
      if (ret == 1) {
        return 1;
      }
    }
  } else {
    while(hap_read(ctx_beta.fp_hap, &h_t) == 0) {
      regidx_t *idx = regidx_init(ctx_beta.fn_bed,NULL,NULL,0,NULL);
      regitr_t *itr = regitr_init(idx);
      while (regitr_loop(itr)) {
        region_t reg_t = region_t {"", 0,0};
        reg_t.chr = itr->seq;
        reg_t.beg =  itr->beg;
        reg_t.end =  itr->end + 1;
        if (h_t.chr == reg_t.chr &&
            h_t.chr_beg >= reg_t.beg &&
            h_t.chr_end <= reg_t.end) {
          ret = process_beta(ctx_beta, h_t);
          if (ret == 1) {
            return 1;
          }
        } else {
          continue;
        }
      }
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

static void help() {
  cout << "Usage: haptools beta -i <in.hap> -c <CpG.gz> [-b bed_file.bed] [-s] [-o name.txt]" << endl;
  cout << "Options:" << endl;
  cout << "  -i  str  input file, hap format" << endl;
  cout << "  -c  str  CpG file, gz format" << endl;
  cout << "  -b  str  bed file, contains query regions" << endl;
  cout << "  -s  flag group results by the direction of hap reads" << endl;
  cout << "  -o  str  output file name [beta.txt]" << endl;
  cout << "Long options:" << endl;
  cout << "  -i  --input" << endl;
  cout << "  -c  --cpg" << endl;
  cout << "  -b  --bed" << endl;
  cout << "  -s  --stranded" << endl;
  cout << "  -o  --output" << endl;
  cout << "Examples:" << endl;
  cout << "- Get beta results:" << endl;
  cout << "  samtools beta -i in.hap -c CpG.gz" << endl << endl;
  cout << "- Get beta results, group results by the direction of hap reads:" << endl;
  cout << "  samtools beta -i in.hap -c CpG.gz -s" << endl << endl;
  cout << "- Get beta results within several regions" << endl;
  cout << "  samtools beta -i in.hap -c CpG.gz -b bed_file.bed" << endl;

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

  int ret = 0;

  cout << "Loding CpG positions..." << endl;
  ret = load_cpg_init_beta_map(ctx_beta);

  if (ret == 1) {
    hts_log_error("Could not load CpG file.");
    return 1;
  }

  ctx_beta.fp_hap = hap_open(ctx_beta.fn_hap, "rb");

  if (ctx_beta.fp_hap == NULL) {
    hts_log_error("Fail to open hap file.");
    return 0;
  }

  cout << "Processing..." << endl;
  ret = get_beta(ctx_beta);
  if (ret == 1) {
    return 1;
  }

  cout << "Saving..." << endl;
  ret = saving_beta(ctx_beta);
  if (ret == 1) {
    return 1;
  }

  return 0;
}

} //namespace std