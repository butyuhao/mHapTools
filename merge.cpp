#include <iostream>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <fcntl.h>
#include <fstream>
#include <unistd.h>
#include <unordered_map>
#include "./htslib-1.10.2/htslib/kseq.h"
#include "./htslib-1.10.2/htslib/sam.h"
#include "./include/utils.h"
#include "./include/merge.h"


namespace std {

extern int _lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos);

vector<mhap_pos_t> get_cpg(ContextMerge &ctx_merge, mhap_t &hap_read) {
  /*
   * get cpg pos for a mhap read and store to mhap_t
   */
  vector<mhap_pos_t>_cpg_pos(hap_read.mhap_str.size(), 0);

  int pos = _lower_bound(ctx_merge.cpg_pos_map[hap_read.chr], hap_read.chr_beg);

  unordered_map<string, vector<hts_pos_t> >::iterator cpg_pos_map_itor;

  cpg_pos_map_itor = ctx_merge.cpg_pos_map.find(hap_read.chr);

  int i = 0;
  if (cpg_pos_map_itor != ctx_merge.cpg_pos_map.end()) {
    while (pos < ctx_merge.cpg_pos_map[hap_read.chr].size() && cpg_pos_map_itor != ctx_merge.cpg_pos_map.end()) {

      if (ctx_merge.cpg_pos_map[hap_read.chr][pos] >= hap_read.chr_beg &&
      ctx_merge.cpg_pos_map[hap_read.chr][pos] <= hap_read.chr_end) {
        if (i < hap_read.mhap_str.size()) {
          _cpg_pos[i] = ctx_merge.cpg_pos_map[hap_read.chr][pos];
        } else {
          _cpg_pos.push_back(ctx_merge.cpg_pos_map[hap_read.chr][pos]);
          hts_log_error("length of cpg pos vec and mhap str doesn't match");
        }
        i++;
        pos++;
      } else {
        break;
      }
    }
  }

  return _cpg_pos;
}

bool is_overlap(vector<mhap_pos_t> &cpg_pos_1, vector<mhap_pos_t> &cpg_pos_2,
                mhap_t &hap_1, mhap_t &hap_2, int *overlap_beg1, int *overlap_beg2,
                int *overlap_end1, int *overlap_end2) {
  if (hap_1.chr != hap_2.chr || hap_1.mhap_direction != hap_2.mhap_direction) {
    return false;
  }
  if (cpg_pos_1.size() != hap_1.mhap_str.size() || cpg_pos_2.size() != hap_2.mhap_str.size()) {
    hts_log_error("length of cpg pos vec and mhap str doesn't match");
  }
  bool break_flag = false;
  *overlap_beg1 = -1;
  *overlap_beg2 = -1;
  for(int i1 = 0; (i1 < cpg_pos_1.size()) && break_flag == false; i1++) {
    for(int i2 = 0; (i2 < cpg_pos_2.size()) && break_flag == false; i2++) {
      if(cpg_pos_1[i1] == cpg_pos_2[i2]) {
        *overlap_beg1 = i1;
        *overlap_beg2 = i2;
        break_flag = true;
      }
    }
  }
  if (*overlap_beg1 == cpg_pos_1.size() && *overlap_beg2 == cpg_pos_2.size()) {
    //can't find same cpg pos in both mhap strings
    return false;
  }
  if ((*overlap_beg1 == -1 || *overlap_beg2 == -1)) {
    //can't find same cpg pos in both mhap strings
    return false;
  }
  int i1 = *overlap_beg1;
  int i2 = *overlap_beg2;
  while(i1 < cpg_pos_1.size() && i2 < cpg_pos_2.size()) {
    if (cpg_pos_1[i1] != cpg_pos_2[i2]) {
      return false;
    }
    if (hap_1.mhap_str[i1] != hap_2.mhap_str[i2]) {
      return false;
    }
    ++i1;
    ++i2;
  }
  *overlap_end1 = i1 - 1;
  *overlap_end2 = i2 - 1;
  return true;
}


bool is_identity(mhap_t &hap_a, mhap_t &hap_b) {
  if (hap_a.chr == hap_b.chr && hap_a.chr_beg == hap_b.chr_beg &&
      hap_a.chr_end == hap_b.chr_end && hap_a.mhap_str == hap_b.mhap_str &&
      hap_a.mhap_direction == hap_b.mhap_direction) {
    return true;
  }
  return false;
}

mhap_t merge(mhap_t &hap_t_1, mhap_t &hap_t_2, int &overlap_beg_a,
             int &overlap_beg_b, int &overlap_end_a, int &overlap_end_b) {
  if (hap_t_1.chr_beg <= hap_t_2.chr_beg && hap_t_1.chr_end >= hap_t_2.chr_end) {
    return hap_t_1;
  }
  if (hap_t_1.chr_beg >= hap_t_2.chr_beg && hap_t_1.chr_end <= hap_t_2.chr_end) {
    return hap_t_2;
  }
  mhap_t hap_merge;
  string hap_str = "";
  if (hap_t_1.chr_beg <= hap_t_2.chr_beg) {
    hap_merge.chr_beg = hap_t_1.chr_beg;
    if (!(overlap_beg_a == 0 && overlap_beg_b == 0)) {
      hap_str += hap_t_1.mhap_str.substr(0, overlap_beg_a);
    }
  } else {
    hap_merge.chr_beg = hap_t_2.chr_beg;
    if (!(overlap_beg_a == 0 && overlap_beg_b == 0)) {
    hap_str += hap_t_2.mhap_str.substr(0, overlap_beg_b);
    }
  }

  hap_str += hap_t_1.mhap_str.substr(overlap_beg_a, overlap_end_a - overlap_beg_a + 1);

  if (hap_t_1.chr_end <= hap_t_2.chr_end) {
    hap_merge.chr_end = hap_t_2.chr_end;
    if (!(overlap_end_a == hap_t_1.mhap_str.size() - 1 &&
        overlap_end_b == hap_t_2.mhap_str.size() - 1)) {

      hap_str += hap_t_2.mhap_str.substr(overlap_end_b + 1);

    }
  } else {
    hap_merge.chr_end = hap_t_1.chr_end;
    if (!(overlap_end_a == hap_t_1.mhap_str.size() - 1 &&
          overlap_end_b == hap_t_2.mhap_str.size() - 1)) {

      hap_str += hap_t_1.mhap_str.substr(overlap_end_a + 1);

    }
  }
  hap_merge.chr = hap_t_1.chr;
  hap_merge.mhap_str = hap_str;
  hap_merge.mhap_count = 1;
  hap_merge.mhap_direction = hap_t_1.mhap_direction;
  return hap_merge;
}

bool load_chr_cpg(ContextMerge &ctx_merge) {
  /*
   * load all the cpg pos to ctx_merge
   */
  ctx_merge.fp_cpg = hts_open(ctx_merge.fn_cpg, "r");
  if (ctx_merge.fp_cpg == NULL) {
    return 1;
  }
  kstring_t cpg_line = {0,0,NULL};
  unordered_map<string, vector<hts_pos_t> >::iterator cpg_pos_map_itor;
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

  return 0;
}

bool merge_opt_check(ContextMerge &ctx_merge) {
  if (ctx_merge.fn_input_vec.size() < 2) {
    hts_log_error("Please specify at least two mhap.gz files to merge");
    return false;
  }
  return true;
}

bool comp_hap(const mhap_t &a, const mhap_t &b)
{
  if (a.chr == b.chr) {
    if (a.chr_beg != b.chr_beg) {
      return a.chr_beg < b.chr_beg;
    } else if (a.chr_end != b.chr_end) {
      return a.chr_end < b.chr_end;
    } else if (a.mhap_str != b.mhap_str) {
      return a.mhap_str < b.mhap_str;
    } else if (a.mhap_count!= b.mhap_count){
      return a.mhap_count < b.mhap_count;
    } else {
      return a.mhap_direction < b.mhap_direction;
    }
  } else {
    return strcmp(a.chr.c_str(), b.chr.c_str()) < 0;
  }
}

void saving_merged_hap(ContextMerge &ctx_merge, vector<mhap_t> &merge_result) {

  string out_stream_name;
  if (ctx_merge.fn_out) {
    out_stream_name = string(ctx_merge.fn_out).substr(0, strlen(ctx_merge.fn_out) - 3);
  } else {
    out_stream_name = "merged.mhap";
  }
  ofstream out_stream(out_stream_name);

  unordered_map<string, bool> is_overlap;
  unordered_map<string, bool>::iterator itor;
  vector<mhap_t>::iterator ht_itor;

  for (ht_itor = merge_result.begin(); ht_itor != merge_result.end(); ht_itor++) {
    string line = (*ht_itor).chr + '\t' + to_string((*ht_itor).chr_beg) + '\t' +
        to_string((*ht_itor).chr_end) + '\t' + (*ht_itor).mhap_str + '\t' +
        to_string((*ht_itor).mhap_count) + '\t' + (*ht_itor).mhap_direction;
    itor = is_overlap.find(line);
    if (itor == is_overlap.end()) {
      out_stream << line << '\n';
    }
    is_overlap[line] = true;
  }
  out_stream.close();

  // compress mhap to gz
  string path_src = out_stream_name;
  if (out_stream_name[0] != '.' && out_stream_name[0] != '/') {
    path_src = "./" + path_src;
  }
  string gz_path = path_src + ".gz";
  int f_src = open(path_src.c_str(), O_RDONLY);
  if (f_src < 0) {
    hts_log_error("bgzip: source file open error.");
  }
  int c = -1;
  static const int WINDOW_SIZE = 64 * 1024;
  void *buffer;
  buffer = malloc(WINDOW_SIZE);
  BGZF *fp;
  char out_mode[3] = "w\0";
  fp = bgzf_open(gz_path.c_str(), out_mode);
  while ((c = read(f_src, buffer, WINDOW_SIZE)) > 0) {
    if (bgzf_write(fp, buffer, c) < 0) hts_log_error("Could not write %d bytes: Error %d\n", c, fp->errcode);
  }
  bgzf_close(fp);
  free(buffer);
  close(f_src);
  // delete temp mhap file
  remove(path_src.c_str());
}

ContextMerge::~ContextMerge() {
  if (fp_hap1_gz) {
    bgzf_close(fp_hap1_gz);
  }
  if (fp_hap2_gz) {
    bgzf_close(fp_hap2_gz);
  }
  if (fp_cpg) {
    hts_close(fp_cpg);
  }
}

static void help() {
  cout << "Usage: mhaptools merge -i <in1.mhap.gz in2.mhap.gz inx.mhap.gz> -c <CpG.gz> [-o name.mhap.gz]" << endl;
  cout << "Options:" << endl;
  cout << "  -i  str  input file, multiple .mhap.gz files to merge" << endl;
  cout << "  -c  str  CpG file, gz format" << endl;
  cout << "  -o  str  output filename [out.mhap.gz]" << endl;
  cout << "Long options:" << endl;
  cout << "  -i  --input" << endl;
  cout << "  -c  --cpg" << endl;
  cout << "  -o  --output" << endl;
  cout << "Examples:" << endl;
  cout << "- Merge two .mhap.gz files:" << endl;
  cout << "  mhaptools merge -i in1.mhap.gz in2.mhap.gz -c CpG.gz" << endl << endl;

}

int convert_fn_suffix_check(ContextMerge &ctx_merge) {
  string gz_suffix = ".gz";
  string mhap_suffix = ".mhap.gz";
  string bed_suffix = ".bed";
  if (ctx_merge.fn_cpg) {
    if (!is_suffix(ctx_merge.fn_cpg, gz_suffix)) {
      hts_log_error("-c opt should be followed by a .gz file.");
      return 1;
    }
  }
  if (ctx_merge.fn_out) {
    if (!is_suffix(ctx_merge.fn_out, mhap_suffix)) {
      hts_log_error("-o opt should be followed by a .mhap.gz file.");
      return 1;
    }
  }
  for (int i = 0; i < ctx_merge.fn_input_vec.size(); i++) {
    if (ctx_merge.fn_input_vec[i].size() < mhap_suffix.size()) {
      return 1;
    }
    string suffix = ctx_merge.fn_input_vec[i].substr(ctx_merge.fn_input_vec[i].size() - mhap_suffix.size(), mhap_suffix.size());
    if (suffix != mhap_suffix) {
      string error_string = "File format error, replace this file with a " + mhap_suffix + " file";
      hts_log_error("%s", error_string.c_str());
      hts_log_error("%s", ctx_merge.fn_input_vec[i].c_str());
      return 1;
    }
  }
  return 0;
}

int main_merge(int argc, char *argv[]) {
  if (argc == optind) {
    help();
    return 0;
  }

  ContextMerge ctx_merge = ContextMerge();

  int long_index;

  static const char *opt_string = "i:o:c:h";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "output", optional_argument, NULL, 'o' },
      { "cpg", required_argument, NULL, 'c' },
      { "help", optional_argument, NULL, 'h' },
      { NULL, no_argument, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        if (argc < 4) {
          return 1;
        } else {
          for (int i = optind - 1; i < argc; i++) {
            if (argv[i][0] == '-') {
              break;
            }
            ctx_merge.fn_input_vec.push_back(argv[i]);
          }
        }
        break;
      }
      case 'c': {
        ctx_merge.fn_cpg = optarg;
        break;
      }
      case 'o': {
        ctx_merge.fn_out = optarg;
        break;
      }
      case 'h': {
        help();
        return 0;
      }
      default: {
        break;
      }
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }

  if (!merge_opt_check(ctx_merge)) {
    hts_log_error("opt error");
    return 1;
  }


  if (convert_fn_suffix_check(ctx_merge) == 1) {
    hts_log_error("filename suffix error.");
    return 1;
  }

  ctx_merge.fn_hap1 = ctx_merge.fn_input_vec[0];

  for (int i_merge = 1; i_merge < ctx_merge.fn_input_vec.size(); i_merge++) {

    ctx_merge.fn_hap2 = ctx_merge.fn_input_vec[i_merge];

    ctx_merge.fp_hap1_gz = bgzf_open(ctx_merge.fn_hap1.c_str(), "r");
    ctx_merge.fp_hap2_gz = bgzf_open(ctx_merge.fn_hap2.c_str(), "r");

    if (ctx_merge.fp_hap1_gz == NULL) {
      string error_string = "Fail to open " + string(ctx_merge.fn_hap1);
      hts_log_error("%s", error_string.c_str());
      return 0;
    }

    if (ctx_merge.fp_hap2_gz == NULL) {
      string error_string = "Fail to open " + string(ctx_merge.fn_hap2);
      hts_log_error("%s", error_string.c_str());
      return 0;
    }
    int ret = 0;

    if (ctx_merge.is_load_cpg) {
      cout << "Loading CpG positions..." << endl;
      ret = load_chr_cpg(ctx_merge);

      if (ret == 1) {
        return 1;
      }
      ctx_merge.is_load_cpg = false;
    }

    mhap_t hap_t_1 = {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};
    mhap_t hap_t_2 = {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};

    vector<mhap_t> merge_result_vec;
    vector<mhap_t> hap_to_merge;

    cout << "Loading mHap files..." << endl;

    kstring_t str = {0, 0, NULL};

    while (bgzf_getline(ctx_merge.fp_hap1_gz, '\n', &str) >= 0) {
      parse_mhap_line(str.s, str.l, &hap_t_1);
      hap_to_merge.push_back(hap_t_1);
    }
    while (bgzf_getline(ctx_merge.fp_hap2_gz, '\n', &str) >= 0) {
      parse_mhap_line(str.s, str.l, &hap_t_2);
      hap_to_merge.push_back(hap_t_2);
    }

    //sort
    cout << "Sorting..." << endl;
    sort(hap_to_merge.begin(), hap_to_merge.end(), comp_hap);

    int i = 1;
    vector<mhap_pos_t> cpg_pos_a;
    vector<mhap_pos_t> cpg_pos_b;
    mhap_t hap_a = {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};
    mhap_t hap_b = {MHAP_NULL_STRING, 0, 0, MHAP_NULL_STRING, 0, MHAP_DEFAULT_DIRECTION};
    mhap_t merge_result;
    cout << "Processing..." << endl;
    if (hap_to_merge.size() == 0) {
      hts_log_error("mhap files are empty");
    } else if (hap_to_merge.size() == 1) {
      merge_result_vec.push_back(hap_to_merge[0]);
    } else if (hap_to_merge.size() >= 2) {
      hap_b = hap_to_merge[0];
      int overlap_beg_a, overlap_beg_b, overlap_end_a, overlap_end_b;
      while (i < hap_to_merge.size()) {
        hap_a = hap_b;
        hap_b = hap_to_merge[i];

        cpg_pos_a.clear();
        cpg_pos_b.clear();
        cpg_pos_a = get_cpg(ctx_merge, hap_a);
        cpg_pos_b = get_cpg(ctx_merge, hap_b);

        if (is_identity(hap_a, hap_b)) {

          hap_b.mhap_count += hap_a.mhap_count;

        } else {
          merge_result_vec.push_back(hap_a);
        }
        i++;
      }
      merge_result_vec.push_back(hap_b);
    }
    cout << "Saving..." << endl;
    saving_merged_hap(ctx_merge, merge_result_vec);

    if (ctx_merge.fn_out) {
      ctx_merge.fn_hap1 = ctx_merge.fn_out;
    } else {
      ctx_merge.fn_hap1 = "merged.mhap.gz";
    }

  }

  return 0;
}
}//namespace std
