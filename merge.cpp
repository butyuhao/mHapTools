#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include "hap.cpp"
#include "merge.h"

namespace std {

extern int _lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos);

vector<hap_pos_t> get_cpg(ContextMerge &ctx_merge, hap_t &hap_read) {
  /*
   * get cpg pos for a hap read and store to hap_t
   */
  vector<hap_pos_t>_cpg_pos;

  int pos = _lower_bound(ctx_merge.cpg_pos_map[hap_read.chr], hap_read.chr_beg);

  unordered_map<string, vector<hts_pos_t>>::iterator cpg_pos_map_itor;

  cpg_pos_map_itor = ctx_merge.cpg_pos_map.find(hap_read.chr);
  cout << "get cpg while" << endl;

  if (cpg_pos_map_itor != ctx_merge.cpg_pos_map.end()) {
    while (pos < ctx_merge.cpg_pos_map[hap_read.chr].size() && cpg_pos_map_itor != ctx_merge.cpg_pos_map.end()) {

      if (ctx_merge.cpg_pos_map[hap_read.chr][pos] >= hap_read.chr_beg &&
      ctx_merge.cpg_pos_map[hap_read.chr][pos] <= hap_read.chr_end) {

        _cpg_pos.push_back(ctx_merge.cpg_pos_map[hap_read.chr][pos]);

        pos++;
      } else {
        break;
      }
    }
  }
  cout << "get cpg end" << endl;

  return _cpg_pos;
}

bool is_overlap(vector<hap_pos_t> &cpg_pos_1, vector<hap_pos_t> &cpg_pos_2,
    hap_t &hap_1, hap_t &hap_2, int *overlap_beg1, int *overlap_beg2,
    int *overlap_end1, int *overlap_end2) {
  if (hap_1.chr != hap_2.chr || hap_1.hap_direction != hap_2.hap_direction) {
    return false;
  }
  if (cpg_pos_1.size() != hap_1.hap_str.size() || cpg_pos_2.size() != hap_2.hap_str.size()) {
    hts_log_error("length of cpg pos vec and hap str doesn't match");
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
  if (*overlap_beg1 == cpg_pos_1.size() && *overlap_beg2 == cpg_pos_2.size() ||
      (*overlap_beg1 == -1 || *overlap_beg2 == -1)) {
    //can't find same cpg pos in both hap strings
    return false;
  }
  int i1 = *overlap_beg1;
  int i2 = *overlap_beg2;
  while(i1 < cpg_pos_1.size() && i2 < cpg_pos_2.size()) {
    if (cpg_pos_1[i1] != cpg_pos_2[i2]) {
      return false;
    }
    if (hap_1.hap_str[i1] != hap_2.hap_str[i2]) {
      return false;
    }
    ++i1;
    ++i2;
  }
  *overlap_end1 = i1 - 1;
  *overlap_end2 = i2 - 1;
  return true;
}


bool is_identity(hap_t &hap_a, hap_t &hap_b) {
  if (hap_a.chr == hap_b.chr && hap_a.chr_beg == hap_b.chr_beg &&
      hap_a.chr_end == hap_b.chr_end && hap_a.hap_str == hap_b.hap_str &&
      hap_a.hap_direction == hap_b.hap_direction) {
    return true;
  }
  return false;
}

hap_t merge(hap_t &hap_t_1, hap_t &hap_t_2, int &overlap_beg_a,
    int &overlap_beg_b, int &overlap_end_a, int &overlap_end_b) {
  cout << "merge" << endl;
  if (hap_t_1.chr_beg <= hap_t_2.chr_beg && hap_t_1.chr_end >= hap_t_2.chr_end) {
    return hap_t_1;
  }
  if (hap_t_1.chr_beg >= hap_t_2.chr_beg && hap_t_1.chr_end <= hap_t_2.chr_end) {
    return hap_t_2;
  }
  hap_t hap_merge;
  string hap_str = "";
  if (hap_t_1.chr_beg <= hap_t_2.chr_beg) {
    hap_merge.chr_beg = hap_t_1.chr_beg;
    if (!(overlap_beg_a == overlap_beg_b == 0)) {
      hap_str += hap_t_1.hap_str.substr(0, overlap_beg_a);
    }
  } else {
    hap_merge.chr_beg = hap_t_2.chr_beg;
    if (!(overlap_beg_a == overlap_beg_b == 0)) {
    hap_str += hap_t_2.hap_str.substr(0, overlap_beg_b);
    }
  }

  hap_str += hap_t_1.hap_str.substr(overlap_beg_a, overlap_end_a - overlap_beg_a + 1);

  if (hap_t_1.chr_end <= hap_t_2.chr_end) {
    hap_merge.chr_end = hap_t_2.chr_end;
    if (!(overlap_end_a == overlap_end_b
        == hap_t_1.hap_str.size() - 1 == hap_t_2.hap_str.size() - 1)) {

      hap_str += hap_t_2.hap_str.substr(overlap_end_b + 1);

    }
  } else {
    hap_merge.chr_end = hap_t_1.chr_end;
    if (!(overlap_end_a == overlap_end_b
        == hap_t_1.hap_str.size() - 1 == hap_t_2.hap_str.size() - 1)) {

      hap_str += hap_t_1.hap_str.substr(overlap_end_a + 1);

    }
  }
  hap_merge.chr = hap_t_1.chr;
  hap_merge.hap_str = hap_str;
  hap_merge.hap_direction = hap_t_1.hap_direction;
  cout << "end merge" << endl;
  return hap_merge;
}

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

bool comp_hap(const hap_t &a, const hap_t &b)
{
  if (a.chr == b.chr) {
    if (a.chr_beg != b.chr_beg) {
      return a.chr_beg < b.chr_beg;
    } else if (a.chr_end != b.chr_end) {
      return a.chr_end < b.chr_end;
    } else if (a.hap_str != b.hap_str) {
      return a.hap_str < b.hap_str;
    } else if (a.hap_count!= b.hap_count){
      return a.hap_count < b.hap_count;
    } else {
      return a.hap_direction < b.hap_direction;
    }
  } else {
    return strcmp(a.chr.c_str(), b.chr.c_str()) < 0;
  }
}

map<hap_pos_t, char> get_cpg_hap_str_map(vector<hap_pos_t> &cpg_pos, string &hap_str) {
  if (cpg_pos.size() != hap_str.size()) {
    cout << cpg_pos.size() << endl;
    cout << hap_str.size() << endl;
    hts_log_error("cpg pos length and hap str length don't match");
  }
  map<hap_pos_t, char> cpg_hap_str_map;
  for(int i = 0; i < cpg_pos.size(); i++) {
    cpg_hap_str_map[cpg_pos[i]] = hap_str[i];
  }
  return cpg_hap_str_map;
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
    return 1;
  }

  load_chr_cpg(ctx_merge);

  hapFile *fp_hap1 = hap_open(ctx_merge.fn_hap1, "rb");
  hapFile *fp_hap2 = hap_open(ctx_merge.fn_hap2, "rb");

  hap_t hap_t_1 = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};
  hap_t hap_t_2 = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};

//    //test
//
//    //get_cpg_pos
//
//  vector<hap_pos_t> mock_cpg_pos_1;
//  mock_cpg_pos_1.push_back(15513);
//  mock_cpg_pos_1.push_back(15526);
//
//  vector<hap_pos_t> mock_cpg_pos_2;
//  mock_cpg_pos_2.push_back(15643);
//
//  hap_t mock_hap_1;
//  mock_hap_1.chr = string("1");
//  mock_hap_1.chr_beg = 1851184;
//  mock_hap_1.chr_end = 1851264;
//  mock_hap_1.hap_str = string("11");
//  mock_hap_1.hap_direction = '+';
//
//  hap_t mock_hap_2;
//  mock_hap_2.chr = string("1");
//  mock_hap_2.chr_beg = 1851201;
//  mock_hap_2.chr_end = 1851283;
//  mock_hap_2.hap_str = string("0");
//  mock_hap_2.hap_direction = '+';
//
//  hap_t m_result;
//
//
//  int overlap_beg_a, overlap_beg_b, overlap_end_a, overlap_end_b;
//
//  if (is_overlap(mock_cpg_pos_1, mock_cpg_pos_2, mock_hap_1,
//      mock_hap_2, &overlap_beg_a, &overlap_beg_b, &overlap_end_a, &overlap_end_b)) {
//    cout << overlap_beg_a << endl;
//    cout << overlap_beg_b << endl;
//    m_result = merge(mock_hap_1, mock_hap_2, overlap_beg_a, overlap_beg_b, overlap_end_a, overlap_end_b);
//  }
//  m_result.print();
//while(1);
//  //test

  vector<hap_t> merge_result_vec;
  vector<hap_t> hap_to_merge;

  while(hap_read(fp_hap1, &hap_t_1) == 0) {
    hap_to_merge.push_back(hap_t_1);
  }
  while(hap_read(fp_hap2, &hap_t_2) == 0) {
    hap_to_merge.push_back(hap_t_2);
  }

  //sort
  sort(hap_to_merge.begin(), hap_to_merge.end(), comp_hap);

  int i = 1;
  vector<hap_pos_t> cpg_pos_a;
  vector<hap_pos_t> cpg_pos_b;
  hap_t hap_a = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};
  hap_t hap_b = {HAP_NULL_STRING, 0, 0, HAP_NULL_STRING, 0, HAP_DEFAULT_DIRECTION};
  hap_t merge_result;
  if (hap_to_merge.size() == 0) {
    hts_log_error("hap files are empty");
  } else if (hap_to_merge.size() == 1) {
    merge_result_vec.push_back(hap_to_merge[0]);
  } else if (hap_to_merge.size() >=2) {

    hap_b = hap_to_merge[0];
    int overlap_beg_a, overlap_beg_b, overlap_end_a, overlap_end_b;
    while(i < hap_to_merge.size()) {
      hap_a = hap_b;
      hap_b = hap_to_merge[i];
      cout << "start" << endl;
      hap_a.print();
      hap_b.print();
      cout << "end print" << endl;

      cout << "get cpg" << endl;
      cpg_pos_a.clear();
      cpg_pos_b.clear();
      cpg_pos_a = get_cpg(ctx_merge, hap_a);
      cpg_pos_b = get_cpg(ctx_merge, hap_b);
      cout << "end cpg" << endl;

      if (is_identity(hap_a, hap_b)) {
        //identity
        cout << "identity" << endl;
        hap_b.hap_count += hap_a.hap_count;

      } else if (is_overlap(cpg_pos_a, cpg_pos_b, hap_a, hap_b,
          &overlap_beg_a, &overlap_beg_b, &overlap_end_a, &overlap_end_b)) {
        //overlap
        cout << "overlap" << endl;
        merge_result = merge(hap_a, hap_b, overlap_beg_a, overlap_beg_b, overlap_end_a, overlap_end_b);
        hap_b = merge_result;
      } else {
        cout << "others" << endl;
        merge_result_vec.push_back(hap_a);
      }
      i++;
    }
    cout << "saving" << endl;
    merge_result_vec.push_back(hap_b);
  }

for(auto h : merge_result_vec) {
  h.print();
}


  return 0;
}
}//namespace std
