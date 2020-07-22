#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <algorithm>
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
  vector<hap_pos_t>_cpg_pos(hap_read.hap_str.size(), 0);

  int pos = _lower_bound(ctx_merge.cpg_pos_map[hap_read.chr], hap_read.chr_beg);

  unordered_map<string, vector<hts_pos_t>>::iterator cpg_pos_map_itor;

  cpg_pos_map_itor = ctx_merge.cpg_pos_map.find(hap_read.chr);

  int i = 0;
  if (cpg_pos_map_itor != ctx_merge.cpg_pos_map.end()) {
    while (pos < ctx_merge.cpg_pos_map[hap_read.chr].size()) {

      if (ctx_merge.cpg_pos_map[hap_read.chr][pos] >= hap_read.chr_beg &&
      ctx_merge.cpg_pos_map[hap_read.chr][pos] <= hap_read.chr_end) {

        _cpg_pos[i] = ctx_merge.cpg_pos_map[hap_read.chr][pos];
        i++;
        pos++;
      } else {
        break;
      }
    }
  }

  return _cpg_pos;
}

bool is_overlap(vector<hap_pos_t> &cpg_pos_1, vector<hap_pos_t> &cpg_pos_2, int *i1, int *i2) {
  for(*i1 = 0; *i1 < cpg_pos_1.size(); (*i1)++) {
    for(*i2 = 0; *i2 < cpg_pos_2.size(); (*i2)++) {
      hap_pos_t c1 = cpg_pos_1[*i1];
      hap_pos_t c2 = cpg_pos_2[*i2];
      if (c1 == c2) {
        return true;
      }
    }
  }
  return false;
}

bool is_identity(hap_t &hap_a, hap_t &hap_b) {
  if (hap_a.chr == hap_b.chr && hap_a.chr_beg == hap_b.chr_beg &&
      hap_a.chr_end == hap_b.chr_end && hap_a.hap_str == hap_b.hap_str &&
      hap_a.hap_direction == hap_b.hap_direction) {
    return true;
  }
  return false;
}

hap_t merge(hap_t &hap_t_1, hap_t &hap_t_2, int &i1, int &i2) {
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
  } else {
    hap_merge.chr_beg = hap_t_2.chr_beg;
  }
  if (hap_merge.chr_beg == hap_t_1.chr_beg) {
    hap_str += hap_t_1.hap_str.substr(0, i1);
  } else {
    hap_str += hap_t_2.hap_str.substr(0, i2);
  }

  int shift = 0 ;

  while(i1 < hap_t_1.hap_str.size() && i2 < hap_t_2.hap_str.size()  &&
  hap_t_1.hap_str[i1 + shift] ==  hap_t_2.hap_str[i2 + shift]) {
    shift++;
  }

  hap_str += hap_t_1.hap_str.substr(i1, shift);

  if(i1 + shift == hap_t_1.hap_str.size()) {
    hap_str += hap_t_2.hap_str.substr(i2 + shift, hap_t_2.hap_str.size() - (i2 + shift));
    hap_merge.chr_end = hap_t_2.chr_end;
  } else if(i2 + shift == hap_t_2.hap_str.size()) {
    hap_str += hap_t_1.hap_str.substr(i1 + shift, hap_t_1.hap_str.size() - (i1 + shift));
    hap_merge.chr_end = hap_t_1.chr_end;
  } else {
    //重叠部分不一致，保留两个hap read,使用chr_beg=-1 与 chr_end=-1来传递此信息
    hap_merge = {HAP_NULL_STRING, -1, -1, HAP_NULL_STRING, HAP_DEFAULT_DIRECTION};
  }

  hap_merge.hap_str = hap_str;
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

    //test

    //get_cpg_pos

  vector<hap_pos_t> mock_cpg_pos_1;
  mock_cpg_pos_1.push_back(10234);
  mock_cpg_pos_1.push_back(10237);
  mock_cpg_pos_1.push_back(10345);
  mock_cpg_pos_1.push_back(10378);
  mock_cpg_pos_1.push_back(12543);
  mock_cpg_pos_1.push_back(12754);
  mock_cpg_pos_1.push_back(13655);
  vector<hap_pos_t> mock_cpg_pos_2;
  mock_cpg_pos_2.push_back(10378);
  hap_t mock_hap_1;
  mock_hap_1.chr = string("1");
  mock_hap_1.chr_beg = 10234;
  mock_hap_1.chr_end = 13655;

  mock_hap_1.hap_str = string("1101001");

  hap_t mock_hap_2;
  mock_hap_2.chr = string("1");
  mock_hap_2.chr_beg = 10378;
  mock_hap_2.chr_end = 10378;
  mock_hap_2.hap_str = string("1");

  hap_t merge_result;

  int i1, i2;

  if (is_overlap(mock_cpg_pos_1, mock_cpg_pos_2, &i1, &i2)) {
    merge_result = merge(mock_hap_1, mock_hap_2, i1, i2);
  }
  while(1);

  //test

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
  if (hap_to_merge.size() == 0) {
    hts_log_error("hap files are empty");
  } else if (hap_to_merge.size() == 1) {
    merge_result_vec.push_back(hap_to_merge[0]);
  } else if (hap_to_merge.size() >=2) {

    hap_b = hap_to_merge[0];
    int pos_a, pos_b;
    while(i < hap_to_merge.size()) {
      hap_a = hap_b;
      hap_b = hap_to_merge[i];
      cout << "start" << endl;

      hap_a.print();
      hap_b.print();

      cpg_pos_a = get_cpg(ctx_merge, hap_a);
      cpg_pos_b = get_cpg(ctx_merge, hap_b);
      if (is_identity(hap_a, hap_b)) {
        //identity
        cout << "identity" << endl;
        hap_b.hap_count += hap_a.hap_count;

      } else if (hap_a.chr == hap_b.chr && hap_a.hap_direction == hap_b.hap_direction &&
                is_overlap(cpg_pos_a, cpg_pos_b, &pos_a, &pos_b)) {
        //overlap
        cout << "overlap" << endl;
        hap_t merge_result;
        merge_result = merge(hap_a, hap_b, pos_a, pos_b);
        if (merge_result.chr_beg == -1 && merge_result.chr_end == -1) {
          cout << "overlap_branch1" << endl;
          hap_a.print();
          merge_result_vec.push_back(hap_a);
        }
        else {
          hap_b = merge_result; //-->问题在这边
        }
      } else {
        cout << "other" << endl;
        hap_a.print();
        merge_result_vec.push_back(hap_a);
      }
      i++;
    }
  }

for(auto h : merge_result_vec) {
  h.print();
}


  return 0;
}
}//namespace std
