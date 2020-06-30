//
// Created by Yuhao Dan on 2020/4/13.
//
#include "convert.h"
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/regidx.h>
#include <chrono>
#include <algorithm>

namespace std {

bool Context::parse_region() {
  //parse the -r chr:start-end
  const char *reg = region;

  int flags = 0;

  while(*reg) {
    reg = sam_parse_region(hdr_bam, reg, &i_tid, &i_beg, &i_end, flags);
    if (!reg) {
      fprintf(stderr, "Failed to parse region\n");
      return false;
    }
  }
  i_chr = hdr_bam->target_name[i_tid];

  return true;
}

void Context::print_region() {
  cout << "Region parsed result: " << hdr_bam->target_name[i_tid] << ":" << i_beg << "-" << i_end << endl;
}

static const char *opt_string = "i:a:b:c:r:o:";

static const struct option long_opts[] = {
    { "input", required_argument, NULL, 'i' },
    { "aligner", optional_argument, NULL, 'a' },
    { "bed_file", optional_argument, NULL, 'b' },
    { "cpg_path", required_argument, NULL, 'c' },
    { "region", optional_argument, NULL, 'r' },
    { "output", optional_argument, NULL, 'o' },
    { NULL, no_argument, NULL, 0 }
};

bool load_get_cpg_with_idx(Context &ctx, char *chr, uint32_t beg, uint32_t end, uint32_t shift = 500) {
  //concat name of the tbi file
  ctx.cpg_pos.clear();

  uint32_t  i_beg;

  i_beg = beg - shift;

  if (i_beg < 0) {
    i_beg = 0;
  }

  uint32_t i_end = end + shift;

  int tbx_tid = tbx_name2id(ctx.idx_cpg, chr);
  ctx.cpg_itr = tbx_itr_queryi(ctx.idx_cpg, tbx_tid, i_beg, i_end);
  kstring_t ksbuf = {0, 0, NULL};
  string cpg_line_sv;//todo:string_view
  u_int32_t cpg_start = 0;
  u_int32_t cpg_end = 0;


  while(tbx_itr_next(ctx.fp_cpg, ctx.idx_cpg, ctx.cpg_itr, &ksbuf) >= 0) {
    cpg_line_sv = string(ksbuf.s); //todo:string_view
    int i = 0;
    for (; i < cpg_line_sv.size(); i++) {
      if (cpg_line_sv[i] == '\t') {
        cpg_line_sv = cpg_line_sv.substr(i + 1);
        //curser清零
        i = 0;
        //找到原line中的第2个tab位置
        for (; i < cpg_line_sv.size(); i++) {
          if (cpg_line_sv[i] == '\t') {
            cpg_start = atoi(string(cpg_line_sv.substr(0, i)).c_str());
            cpg_end = atoi(string(cpg_line_sv.substr(i, cpg_line_sv.size() - i)).c_str());
            //确保cpg位点在用户指定的范围内。
            ctx.cpg_pos.push_back(cpg_start);
          }
        }
      }
    }
  }
  return true;
}

inline SamRead::~SamRead() {
}

bool SamRead::init(Context &ctx) {
  bool ret;

  this->ctx = &ctx;

  read_len = ctx.aln->core.l_qseq; //length of the read.

  read_chr = ctx.hdr_bam->target_name[ctx.aln->core.tid] ; //checked
  read_start = ctx.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
  read_end = ctx.aln->core.pos + read_len;

  read_qual = bam_get_qual(ctx.aln); //quality string
  read_name = bam_get_qname(ctx.aln);
  flag = ctx.aln->core.flag; //

  read_map_quality = int(ctx.aln->core.qual); //mapping quality
  read_cigar = bam_get_cigar(ctx.aln);
  //get the seq pointer and put it into a vector
  uint8_t *seq_p = bam_get_seq(ctx.aln);
  if (read_len == 0) {
    return false;
  }

  if (strcmp(ctx.aligner, "BISMARK") == 0) {
    ret = _get_bismark_std();

    if (!ret) {
      hts_log_trace("_get_bismark_std(): fail to get bismark std.");
      return false;
    }

    if (_get_XM_tag(ctx)) {
      ret = _get_bismark_QC(ctx);
      if (!ret) {
        hts_log_trace("_get_bismark_QC(): fail to get bismark QC.");
        return false;
      }
    } else {
      hts_log_trace("XM tag is required in SAM");
      return false;
    }
  } else if(strcmp(ctx.aligner, "BSMAP") == 0) {
    if (!_get_ZS_tag(ctx)) {
      hts_log_trace("_get_ZS_tag(): fail to get ZS tag");
      return false;
    }

  } else if(strcmp(ctx.aligner, "UNKNOWN") == 0) {
    //TODO(butyuhao@foxmail.com) 使用unknown的时候，与BISMARK比较可能多出一些结果，需要排查一下原因。
    read_WC = DIRECTION_UNKNOWN;
  } else {
    hts_log_error("Only BSMAP, BISMARK and UNKNOWN are supported.");
    return false;
  }

  seq = new char[read_len];
  for(int i = 0; i < read_len; i++) {
    *(seq + i) = kbase[bam_seqi(seq_p, i)];
  }

  return true;
}

bool SamRead::haplo_type() {
  uint32_t r_pos;
  vector<uint32_t> cpg;
  vector<int8_t> hap_qual;
  string hap_seq = "";
  uint32_t pos;
  for (int i = 0; i < ctx->cpg_pos.size(); i++) {
    pos = ctx->cpg_pos[i];
    if (pos < read_start) {
      continue;
    }
    if (pos > read_end) {
      //因为cpg位点是有序的，所以超过end就说明之后不可能再有在范围内的位点了。
      break;
    }
    if (read_WC == DIRECTION_PLUS || read_WC == DIRECTION_UNKNOWN) {
      r_pos = pos - read_start;
    } else {
      r_pos = pos - read_start + 1;
    }
    if (r_pos >= read_len || r_pos < 0) {
      continue;
    }
    cpg.push_back(pos);

    hap_seq += seq[r_pos];

    hap_qual.push_back(read_qual[r_pos]);

  }

  if (seq != NULL) {
    delete [] seq;
    seq = NULL;
  }

  if (cpg.size() == 0) {
    QC = false;
  }
  _hap_seq = hap_seq;
  _hap_qual = hap_qual;
  _cpg = cpg;
  string hap_met = "";
  for (int i = 0; i < hap_seq.size(); i++) {
    char nucleobases = hap_seq[i];
    if (read_WC == DIRECTION_PLUS || read_WC == DIRECTION_UNKNOWN) {
      if (nucleobases == 'C') {
        hap_met += '1';
      } else if (nucleobases == 'T') {
        hap_met += '0';
      } else {
        hap_met += nucleobases;
        //cout << hap_met << endl;
        QC = false;
      }
    } else if (read_WC == DIRECTION_MINUS) {
      if (nucleobases == 'G') {
        hap_met += '1';
      } else if (nucleobases == 'A') {
        hap_met += '0';
      } else {
        hap_met += nucleobases;
        QC = false;
      }
    } else {
      cout << "Error! Strand undefined" << endl;
    }
  }
  _hap_met = hap_met;
  if (_cpg.size() > 0) {
    HT = HT_s(read_chr, _cpg[0], _cpg[_cpg.size() - 1], _hap_met, 1, read_WC);
  }
  return true;
}

bool SamRead::_get_bismark_std() {
  if (ctx->aln->core.flag == 99 || ctx->aln->core.flag == 147) {
    read_WC = DIRECTION_PLUS;
  } else if (ctx->aln->core.flag == 83 || ctx->aln->core.flag == 163) {
    read_WC = DIRECTION_MINUS;
  } else {
    cout << "sam_read::_get_bismark_std() Unknown flag" << endl;
    return false;
  }
  return true;
}

bool SamRead::_get_XM_tag(Context &ctx) {

  ctx.bam_aux_p = bam_aux_get(ctx.aln, "XM");
  if (!ctx.bam_aux_p) {
    return false;
  }
  return true;
}

bool SamRead::_get_ZS_tag(Context &ctx) {

  ctx.bam_aux_p = bam_aux_get(ctx.aln, "ZS");
  if (!ctx.bam_aux_p) {
    return false;
  }
  ZS_tag = bam_aux2Z(ctx.bam_aux_p);
  if (!ZS_tag) {
    hts_log_error("has no ZS tag");
    return false;
  }
  if (*ZS_tag == '+') {
    read_WC = DIRECTION_PLUS;
  } else if (*ZS_tag == '-'){
    read_WC = DIRECTION_MINUS;
  } else {
    hts_log_trace("Direction tag error");
    return false;
  }
  return true;
}

bool SamRead::_get_bismark_QC(Context &ctx) {
  XM_tag = bam_aux2Z(ctx.bam_aux_p);
  if (!XM_tag) {
    hts_log_error("has no XM tag");
    return false;
  }
  string XM_tag_sv  = string(XM_tag); //todo:string_view
  for (auto c : XM_tag_sv) {
    if (c == 'X' || c == 'H' || c == 'U') {
      hts_log_trace("XM tag QC check failed");
      return false;
    }
  }
  return true;
}

bool paired_end_check(SamRead &samF, SamRead &samR) {
  if (strcmp(samF.read_name, samR.read_name) != 0) {
    return false;
  }
  if (strcmp(samF.read_chr, samR.read_chr) != 0) {
    return false;
  }
  bool checkF = samF._cpg[samF._cpg.size() - 1] >= samR._cpg[0];
  bool checkR = samR._cpg[samR._cpg.size() - 1] >= samF._cpg[0];
  return (checkF && checkR) || (!checkF && !checkR);
}

HT_s paired_end_merge(SamRead &samF, SamRead &samR) {
  map<uint32_t, char> merged_seq;
  map<uint32_t, int8_t> merged_qual;
  map<uint32_t, char> merged_met;
  map<uint32_t, char>::iterator merged_met_itor;
  for (int i = 0; i < samF._cpg.size(); i++) {
    uint32_t pos = samF._cpg[i];
    merged_seq[pos] = samF._hap_seq[i];
    merged_qual[pos] = samF._hap_qual[i];
    merged_met[pos] = samF._hap_met[i];
  }
  for (int i = 0; i < samR._cpg.size(); i++) {
    uint32_t pos = samR._cpg[i];
    merged_met_itor = merged_met.find(pos);
    if (merged_met_itor == merged_met.end() || samR._hap_qual[i] > merged_qual[pos]) {
      merged_seq[pos] = samR._hap_seq[i];
      merged_qual[pos] = samR._hap_qual[i];
      merged_met[pos] = samR._hap_met[i];
    }
  }
  //map创建后，其key直接就是降序排列的
  string hap_seq = "";
  string hap_met = "";
  merged_met_itor = merged_met.begin();
  vector<uint32_t> cpg;
  while(merged_met_itor != merged_met.end()) {
    cpg.push_back(merged_met_itor->first);
    hap_seq += merged_seq[merged_met_itor->first];
    hap_met += merged_met[merged_met_itor->first];
    merged_met_itor++;
  }
  HT_s merged_HT = HT_s(samF.read_chr, cpg[0], cpg[cpg.size() - 1], hap_met, 1, samF.read_WC);
  return merged_HT;
}

//struct cmp
//{
//  bool operator()(const pair<string,u_int32_t> &p1,const pair<string,u_int32_t> &p2)
//  {
//    return p1.second < p2.second;
//  }
//};

bool region_to_parse(Context &ctx) {
  /* Recognize the way the user specifies the region
   * 1.use command line -r option to specify single region --> SINGLE_REGION
   * 2.use bed file --> MULTI_REGION
   * 3.not specify, process whole file --> WHOLE_FILE
   * */
  if (ctx.region) {
    ctx.region_to_parse = SINGLE_REGION;
    hts_log_trace("SINGLE_REGION");
  } else if (ctx.fn_bed) {
    ctx.region_to_parse = MULTI_REGION;
    hts_log_trace("MULTI_REGION");
  } else if (ctx.fn_cpg && ctx.fn_bam) {
    ctx.region_to_parse = WHOLE_FILE;
    hts_log_trace("WHOLE_FILE");
  } else {
    return false;
  }
  return true;
}

bool comp(const HT_s &a, const HT_s &b)
{
  if (strcmp(a.h_chr, b.h_chr) == 0) {
    if (a.h_start != b.h_start) {
      return a.h_start < b.h_start;
    } else if (a.h_end != b.h_end) {
      return a.h_end < b.h_end;
    } else if (a.hap_met != b.hap_met) {
      return a.hap_met < b.hap_met;
    } else if (a.ht_count != b.ht_count){
      return a.ht_count < b.ht_count;
    } else {
      return a.WC < b.WC;
    }
  } else {
    return strcmp(a.h_chr, b.h_chr) < 0;
  }
}

bool load_cpg_no_idx(Context &ctx) {
  kstring_t cpg_line = {0,0,NULL};
  unordered_map<int, vector<hts_pos_t>>::iterator cpg_pos_map_itor;
  while (hts_getline(ctx.fp_cpg, KS_SEP_LINE, &cpg_line) > 0) {

    char *p ,*q;
    int tid;
    hts_pos_t cpg_start;
    p = q = cpg_line.s;
    while(*q && *q != '\t') {
      q++;
    }
    *q = '\0';
    tid = bam_name2id(ctx.hdr_bam, p);
    *q = '\t';
    p = q + 1;
    q = p;
    while(*q && *q != '\t') {
      q++;
    }
    *q = '\0';
    cpg_start = atoll(p);

    cpg_pos_map_itor = ctx.cpg_pos_map.find(tid);

    if (cpg_pos_map_itor == ctx.cpg_pos_map.end()) {
      vector<hts_pos_t> v;
      v.push_back(cpg_start);
      ctx.cpg_pos_map[tid] = v;
    } else {
      ctx.cpg_pos_map[tid].push_back(cpg_start);
    }

  }
  return true;
}

int _lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos)//二分查找求下界
//找到cpg_pos这个数字在v中的下标

{
  int low = 0, high = v.size() - 1;
  while(low < high)
  {
    int mid = low + (high - low)/2;
    if(v[mid] >= cpg_pos) high = mid;
    else low = mid + 1;
  }
  return low;
}

bool get_cpg_no_idx(Context &ctx, char *chr, hts_pos_t &beg, hts_pos_t &end, hts_pos_t shift = 500) {
  ctx.cpg_pos.clear();

  hts_pos_t  i_beg;

  i_beg = beg - shift;

  if (i_beg < 0) {
    i_beg = 0;
  }

  hts_pos_t i_end = end + shift;

  int tid = bam_name2id(ctx.hdr_bam, chr);

  int pos = _lower_bound(ctx.cpg_pos_map[tid], i_beg);

  unordered_map<int, vector<hts_pos_t>>::iterator cpg_pos_map_itor;

  cpg_pos_map_itor = ctx.cpg_pos_map.find(tid);

  if (cpg_pos_map_itor != ctx.cpg_pos_map.end()) {
    while (pos < ctx.cpg_pos_map[tid].size()) {
      if (ctx.cpg_pos_map[tid][pos] >= i_beg && ctx.cpg_pos_map[tid][pos] <= i_end) {
        ctx.cpg_pos.push_back(ctx.cpg_pos_map[tid][pos]);
        pos++;
      } else {
        break;
      }
    }
  }
}

vector<HT_s> itor_sam(Context &ctx) {

  map<string, vector<SamRead>> sam_map;
  map<string, vector<SamRead>>::iterator iter;
  vector<HT_s> HT_vec;

  if (ctx.region_to_parse == SINGLE_REGION) {
    //load tbi index outside the load_get_cpg_with_idx()
    string cpg_idx_fn = string(ctx.fn_cpg) + string(".tbi");
    ctx.idx_cpg = tbx_index_load(cpg_idx_fn.c_str());
    ctx.has_idx_cpg = true;
  }


  if (ctx.region_to_parse == SINGLE_REGION) {
    auto single_start = std::chrono::high_resolution_clock::now(); //stop_watch
    //load bai index
    static string bam_idx_fn = string(ctx.fn_bam) + ".bai";
    ctx.idx_bam = sam_index_load(ctx.fp_bam, bam_idx_fn.c_str());

    hts_itr_t *sam_itr = sam_itr_queryi(ctx.idx_bam, ctx.i_tid, ctx.i_beg, ctx.i_end);

    while(sam_itr_next(ctx.fp_bam, sam_itr, ctx.aln) >= 0) {

      if (ctx.aln->core.flag & BAM_FQCFAIL || ctx.aln->core.flag & BAM_FUNMAP || ctx.aln->core.flag & BAM_FDUP
          || ctx.aln->core.flag & BAM_FSECONDARY || ctx.aln->core.flag & BAM_FSUPPLEMENTARY) {
        continue;
      }

      SamRead sam_r = SamRead();

      int ret = sam_r.init(ctx);

      if (!ret) {
        hts_log_trace("");
        if (sam_r.seq != NULL) {
          delete [] sam_r.seq;
        }
        continue;
      }

      string qname = string(sam_r.read_name);

      load_get_cpg_with_idx(ctx, sam_r.read_chr, sam_r.read_start, sam_r.read_end);

      sam_r.haplo_type();

      if (!sam_r.QC) {
        continue;
      }

      iter = sam_map.find(qname);

      if (iter == sam_map.end()) {
        vector<SamRead> v;
        v.push_back(sam_r);
        sam_map[qname] = v;
      } else {
        vector<SamRead> v;
        v = sam_map[qname];
        v.push_back(sam_r);
        sam_map[qname] = v;
      }

    }

  } else if (ctx.region_to_parse == WHOLE_FILE) {
    load_cpg_no_idx(ctx);
    while(sam_read1(ctx.fp_bam, ctx.hdr_bam, ctx.aln) >= 0) {
      if (ctx.aln->core.flag & BAM_FQCFAIL || ctx.aln->core.flag & BAM_FUNMAP || ctx.aln->core.flag & BAM_FDUP
          || ctx.aln->core.flag & BAM_FSECONDARY || ctx.aln->core.flag & BAM_FSUPPLEMENTARY) {
        continue;
      }

      SamRead sam_r = SamRead();

      int ret = sam_r.init(ctx);

      if (!ret) {
        hts_log_trace("");
        continue;
      }

      string qname = string(sam_r.read_name);

      get_cpg_no_idx(ctx, sam_r.read_chr, sam_r.read_start, sam_r.read_end);

      sam_r.haplo_type();
      if (!sam_r.QC) {
        continue;
      }

      iter = sam_map.find(qname);
      if (iter == sam_map.end()) {
        vector<SamRead> v;
        v.push_back(sam_r);
        sam_map[qname] = v;
      } else {
        vector<SamRead> v;
        v = sam_map[qname];
        v.push_back(sam_r);
        sam_map[qname] = v;
      }
    }
  } else if (ctx.region_to_parse == MULTI_REGION) {
    load_cpg_no_idx(ctx);
    regidx_t *idx = regidx_init(ctx.fn_bed,NULL,NULL,0,NULL);
    regitr_t *itr = regitr_init(idx);

    //load bai index
    static string bam_idx_fn = string(ctx.fn_bam) + ".bai";
    ctx.idx_bam = sam_index_load(ctx.fp_bam, bam_idx_fn.c_str());

    while ( regitr_loop(itr) ) {
      int tid = bam_name2id(ctx.hdr_bam, itr->seq);

      if (tid == -1 || tid == -2) {
        hts_log_trace("tid not found");
        continue;
      }
      hts_itr_t *sam_itr = sam_itr_queryi(ctx.idx_bam, tid, itr->beg+1, itr->end+1);
      if (sam_itr == NULL) {
        hts_log_debug("%d:%lld-%lld", tid, itr->beg+1, itr->end+1);
        continue;
      }
      while(sam_itr_next(ctx.fp_bam, sam_itr, ctx.aln) >= 0) {
        if (ctx.aln->core.flag & BAM_FQCFAIL || ctx.aln->core.flag & BAM_FUNMAP || ctx.aln->core.flag & BAM_FDUP
            || ctx.aln->core.flag & BAM_FSECONDARY || ctx.aln->core.flag & BAM_FSUPPLEMENTARY) {
          continue;
        }

        SamRead sam_r = SamRead();

        int ret = sam_r.init(ctx);

        if (!ret) {
          hts_log_trace("");
          continue;
        }

        string qname = string(sam_r.read_name);

        get_cpg_no_idx(ctx, sam_r.read_chr, sam_r.read_start, sam_r.read_end);

        sam_r.haplo_type();
        if (!sam_r.QC) {
          continue;
        }

        iter = sam_map.find(qname);
        if (iter == sam_map.end()) {
          vector<SamRead> v;
          v.push_back(sam_r);
          sam_map[qname] = v;
        } else {
          vector<SamRead> v;
          v = sam_map[qname];
          v.push_back(sam_r);
          sam_map[qname] = v;
        }
      }
    }
    regidx_destroy(idx);
    regitr_destroy(itr);
  } else {
    hts_log_error("region_error");
    exit(1);
  }





  //merge
  for (auto sam_l :  sam_map) {
    if (sam_l.second.size() == 2) {
      SamRead samF = sam_l.second[0];
      SamRead samR = sam_l.second[1];
      if (paired_end_check(samF, samR)) {
        HT_s ht = paired_end_merge(samF, samR);
        HT_vec.push_back(ht);
      } else {
        for (int i = 0; i < sam_l.second.size(); i++) {
          HT_vec.push_back(sam_l.second[i].HT);
        }
      }
    } else {
      for (int i = 0; i < sam_l.second.size(); i++) {
        HT_vec.push_back(sam_l.second[i].HT);
      }
    }
  }

  vector<HT_s>::iterator ht_itor;
  for (ht_itor = HT_vec.begin(); ht_itor != HT_vec.end(); ht_itor++) {//auto _ht: HT_vec
    string ht_id = (*ht_itor).to_str();

    map<string, int>::iterator res_map_itor;

    res_map_itor = ctx.res_map.find(ht_id);
    if (res_map_itor == ctx.res_map.end()){
      ctx.res_map[ht_id] = 1;
    } else {
      ctx.res_map[ht_id] += 1;
    }
  }

  return HT_vec;
}

inline Context::~Context() {
  //to_do明确一下哪些指针需要被关掉。
  if (fp_bam) {
    hts_close(fp_bam);
  }
  if (fp_cpg) {
    hts_close(fp_cpg);
  }
  if (aln) {
    bam_destroy1(aln);
  }

  //当用到load_get_no_idx()的时候，生成以下两个指针
  if (has_idx_cpg) {
    tbx_destroy(idx_cpg);
  }
  if (has_idx_cpg) {
    hts_itr_destroy(cpg_itr);
  }

  if (idx_bam) {
    hts_idx_destroy(idx_bam);
  }
  if (sam_itr) {
    hts_itr_destroy(sam_itr);
  }
  if (hdr_bam) {
    bam_hdr_destroy(hdr_bam);
  }
}

bool open_bam_file(Context &ctx) {

  ctx.fp_bam = hts_open(ctx.fn_bam, "r");
  ctx.hdr_bam = sam_hdr_read(ctx.fp_bam); //read header

  if (ctx.hdr_bam == NULL) {
    return false;
  }
  ctx.aln = bam_init1();
  if (ctx.aln == NULL) {
    return false;
  }
  string fn_bam_idx = string(ctx.fn_bam) + ".bai";
  ctx.idx_bam = sam_index_load2(ctx.fp_bam, ctx.fn_bam, fn_bam_idx.c_str());

  return true;
}

bool open_cpg_file(Context &ctx) {
  //open the cpg file (.gz file)
  ctx.fp_cpg = hts_open(ctx.fn_cpg, "r");
  if (ctx.fp_cpg == NULL) {
    return false;
  }
  return true;
}

inline void HT_s::get_WC_symbol() {
  {
    if (WC == DIRECTION_PLUS) {
      WC_symbol = '+';
    } else if (WC == DIRECTION_MINUS) {
      WC_symbol = '-';
    } else if (WC == DIRECTION_UNKNOWN) {
      WC_symbol = '*';
    }
  }
}

int main_convert(int argc, char *argv[]) {
  hts_log_trace("enter main_convert");
//TODO(butyuhao@foxmail.com): 增加检查option合法性的部分
//TODO(butyuhao@foxmail.com): 增加日志

  hts_log_trace("create Context()");
  Context ctx = Context();

  hts_log_trace("parse options");
  int long_index;

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        ctx.fn_bam = optarg;
        break;
      }
      case 'a': {
        ctx.aligner = optarg;
        break;
      }
      case 'b': {
        ctx.fn_bed = optarg;
        break;
      }
      case 'c': {
        ctx.fn_cpg = optarg;
        break;
      }
      case 'r': {
        ctx.region = optarg;
        break;
      }
      case 'o': {
        ctx.output_path = optarg;
        break;
      }
      default: {
        break;
      }
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }

    int ret;

    hts_log_trace("open_bam_file(ctx)");
    ret = open_bam_file(ctx);
    if (!ret) {
      hts_log_error("open_bam_file():fail to open bam file");
      return EXIT_FAILURE;
    }
    hts_log_trace("region_to_parse(ctx)");
    ret = region_to_parse(ctx);
    if (!ret) {
      hts_log_error("region_to_parse()");
    }
    if (ctx.region_to_parse == SINGLE_REGION) {
      ret = ctx.parse_region();
      hts_log_info("parse_region(): %s:%lld-%lld", ctx.hdr_bam->target_name[ctx.i_tid], ctx.i_beg, ctx.i_end);
    }
    if (!ret) {
      hts_log_error("parse_region():fail to parse region");
      return EXIT_FAILURE;
    }

    hts_log_trace("open_cpg_file(ctx)");
    ret = open_cpg_file(ctx);
    if (!ret) {
      hts_log_error("open_cpg_file():fail to open cpg file");
      return EXIT_FAILURE;
    }

    vector<HT_s> HT_vec;
    hts_log_trace("itor_sam(ctx).");
    HT_vec = itor_sam(ctx);

    hts_log_trace("Output.");

    string out_stream_name;
    if (ctx.output_path) {
      out_stream_name = ctx.output_path;
    } else {
      out_stream_name = "out.hap";
    }
    ofstream out_stream(out_stream_name);

    unordered_map<string, bool> is_overlap;
    unordered_map<string, bool>::iterator itor;
    vector<HT_s>::iterator ht_itor;
    for (ht_itor = HT_vec.begin(); ht_itor != HT_vec.end(); ht_itor++) {
      (*ht_itor).ht_count = ctx.res_map[(*ht_itor).to_str()];
    }
  //sort
  sort(HT_vec.begin(), HT_vec.end(), comp);

  for (ht_itor = HT_vec.begin(); ht_itor != HT_vec.end(); ht_itor++) {
    (*ht_itor).get_WC_symbol();
    string line = string((*ht_itor).h_chr) + '\t' + to_string((*ht_itor).h_start) + '\t' +
        to_string((*ht_itor).h_end) + '\t' + (*ht_itor).hap_met + '\t' +
        to_string((*ht_itor).ht_count) + '\t' + (*ht_itor).WC_symbol;
    itor = is_overlap.find(line);
    if (itor == is_overlap.end()) {
      out_stream << line << '\n';
    }
    is_overlap[line] = true;
  }

    out_stream.close();

  return EXIT_SUCCESS;
}
} // namespace std
