//
// Created by Yuhao Dan on 2020/4/13.
//

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <fcntl.h>
#include <stdio.h>
#include "./htslib-1.10.2/htslib/kseq.h"
#include "./htslib-1.10.2/htslib/bgzf.h"
#include "./htslib-1.10.2/htslib/hfile.h"
#include "./htslib-1.10.2/htslib/regidx.h"
#include "./include/convert.h"
#include "./include/utils.h"
namespace std {

bool ContextConvert::parse_region() {
  //parse the -r chr:start-end
  const char *reg = region;

  int flags = 0;

  while(*reg) {
    reg = sam_parse_region(hdr_bam, reg, &i_tid, &i_beg, &i_end, flags);
    if (!reg) {
      return false;
    }
  }
  i_chr = hdr_bam->target_name[i_tid];

  return true;
}

bool load_get_cpg_with_idx(ContextConvert &ctx, char *chr, hts_pos_t beg, hts_pos_t end, hts_pos_t shift = 500) {
  //concat name of the tbi file
  ctx.cpg_pos.clear();

  hts_pos_t  i_beg;

  i_beg = beg - shift;

  if (i_beg < 0) {
    i_beg = 0;
  }

  hts_pos_t i_end = end + shift;

  int tbx_tid = tbx_name2id(ctx.idx_cpg, chr);
  ctx.cpg_itr = tbx_itr_queryi(ctx.idx_cpg, tbx_tid, i_beg, i_end);
  kstring_t ksbuf = {0, 0, NULL};
  string cpg_line_sv;
  hts_pos_t cpg_start = 0;
  hts_pos_t cpg_end = 0;


  while(tbx_itr_next(ctx.fp_cpg, ctx.idx_cpg, ctx.cpg_itr, &ksbuf) >= 0) {
    cpg_line_sv = string(ksbuf.s);
    int i = 0;
    for (; i < cpg_line_sv.size(); i++) {
      if (cpg_line_sv[i] == '\t') {
        cpg_line_sv = cpg_line_sv.substr(i + 1);
        //curser清零
        i = 0;
        //找到原line中的第2个tab位置
        for (; i < cpg_line_sv.size(); i++) {
          if (cpg_line_sv[i] == '\t') {
            cpg_start = atoll(string(cpg_line_sv.substr(0, i)).c_str());
            cpg_end = atoll(string(cpg_line_sv.substr(i, cpg_line_sv.size() - i)).c_str());
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

bool SamRead::init(ContextConvert &ctx) {
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
  //quality check
  if (_get_XM(ctx)) {
    ret = _get_bismark_QC(ctx);
    if (!ret) {
      return false;
    }
  }
  //get direction
  if (ctx.non_directional) {
    read_WC = DIRECTION_UNKNOWN;
  } else {
    if (ctx.aln->core.flag & BAM_FPAIRED && ctx.aln->core.flag & BAM_FPROPER_PAIR) {
      if (ctx.aln->core.flag & BAM_FREAD1) {
        if (ctx.aln->core.flag & BAM_FREVERSE) {
          read_WC = DIRECTION_MINUS;
        } else if (ctx.aln->core.flag & BAM_FMREVERSE) {
          read_WC = DIRECTION_PLUS;
        }
      } else if (ctx.aln->core.flag & BAM_FREAD2 ) {
        if (ctx.aln->core.flag & BAM_FREVERSE) {
          read_WC = DIRECTION_PLUS;
        } else if (ctx.aln->core.flag & BAM_FMREVERSE) {
          read_WC = DIRECTION_MINUS;
        }
      }
    } else if (!(ctx.aln->core.flag & BAM_FPAIRED) && !(ctx.aln->core.flag & BAM_FPROPER_PAIR)){
      if (ctx.aln->core.flag & BAM_FREVERSE) {
        read_WC = DIRECTION_MINUS;
      } else {
        read_WC = DIRECTION_PLUS;
      }
    } else {
      return false;
    }
  }

  seq = new char[read_len];
  for(int i = 0; i < read_len; i++) {
    *(seq + i) = kbase[bam_seqi(seq_p, i)];
  }

  return true;
}

bool SamRead::haplo_type() {
  hts_pos_t r_pos;
  vector<hts_pos_t> cpg;
  vector<int8_t> hap_qual;
  _hap_seq = "";
  hts_pos_t pos;
  for (int i = 0; i < ctx->cpg_pos.size(); i++) {
    pos = ctx->cpg_pos[i];

    if (read_WC == DIRECTION_PLUS || read_WC == DIRECTION_UNKNOWN) {
      if (pos < read_start) {
        continue;
      }
      if (pos > read_end) {
        break;
      }
      r_pos = pos - read_start;
    } else {
      if (pos < read_start - 1) {
        continue;
      }
      if (pos > read_end - 1) {
        break;
      }
      r_pos = pos - read_start + 1;
    }

    if (r_pos >= read_len || r_pos < 0) {
      continue;
    }
    cpg.push_back(pos);

    _hap_seq += seq[r_pos];

    hap_qual.push_back(read_qual[r_pos]);

  }

  if (seq != NULL) {
    delete [] seq;
    seq = NULL;
  }

  if (cpg.size() == 0) {
    hts_log_trace("remove this read due to the size of cpg is 0.");
    QC = false;
  }
  _hap_qual = hap_qual;
  _cpg = cpg;
  _hap_met = "";
  for (int i = 0; i < _hap_seq.size(); i++) {
    char nucleobases = _hap_seq[i];
    if (read_WC == DIRECTION_PLUS || read_WC == DIRECTION_UNKNOWN) {
      if (nucleobases == 'C') {
        _hap_met += '1';
      } else if (nucleobases == 'T') {
        _hap_met += '0';
      } else {
        _hap_met += nucleobases;
        hts_log_trace("Direction + or * : beg: %lld, end: %lld, nucleobases error: %s", read_start,read_end, _hap_met.c_str());
        QC = false;
      }
    } else if (read_WC == DIRECTION_MINUS) {
      if (nucleobases == 'G') {
        _hap_met += '1';
      } else if (nucleobases == 'A') {
        _hap_met += '0';
      } else {
        _hap_met += nucleobases;
        hts_log_trace("Direction - : beg: %lld, end: %lld, nucleobases error: %s", read_start,read_end, _hap_met.c_str());
        QC = false;
      }
    } else {
      hts_log_trace("Strand undefined");
      QC = false;
    }
  }
  if (_cpg.size() > 0) {
    HT = HT_s(read_chr, _cpg[0], _cpg[_cpg.size() - 1], _hap_met, 1, read_WC);
  }
  return true;
}

bool convert_opt_check(ContextConvert &ctx) {
  if (ctx.fn_bam == NULL || ctx.fn_cpg == NULL) {
    hts_log_error("Please specify -i and -c options.");
    return false;
  }
  if (ctx.region != NULL && ctx.fn_bed != NULL) {
    hts_log_error("You can not specify both -r and -b options.");
    return false;
  }
  return true;
}

bool SamRead::_get_XM(ContextConvert &ctx) {

  ctx.bam_aux_p = bam_aux_get(ctx.aln, "XM");
  if (!ctx.bam_aux_p) {
    return false;
  }
  return true;
}

bool SamRead::_get_bismark_QC(ContextConvert &ctx) {
  XM_string = bam_aux2Z(ctx.bam_aux_p);
  if (!XM_string) {
    hts_log_error("has no XM string");
    return false;
  }
  string XM_str  = string(XM_string); //todo:string_view
  for (auto c : XM_str) {
    if (c == 'X' || c == 'H' || c == 'U') {
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
  map<hts_pos_t, char> merged_seq;
  map<hts_pos_t, int8_t> merged_qual;
  map<hts_pos_t, char> merged_met;
  map<hts_pos_t, char>::iterator merged_met_itor;
  for (int i = 0; i < samF._cpg.size(); i++) {
    hts_pos_t pos = samF._cpg[i];
    merged_seq[pos] = samF._hap_seq[i];
    merged_qual[pos] = samF._hap_qual[i];
    merged_met[pos] = samF._hap_met[i];
  }
  for (int i = 0; i < samR._cpg.size(); i++) {
    hts_pos_t pos = samR._cpg[i];
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
  vector<hts_pos_t> cpg;
  while(merged_met_itor != merged_met.end()) {
    cpg.push_back(merged_met_itor->first);
    hap_seq += merged_seq[merged_met_itor->first];
    hap_met += merged_met[merged_met_itor->first];
    merged_met_itor++;
  }
  HT_s merged_HT = HT_s(samF.read_chr, cpg[0], cpg[cpg.size() - 1], hap_met, 1, samF.read_WC);
  return merged_HT;
}

bool region_to_parse(ContextConvert &ctx) {
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

bool comp_HT_vec(const HT_s &a, const HT_s &b)
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

bool load_cpg_no_idx(ContextConvert &ctx) {
  kstring_t cpg_line = {0,0,NULL};
  unordered_map<int, vector<hts_pos_t> >::iterator cpg_pos_map_itor;
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

void get_cpg_no_idx(ContextConvert &ctx, char *chr, hts_pos_t &beg, hts_pos_t &end, hts_pos_t shift = 500) {
  ctx.cpg_pos.clear();

  hts_pos_t  i_beg;

  i_beg = beg - shift;

  if (i_beg < 0) {
    i_beg = 0;
  }

  hts_pos_t i_end = end + shift;

  int tid = bam_name2id(ctx.hdr_bam, chr);

  int pos = _lower_bound(ctx.cpg_pos_map[tid], i_beg);

  unordered_map<int, vector<hts_pos_t> >::iterator cpg_pos_map_itor;

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

vector<HT_s> itor_sam(ContextConvert &ctx) {

  map<string, vector<SamRead> > sam_map;
  map<string, vector<SamRead> >::iterator iter;
  vector<HT_s> HT_vec;

  if (ctx.region_to_parse == SINGLE_REGION) {
    //load tbi index outside the load_get_cpg_with_idx()
    string cpg_idx_fn = string(ctx.fn_cpg) + string(".tbi");
    ctx.idx_cpg = tbx_index_load(cpg_idx_fn.c_str());
    ctx.has_idx_cpg = true;

    if (ctx.idx_cpg == NULL) {
      string error_message = "Please run the command to generate the index file (.tbi):\ntabix -b 2 -e 3 -p bed " +
          string(ctx.fn_cpg) + "\nand place the index file with the CpG file in the same folder.";
      hts_log_error("%s", error_message.c_str());
      vector<HT_s> null_HT_s;
      return null_HT_s;
    }
  }

  uint64_t cnt = 0;

  if (ctx.region_to_parse == SINGLE_REGION) {
    auto single_start = std::chrono::high_resolution_clock::now(); //stop_watch
    //load bai index
    static string bam_idx_fn = string(ctx.fn_bam) + ".bai";
    ctx.idx_bam = sam_index_load(ctx.fp_bam, bam_idx_fn.c_str());

    if (ctx.idx_bam == NULL) {
      string error_message = "Please run the command to generate the index file (.bai):\nsamtools index " +
      string(ctx.fn_bam) + "\nand place the index file with the input file in the same folder.";
      hts_log_error("%s", error_message.c_str());
      vector<HT_s> null_HT_s;
      return null_HT_s;
    }

    hts_itr_t *sam_itr = sam_itr_queryi(ctx.idx_bam, ctx.i_tid, ctx.i_beg, ctx.i_end);

    while(sam_itr_next(ctx.fp_bam, sam_itr, ctx.aln) >= 0) {
      ++cnt;
      if (cnt % 1000000 == 0) {
        cout << cnt << " reads processed." << endl;
      }

      if (ctx.aln->core.flag & BAM_FQCFAIL || ctx.aln->core.flag & BAM_FUNMAP || ctx.aln->core.flag & BAM_FDUP
          || ctx.aln->core.flag & BAM_FSECONDARY || ctx.aln->core.flag & BAM_FSUPPLEMENTARY) {
        continue;
      }

      SamRead sam_r = SamRead();

      int ret = sam_r.init(ctx);

      if (!ret) {
        hts_log_trace("san_r init error.");
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

      ++cnt;
      if (cnt % 1000000 == 0) {
        cout << cnt << " reads processed." << endl;
      }

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
        hts_log_trace("--> QC check fail, remove the read.");
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

    if (ctx.idx_bam == NULL) {
      string error_message = "Please run the command to generate the index file (.bai):\nsamtools index " +
          string(ctx.fn_bam) + "\nand place the index file with the input file in the same folder.";
      hts_log_error("%s", error_message.c_str());
      vector<HT_s> null_HT_s;
      return null_HT_s;
    }

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

        ++cnt;
        if (cnt % 1000000 == 0) {
          cout << cnt << " reads processed." << endl;
        }

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
    idx = NULL;
    regitr_destroy(itr);
    itr = NULL;
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

inline ContextConvert::~ContextConvert() {
  //to_do明确一下哪些指针需要被关掉。
  if (fp_bam) {
    hts_close(fp_bam);
    fp_bam = NULL;
  }
  if (fp_cpg) {
    hts_close(fp_cpg);
    fp_cpg = NULL;
  }
  if (aln) {
    bam_destroy1(aln);
    aln = NULL;
  }

  //当用到load_get_no_idx()的时候，生成以下两个指针
  if (idx_cpg) {
    tbx_destroy(idx_cpg);
    idx_cpg = NULL;
  }
  if (cpg_itr) {
    hts_itr_destroy(cpg_itr);
    cpg_itr = NULL;
  }

  if (idx_bam) {
    hts_idx_destroy(idx_bam);
    idx_bam = NULL;
  }
  if (sam_itr) {
    hts_itr_destroy(sam_itr);
    sam_itr = NULL;
  }
  if (hdr_bam) {
    bam_hdr_destroy(hdr_bam);
    hdr_bam = NULL;
  }
}

bool open_bam_file(ContextConvert &ctx) {

  ctx.fp_bam = hts_open(ctx.fn_bam, "r");
  ctx.hdr_bam = sam_hdr_read(ctx.fp_bam); //read header

  if (ctx.hdr_bam == NULL) {
    return false;
  }
  ctx.aln = bam_init1();
  if (ctx.aln == NULL) {
    return false;
  }

  return true;
}

bool open_cpg_file(ContextConvert &ctx) {
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

void saving_hap(ContextConvert &ctx, vector<HT_s> &HT_vec) {

  string out_stream_name;
  if (ctx.fn_out) {
    out_stream_name = string(ctx.fn_out).substr(0, strlen(ctx.fn_out) - 3);
  } else {
    out_stream_name = "out.mhap";
  }
  ofstream out_stream(out_stream_name);

  unordered_map<string, bool> is_overlap;
  unordered_map<string, bool>::iterator itor;
  vector<HT_s>::iterator ht_itor;
  for (ht_itor = HT_vec.begin(); ht_itor != HT_vec.end(); ht_itor++) {
    (*ht_itor).ht_count = ctx.res_map[(*ht_itor).to_str()];
  }
  //sort
  sort(HT_vec.begin(), HT_vec.end(), comp_HT_vec);

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

  // compress mhap to gz
  string path_src = "./" + out_stream_name;
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

static void help() {
  cout << "Usage: mhaptools convert -i <in.bam>|<in.sam> -c <CpG.gz> [-r chr:beg-end | -b bed_file.bed ] [-n] [-o name.mhap.gz]" << endl;
  cout << "Options:" << endl;
  cout << "  -i  str  input file, SAM/BAM format, sorted by samtools" << endl;
  cout << "  -c  str  CpG file, gz format" << endl;
  cout << "  -r  str  region, e.g. chr1:2000-200000" << endl;
  cout << "  -b  str  bed file, one query region per line" << endl;
  cout << "  -n  flag non-directional, not group results by the direction of reads." << endl;
  cout << "  -o  str  output filename [out.mhap.gz]" << endl;
  cout << "Long options:" << endl;
  cout << "  -i  --input" << endl;
  cout << "  -c  --cpg" << endl;
  cout << "  -r  --region" << endl;
  cout << "  -b  --bed" << endl;
  cout << "  -n  --non-directional" << endl;
  cout << "  -o  --output" << endl;
  cout << "Examples:" << endl;
  cout << "- Convert the entire SAM/BAM file to mhap format:" << endl;
  cout << "  mhaptools convert -i in.bam -c CpG.gz" << endl << endl;
  cout << "- Convert the SAM/BAM file to mhap format within a region" << endl;
  cout << "  mhaptools convert -i in.bam -c CpG.gz -r chr1:2000-200000" << endl << endl;
  cout << "- Convert the SAM/BAM file to mhap format within several regions" << endl;
  cout << "  mhaptools convert -i in.bam -c CpG.gz -b bed_file.bed" << endl << endl;
}

int convert_fn_suffix_check(ContextConvert &ctx_cvt) {
  string gz_suffix = ".gz";
  string output_suffix = ".mhap.gz";
  string bed_suffix = ".bed";
  if (ctx_cvt.fn_cpg) {
    if (!is_suffix(ctx_cvt.fn_cpg, gz_suffix)) {
      hts_log_error("-c opt should be followed by a .gz file.");
      return 1;
    }
  }
  if (ctx_cvt.fn_out) {
    if (!is_suffix(ctx_cvt.fn_out, output_suffix)) {
      hts_log_error("-o opt should be followed by a .mhap.gz file.");
      return 1;
    }
  }
  if (ctx_cvt.fn_bed) {
    if (!is_suffix(ctx_cvt.fn_bed, bed_suffix)) {
      hts_log_error("-o opt should be followed by a .bed file.");
      return 1;
    }
  }
  return 0;
}

int main_convert(int argc, char *argv[]) {
  if (argc == optind) {
    help();
    return 0;
  }

  hts_log_info("create ContextConvert()");
  ContextConvert ctx = ContextConvert();

  hts_log_info("parse options");

  int long_index;

  static const char *opt_string = "i:b:c:r:o:nh";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "bed", optional_argument, NULL, 'b' },
      { "cpg", required_argument, NULL, 'c' },
      { "region", optional_argument, NULL, 'r' },
      { "output", optional_argument, NULL, 'o' },
      { "non-directional", optional_argument, NULL, 'n' },
      { "help", optional_argument, NULL, 'h' },
      { NULL, no_argument, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        ctx.fn_bam = optarg;
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
        ctx.fn_out = optarg;
        break;
      }
      case 'n': {
        ctx.non_directional = true;
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

  int ret;

  if (!convert_opt_check(ctx)) {
    hts_log_error("opt error");
    return 1;
  }

  if (convert_fn_suffix_check(ctx) == 1) {
    hts_log_error("filename suffix error.");
    return 1;
  }

    hts_log_info("open_bam_file(ctx)");
    ret = open_bam_file(ctx);
    if (!ret) {
      hts_log_error("open_bam_file():fail to open bam file");
      return EXIT_FAILURE;
    }
    hts_log_info("region_to_parse(ctx)");
    ret = region_to_parse(ctx);
    if (!ret) {
      hts_log_error("region_to_parse()");
    }
    if (ctx.region_to_parse == SINGLE_REGION) {
      ret = ctx.parse_region();
      hts_log_info("parse_region(): %s:%lld-%lld", ctx.hdr_bam->target_name[ctx.i_tid], ctx.i_beg, ctx.i_end);
    }
    if (!ret) {
      hts_log_error("Failed to parse region, check the -r specified.");
      return EXIT_FAILURE;
    }

    hts_log_info("open_cpg_file(ctx)");
    ret = open_cpg_file(ctx);
    if (!ret) {
      hts_log_error("open_cpg_file():fail to open cpg file");
      return EXIT_FAILURE;
    }

    vector<HT_s> HT_vec;

    hts_log_info("itor_sam(ctx).");
    cout << "Start processing..." << endl;
    HT_vec = itor_sam(ctx);
    if (HT_vec.size() == 0) {
      return 0;
    }

    hts_log_info("saving mhap");
    cout << "Saving..." << endl;
    saving_hap(ctx, HT_vec);

  return 0;
}
} // namespace std
