//
// Created by Yuhao Dan on 2020/4/13.
//
#include "convert.h"
#include <sstream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include <string_view>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>

namespace std {

bool context::parse_region() {
  string_view arg_sv = string_view (region);
  bool is_get_chr = false;
  bool is_get_start = false;
  bool is_get_end = false;
  int i;
  for (i = 0; i < arg_sv.size(); i++) {
    if (arg_sv[i] == ':') {
      //get the chr name
      i_chr = string(arg_sv.substr(0, i));
      is_get_chr = true;
      break;
    }
  }
  if (is_get_chr) {
    for (i++; i < arg_sv.size(); i++) {
      if (arg_sv[i] >= '0' && arg_sv[i] <= '9') {
        //Initialize i_start
        if (i_start == -1) {
          i_start = 0;
        }
        i_start = i_start * 10 + (arg_sv[i] - '0');
      } else if (arg_sv[i] == '-') {
        is_get_start = true;
        i_start += 1;
        break;
      } else {
        return false;
      }
    }
  }
  if (is_get_start) {
    for (i++; i < arg_sv.size(); i++) {
      if (arg_sv[i] >= '0' && arg_sv[i] <= '9') {
        //Initialize i_end
        if (i_end == -1) {
          i_end = 0;
        }
        i_end = i_end * 10 + (arg_sv[i] - '0');
      } else {
        return false;
      }
    }
  }
  if (i_end != -1) {
    is_get_end = true;
  }
  if (is_get_end && is_get_start && is_get_chr) {
    return true;
  }
  return false;
}

void context::print_region() {
  cout << "Chr:" << i_chr << " Start:" << i_start << " End:" << i_end << endl;
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

void parse_cpg_line(context &ctx, uint32_t shift = 500) {
  //匹配返回start位置，不匹配返回-1
  int tab_cnt = 0;
  uint32_t i_start;
  string_view cpg_line_sv = string_view(ctx.fp_cpg->line.s);
  //todo 现在限定了start的下限，避免减去500后溢出，之后应该确保所有输入都能被正确处理
  if (ctx.i_start >= 500) {
    i_start = ctx.i_start - shift;
  } else {
    i_start = 0;
  }
  uint32_t i_end = ctx.i_end + shift;
  int i = 0;
  uint32_t cpg_start = 0;
  uint32_t cpg_end = 0;
  for (; i < cpg_line_sv.size(); i++) {
    if (cpg_line_sv[i] == '\t') {
      //指定的ichr与当前读取到的cpg的ichr相同
      if (cpg_line_sv.substr(0, i).compare(ctx.i_chr) == 0) {
        //将ichr部分从line中去除
        cpg_line_sv = cpg_line_sv.substr(i + 1);
        //curser清零
        i = 0;
        //找到原line中的第2个tab位置
        for (; i < cpg_line_sv.size(); i++) {
          if (cpg_line_sv[i] == '\t') {
            cpg_start = atoi(string(cpg_line_sv.substr(0, i)).c_str());
            cpg_end = atoi(string(cpg_line_sv.substr(i, cpg_line_sv.size() - i)).c_str());
            //确保cpg位点在用户指定的范围内。
            if (cpg_start >= i_start && cpg_end <= i_end) {
              ctx.cpg_pos.push_back(cpg_start);
            }
          }
        }
      }
    }
  }
}

bool load_cpg_with_idx(context &ctx, uint32_t shift = 500) {
  //concat name of the tbi file
  uint32_t  i_start;

  if (ctx.i_start >= 500) {
    i_start = ctx.i_start - shift;
  } else {
    i_start = 0;
  }
  uint32_t i_end = ctx.i_end + shift;

  string cpg_idx_fn = string(ctx.fn_cpg) + string(".tbi");
  if (cpg_idx_fn.size() == 0) {
    return false;
  }
  //get the tbi index
  tbx_t *idx_cpg = tbx_index_load(cpg_idx_fn.c_str());
  int tid = tbx_name2id(idx_cpg, ctx.i_chr.c_str());
  hts_itr_t *cpg_itr = tbx_itr_queryi(idx_cpg, tid, i_start, i_end);
  kstring_t ksbuf = {0, 0, NULL};
  string_view cpg_line_sv;
  u_int32_t cpg_start = 0;
  u_int32_t cpg_end = 0;


  while(tbx_itr_next(ctx.fp_cpg, idx_cpg, cpg_itr, &ksbuf) >= 0) {
    cpg_line_sv = string_view(ksbuf.s);
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



bool get_cpg_pos(context &ctx) {
  //读取cpg文件，如果cpg的开始位点在用户指定的范围内，将其放到cpg_pos中
  int  ret;
  //cpg文件为XXX.gz,则对应文件为XXX.gz.tbi
  if (load_cpg_with_idx(ctx)) {
    //使用带有index的方法读取
    return true;
  } else {
    //使用不带有index的方法读取
    ret = hts_getline(ctx.fp_cpg, KS_SEP_LINE, &ctx.fp_cpg->line);
    parse_cpg_line(ctx);
    while(ret >= 0) {
      ret = hts_getline(ctx.fp_cpg, KS_SEP_LINE, &ctx.fp_cpg->line);
      parse_cpg_line(ctx);
    }
  }

  return true;
}

bool SamRead::init(context &ctx) {
  bool ret;

  this->ctx = &ctx;

  read_len = ctx.aln->core.l_qseq; //length of the read.

  read_chr = ctx.hdr_bam->target_name[ctx.aln->core.tid] ; //checked
  read_start = ctx.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
  read_end = read_start + read_len - 1;
  //判断read在用户指定的范围内
  if (strcmp(ctx.i_chr.c_str(), read_chr) != 0) {
    return false;
  }
  if (!((read_start >= ctx.i_start && read_start <= ctx.i_end) || (read_end > ctx.i_start && read_end <= ctx.i_end))) {
    return false;
  }
  read_qual = bam_get_qual(ctx.aln); //quality string
  read_name = bam_get_qname(ctx.aln);
  flag = ctx.aln->core.flag; //

  read_map_quality = int(ctx.aln->core.qual); //mapping quality
  read_cigar = bam_get_cigar(ctx.aln);
  //get the seq pointer and put it into a vector
  uint8_t *seq_p = bam_get_seq(ctx.aln);
  for (int i=0; i < read_len; ++i){
    seq.push_back(kbase[bam_seqi(seq_p, i)]);
  }
  if (strcmp(ctx.aligner, "BISMARK") == 0) {
    ret = _get_bismark_std();

    if (!ret) {
      cout << "Error:_get_bismark_std()." << endl;
      return false;
    }

    if (_get_XM_tag(ctx)) {
      ret = _get_bismark_QC(ctx);
      if (!ret) {
        //cout << "ERROR:_get_bismark_QC()." << endl;
        return false;
      }
    } else {
      cout << "Error:XM tag is required in SAM." << endl;
      return false;
    }
  } else if(strcmp(ctx.aligner, "BSMAP") == 0) {
    if (!_get_ZS_tag(ctx)) {
      return false;
    }

  } else if(strcmp(ctx.aligner, "UNKNOWN") == 0) {
    //todo 使用unknown的时候，与BISMARK比较可能多出一些结果，需要排查一下原因。
    read_WC = DIRECTION_UNKNOWN;
  } else {
    cout << "Only BSMAP, BISMARK and UNKNOWN are supported." << endl;
    return false;
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
      break;
    }
    if (read_WC == DIRECTION_PLUS || read_WC == DIRECTION_UNKNOWN) {
      r_pos = pos - read_start;
    } else {
      r_pos = pos - read_start + 1;
    }
    if (r_pos >= read_len) {
      continue;
    }
    cpg.push_back(pos);
    if (r_pos > read_len) {
      cout << "r_pos > len" << endl;
      continue;
    }
    hap_seq += seq[r_pos];
    hap_qual.push_back(read_qual[r_pos]);
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

bool SamRead::_get_XM_tag(context &ctx) {

  ctx.bam_aux_p = bam_aux_get(ctx.aln, "XM");
  if (!ctx.bam_aux_p) {
    return false;
  }
  return true;
}

bool SamRead::_get_ZS_tag(context &ctx) {

  ctx.bam_aux_p = bam_aux_get(ctx.aln, "ZS");
  if (!ctx.bam_aux_p) {
    return false;
  }
  ZS_tag = bam_aux2Z(ctx.bam_aux_p);
  if (!ZS_tag) {
    return false;
  }
  if (*ZS_tag == '+') {
    read_WC = DIRECTION_PLUS;
  } else if (*ZS_tag == '-'){
    read_WC = DIRECTION_MINUS;
  } else {
    return false;
  }
  return true;
}

bool SamRead::_get_bismark_QC(context &ctx) {
  XM_tag = bam_aux2Z(ctx.bam_aux_p);
  if (!XM_tag) {
    return false;
  }
  string_view XM_tag_sv  = string_view(XM_tag);
  for (auto c : XM_tag_sv) {
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

struct cmp
{
  bool operator()(const pair<string,u_int32_t> &p1,const pair<string,u_int32_t> &p2)
  {
    return p1.second < p2.second;
  }
};

map<string, int> itor_sam(context &ctx) {
  //cout << "enter itor_sam()" << endl;

  map<string, vector<SamRead>> sam_map;
  map<string, int> res_map;
  map<string, vector<SamRead>>::iterator iter;
  vector<HT_s> res_l;
  static string bam_idx_fn = string(ctx.fn_bam) + ".bai";
  ctx.idx_bam = sam_index_load(ctx.fp_bam, bam_idx_fn.c_str());
  int chr = sam_hdr_name2tid(ctx.hdr_bam, ctx.i_chr.c_str());
  hts_itr_t *sam_itr = sam_itr_queryi(ctx.idx_bam, chr, ctx.i_start, ctx.i_end);
  while(sam_itr_next(ctx.fp_bam, sam_itr, ctx.aln) >= 0) {
    if (ctx.aln->core.flag & BAM_FQCFAIL || ctx.aln->core.flag & BAM_FUNMAP || ctx.aln->core.flag & BAM_FDUP
        || ctx.aln->core.flag & BAM_FSECONDARY || ctx.aln->core.flag & BAM_FSUPPLEMENTARY) {
      //cout << "flag continue" << endl;
      continue;
    }
    SamRead sam_r = SamRead();
    int ret = sam_r.init(ctx);
    //cout << "start:" << sam_r.read_start << " end:" << sam_r.read_end << endl;
    if (!ret) {
      continue;
    }

    string qname = string(sam_r.read_name);

    sam_r.haplo_type();
    if (!sam_r.QC) {
      //cout << "!QC" << endl;
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
  for (auto sam_l :  sam_map) {
    if (sam_l.second.size() == 2) {
      SamRead samF = sam_l.second[0];
      SamRead samR = sam_l.second[1];
      if (paired_end_check(samF, samR)) {
        HT_s ht = paired_end_merge(samF, samR);
        res_l.push_back(ht);
        if (test_mode) {
          //cout << "Merged:" << samF.read_name << " and " << samR.read_name << endl;
        }
      } else {
        for (int i = 0; i < sam_l.second.size(); i++) {
          res_l.push_back(sam_l.second[i].HT);
          //cout << "overlap:";
          //cout << sam_l.second[i].read_name << endl;
        }
      }
    } else {
      for (int i = 0; i < sam_l.second.size(); i++) {
        res_l.push_back(sam_l.second[i].HT);
        //cout << "single" << sam_l.second[i].read_name << endl;
      }
    }
  }

  for (auto _ht: res_l) {
    string ht_id = _ht.to_str();

    map<string, int>::iterator res_map_itor;
    ctx.res_map_sort[ht_id] = _ht.get_h_start();
    res_map_itor = res_map.find(ht_id);
    if (res_map_itor == res_map.end()){
      res_map[ht_id] = 1;
    } else {
      res_map[ht_id] += 1;
    }
  }
  for (auto rms : ctx.res_map_sort) {
    ctx.vt.push_back(pair<string, u_int32_t >(rms.first, rms.second));
  }
  sort(ctx.vt.begin(),ctx.vt.end(),cmp());

  return res_map;
}


inline context::~context() {
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
  if (idx_cpg) {
    tbx_destroy(idx_cpg);
  }
}

bool open_bam_file(context &ctx) {

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

bool open_cpg_file(context &ctx) {
  ctx.fp_cpg = hts_open(ctx.fn_cpg, "r");
  if (ctx.fp_cpg == NULL) {
    return false;
  }
  return true;
}


int main_convert(int argc, char *argv[]) {
//TODO(butyuhao@foxmail.com): 增加检查option合法性的部分
//TODO(butyuhao@foxmail.com): 增加日志

  context ctx = context();

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
        cout << optarg << endl;
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

  if (ctx.region) {
    if (test_mode) {
      cout << ctx.region << endl;
    }
    bool ret;
    ret = ctx.parse_region();
    if (test_mode) {
      ctx.print_region();
    }

    if (!ret) {
      //fail to parse the query
      ctx.print_region();
      return EXIT_FAILURE;
    }
    ret = open_bam_file(ctx);

    if (!ret) {
      cout << "Fail to open bam file" << endl;
      return EXIT_FAILURE;
    }

    ret = open_cpg_file(ctx);
    if (!ret) {
      cout << "Fail to open cpg file" << endl;
      return EXIT_FAILURE;
    }

    get_cpg_pos(ctx);
    map<string, int> res_map;
    res_map = itor_sam(ctx);
    string out_stream_name;
    if (ctx.output_path) {
      out_stream_name = ctx.output_path;
    } else {
      out_stream_name = "out.hap";
    }
    ofstream out_stream(out_stream_name);

    for (auto p: ctx.vt) {
      char direction = '+';
      if (p.first[p.first.size() - 1] == '1') {
        direction = '-';
      } else if (p.first[p.first.size() - 1] == DIRECTION_UNKNOWN + '0') {
        direction = '*';
      }
      string out_string = p.first.substr(0, p.first.size() - 1);
      out_string = out_string + '\t' + to_string(res_map[p.first]) + '\t' + direction;
      out_stream << out_string << endl;

    }

    out_stream.close();
  }

  return EXIT_SUCCESS;
}
} // namespace std
