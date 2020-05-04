//
// Created by Yuhao Dan on 2020/4/13.
//

#include "convert.h"



using namespace std;

bool context::parse_region() {
  string_view arg_sv = string_view (region);
  bool is_get_chr = false;
  bool is_get_start = false;
  bool is_get_end = false;
  int i;
  for(i = 0; i < arg_sv.size(); i++) {
    if(arg_sv[i] == ':') {
      //get the chr name
      i_chr = string(arg_sv.substr(0, i));
      is_get_chr = true;
      break;
    }
  }
  if(is_get_chr) {
    for(i++; i < arg_sv.size(); i++) {
      if(arg_sv[i] >= '0' && arg_sv[i] <= '9') {
        //Initialize i_start
        if(i_start == -1) {
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
  if(is_get_start) {
    for(i++; i < arg_sv.size(); i++) {
      if(arg_sv[i] >= '0' && arg_sv[i] <= '9') {
        //Initialize i_end
        if(i_end == -1) {
          i_end = 0;
        }
        i_end = i_end * 10 + (arg_sv[i] - '0');
      } else {
        return false;
      }
    }
  }
  if(i_end != -1) {
    is_get_end = true;
  }
  if(is_get_end && is_get_start && is_get_chr) {
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
    { NULL, no_argument, NULL, 0 }
};

void parse_cpg_line(context &ctx, int shift = 500) {
  //匹配返回start位置，不匹配返回-1
  int tab_cnt = 0;
  string_view cpg_line_sv = string_view(ctx.fp_cpg->line.s);
  uint32_t i_start = ctx.i_start - shift;
  uint32_t i_end = ctx.i_end + shift;
  int i = 0;
  int i_ = 0;
  uint32_t cpg_start = 0;
  uint32_t cpg_end = 0;
  for(; i < cpg_line_sv.size(); i++) {
    if(cpg_line_sv[i] == '\t') {
      //指定的ichr与当前读取到的cpg的ichr相同
      if(cpg_line_sv.substr(0, i).compare(ctx.i_chr) == 0) {
        //将ichr部分从line中去除
        cpg_line_sv = cpg_line_sv.substr(i + 1);
        //curser清零
        i = 0;
        //找到原line中的第2个tab位置
        for(; i < cpg_line_sv.size(); i++) {
          if(cpg_line_sv[i] == '\t') {
            cpg_start = atoi(string(cpg_line_sv.substr(0, i)).c_str());
            cpg_end = atoi(string(cpg_line_sv.substr(i, cpg_line_sv.size() - i)).c_str());
            //确保cpg位点在用户指定的范围内。
            if(cpg_start >= i_start && cpg_end <= i_end) {
              ctx.cpg_pos.push_back(cpg_start);
            }
          }
        }
      }
    }

  }
}

bool get_cpg_pos(context &ctx) {
  //读取cpg文件，如果cpg的开始位点在用户指定的范围内，将其放到cpg_pos中
  int  ret;
  ret = hts_getline(ctx.fp_cpg, KS_SEP_LINE, &ctx.fp_cpg->line);
  parse_cpg_line(ctx);
  while(ret >= 0) {
    ret = hts_getline(ctx.fp_cpg, KS_SEP_LINE, &ctx.fp_cpg->line);
    parse_cpg_line(ctx);
  }
  return ret;
}

bool sam_read::init(context &ctx) {
  bool ret;

  this->ctx = &ctx;

  read_len = ctx.aln->core.l_qseq; //length of the read.

  read_chr = ctx.hdr_bam->target_name[ctx.aln->core.tid] ; //checked
  read_start = ctx.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
  read_end = read_start + read_len - 1;
  //判断read在用户指定的范围内
  if(strcmp(ctx.i_chr.c_str(), read_chr) != 0) {
    return false;
  }
  if(read_start < ctx.i_start || read_end > ctx.i_end) {
    return false;
  }
  read_qual = bam_get_qual(ctx.aln); //quality string
  read_name = bam_get_qname(ctx.aln);
  flag = ctx.aln->core.flag; //

  read_map_quality = int(ctx.aln->core.qual); //mapping quality
  read_cigar = bam_get_cigar(ctx.aln);
//  for(int i=0; i < ctx.aln->core.n_cigar;++i){
//    int icigar = cigar[i];
//    printf("%d%d\n",bam_cigar_op(icigar),bam_cigar_oplen(icigar));
//  }

  //get the seq pointer and put it into a vector
  uint8_t *seq_p = bam_get_seq(ctx.aln);
  for(int i=0; i < read_len; ++i){
    seq.push_back(_base[bam_seqi(seq_p, i)]);
  }
  if(strcmp(ctx.aligner, "BISMARK") == 0) {
    ret = _get_bismark_std();

    if(!ret) {
      return false;
    }

    if(_get_XM_tag(ctx)) {
      ret = _get_bismark_QC(ctx);
      if(!ret) {
        return false;
      }
    } else {
      cout << "Error:XM tag is required in SAM." << endl;
      return false;
    }
  } else if(strcmp(ctx.aligner, "BSMAP") == 0) {
    //todo

  } else if(strcmp(ctx.aligner, "MAQ") == 0) {
    //todo
  } else {
    cout << "Only BSMAP, BISMARK and MAQ are supported." << endl;
    return false;
  }
  return true;
}

bool sam_read::haplo_type() {
  uint32_t r_pos;
  vector<uint32_t> cpg;
  vector<int8_t> hap_qual;
  string hap_seq = "";
  uint32_t pos;
  for(int i = 0; i < ctx->cpg_pos.size(); i++) {
    pos = ctx->cpg_pos[i];
    if(pos < read_start) {
      continue;
    }
    if(pos > read_end) {
      break;
    }
    if(read_WC == DIRECTION_PLUS) {
      r_pos = pos - read_start;
    } else {
      r_pos = pos - read_start + 1;
    }
    if(r_pos >= read_len) {
      continue;
    }
    cpg.push_back(pos);
    if(r_pos > read_len) {
      cout << "r_pos > len" << endl;
      continue;
    }
    hap_seq += seq[r_pos];
    hap_qual.push_back(read_qual[r_pos]);
  }
  if(cpg.size() == 0) {
    QC = false;
  }
  _hap_seq = hap_seq;
  _hap_qual = hap_qual;
  _cpg = cpg;
  string hap_met = "";
  for(int i = 0; i < hap_seq.size(); i++) {
    char nucleobases = hap_seq[i];
    if(read_WC == DIRECTION_PLUS) {
      if(nucleobases == 'C') {
        hap_met += '1';
      } else if (nucleobases == 'T') {
        hap_met += '0';
      } else {
        hap_met += nucleobases;
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

bool sam_read::_get_bismark_std() {
  if(ctx->aln->core.flag == 99 || ctx->aln->core.flag == 147) {
    read_WC = DIRECTION_PLUS;
  } else if (ctx->aln->core.flag == 83 || ctx->aln->core.flag == 163) {
    read_WC = DIRECTION_MINUS;
  } else {
    cout << "sam_read::_get_bismark_std() Unknown flag" << endl;
    return false;
  }
  return true;
}

bool sam_read::_get_XM_tag(context &ctx) {

  ctx.bam_aux_p = bam_aux_get(ctx.aln, "XM");
  if(!ctx.bam_aux_p) {
    return false;
  }
  return true;
}

bool sam_read::_get_ZS_tag() {
  for(int i = 0; i < ctx->aln->l_data - 1; i++) {
    //get XM tag. if success, return EXIT_SUCCESS, else, return EXIT_FAILURE
    if(ctx->aln->data[i] == 'Z' && ctx->aln->data[i + 1] == 'S') {
      ZS_tag = ctx->aln->data + i + 2;
      return true;
    }
  }
  return false;
}

bool sam_read::_get_bismark_QC(context &ctx) {
  XM_tag = bam_aux2Z(ctx.bam_aux_p);
  if(!XM_tag) {
    return false;
  }
  string_view XM_tag_sv  = string_view(XM_tag);
  for(auto c : XM_tag_sv) {
    if(c == 'X' || c == 'H' || c == 'U') {
      return false;
    }
  }
  return true;
}

bool paired_end_check(sam_read &samF, sam_read &samR) {
  if(strcmp(samF.read_name, samR.read_name) != 0) {
    return false;
  }
  if (strcmp(samF.read_chr, samR.read_chr) != 0) {
    return false;
  }
  bool checkF = samF._cpg[samF._cpg.size() - 1] >= samR._cpg[0];
  bool checkR = samR._cpg[samR._cpg.size() - 1] >= samF._cpg[0];
  return (checkF && checkR) || (!checkF && !checkR);
}

HT_s paired_end_merge(sam_read &samF, sam_read &samR) {
  map<uint32_t, char> merged_seq;
  map<uint32_t, int8_t> merged_qual;
  map<uint32_t, char> merged_met;
  map<uint32_t, char>::iterator merged_met_itor;
  for(int i = 0; i < samF._cpg.size(); i++) {
    uint32_t pos = samF._cpg[i];
    merged_seq[pos] = samF._hap_seq[i];
    merged_qual[pos] = samF._hap_qual[i];
    merged_met[pos] = samF._hap_met[i];
  }
  for(int i = 0; i < samR._cpg.size(); i++) {
    uint32_t pos = samR._cpg[i];
    merged_met_itor = merged_met.find(pos);
    if(merged_met_itor == merged_met.end() || samR._hap_qual[i] > merged_qual[pos]) {
      merged_seq[pos] = samR._hap_seq[i];
      merged_qual[pos] = samR._hap_qual[i];
      merged_met[pos] = samR._hap_met[i];
    }
  }
  cout << "cpg:";
  merged_met_itor = merged_met.begin();
  while(merged_met_itor != merged_met.end()) {
    cout << " " << merged_met_itor->first;
  }
  cout << endl;

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



map<string, int> itor_sam(context &ctx) {
  cout << "enter itor_sam()" << endl;
  map<string, vector<sam_read>> sam_map;
  map<string, int> res_map;
  map<string, vector<sam_read>>::iterator iter;
  vector<HT_s> res_l;
  while(sam_read1(ctx.fp_bam, ctx.fp_bam->bam_header, ctx.aln) > 0){
    //todo 在生成sam对象之前先来一个is_valid()来判断，如果不符合要求直接过。
    if(ctx.aln->core.flag & BAM_FQCFAIL || ctx.aln->core.flag & BAM_FUNMAP || ctx.aln->core.flag & BAM_FDUP
      || ctx.aln->core.flag & BAM_FSECONDARY || ctx.aln->core.flag & BAM_FSUPPLEMENTARY) {
      cout << "flag continue" << endl;
      continue;
    }
    sam_read sam_r = sam_read();
    int ret = sam_r.init(ctx);
    if(!ret) {
      continue;
    }

    string qname = string(sam_r.read_name);

    //结果限定在用户指定的范围内

    sam_r.haplo_type();
    if(!sam_r.QC) {
      continue;
    }

    iter = sam_map.find(qname);
    if(iter == sam_map.end()) {
      cout << "no" << endl;
      vector<sam_read> v;
      v.push_back(sam_r);
      sam_map[qname] = v;
    } else {
      cout << "yes" << endl;
      vector<sam_read> v;
      v = sam_map[qname];
      v.push_back(sam_r);
      sam_map[qname] = v;
    }
  }
  iter = sam_map.begin();
  while(iter != sam_map.end()) {
    if(iter->second.size() == 2) {
      sam_read samF = iter->second[0];
      sam_read samR = iter->second[1];
      if(paired_end_check(samF, samR)) {
        HT_s ht = paired_end_merge(samF, samR);
        res_l.push_back(ht);
        if(test_mode) {
          cout << "Merged:" << samF.read_name << " and " << samR.read_name << endl;
        }
      } else {
        for(int i = 0; i < iter->second.size(); i++) {
          res_l.push_back(iter->second[i].HT);
          cout << "overlap:";
          cout << iter->second[i].read_name << endl;
        }
      }
    } else {
      for(int i = 0; i < iter->second.size(); i++) {
        res_l.push_back(iter->second[i].HT);
        cout << "single" << iter->second[i].read_name << endl;
      }
    }
    iter++;
  }
  for(auto _ht: res_l) {
    string ht_id = _ht.to_str();
    map<string, int>::iterator res_map_itor;
    res_map_itor = res_map.find(ht_id);
    if(res_map_itor == res_map.end()){
      res_map[ht_id] = 1;
    } else {
      res_map[ht_id] += 1;
    }
  }
  cout << "leave itor_sam()" << endl;
  return res_map;
}


inline context::~context() {
  //to_do明确一下哪些指针需要被关掉。
  if(fp_bam) {
    hts_close(fp_bam);
  }
  if(fp_cpg) {
    hts_close(fp_cpg);
  }
  if(aln) {
    bam_destroy1(aln);
  }
}

inline void context::init_ctx() {
  //Initialize ctx
  fp_bam = NULL;
  fp_cpg = NULL;
  hdr_bam = NULL;/* -i bam文件头部的指针
  aln = NULL;
  fn_bam = NULL;
  fn_cpg = NULL;
  bam_path = NULL;   /* -i option */
  cpg_path = NULL;      /* -c option */
  output_path = NULL;  /* -o option */
  aligner = NULL;      /* -a option */
  bed_file = NULL;     /* -b option
  region = NULL;       /* -r option */
}

bool open_bam_file(context &ctx) {

  ctx.fp_bam = hts_open(ctx.bam_path, "r");
  ctx.hdr_bam = sam_hdr_read(ctx.fp_bam); //read header
  if(ctx.hdr_bam == NULL) {
    return false;
  }
  ctx.aln = bam_init1();
  if(ctx.aln == NULL) {
    return false;
  }
  return true;
}

bool open_cpg_file(context &ctx) {
  ctx.fp_cpg = hts_open(ctx.cpg_path, "r");
  if(ctx.fp_cpg == NULL) {
    return false;
  }
  return true;
}

int main_convert(int argc, char *argv[]) {


  context ctx = context();
  //初始化context中的变量
  ctx.init_ctx();

  int long_index;

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while(opt != -1) {
    switch(opt) {
      case 'i':
        ctx.bam_path = optarg;
        break;

      case 'a':
        ctx.aligner = optarg;
        break;

      case 'b':
        cout << optarg << endl;
        break;

      case 'c':
        ctx.cpg_path = optarg;
        break;

      case 'r':
        //cout << optarg << endl;
        ctx.region = optarg;
        break;

      default:
        break;
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }

  if (ctx.region) {
    if(test_mode) {
      cout << "start parsing region" << endl;
      cout << ctx.region << endl;
    }
    bool ret;
    ret = ctx.parse_region();

    if(test_mode) {
      ctx.print_region();
    }

    if(!ret) {
      //fail to parse the query
      ctx.print_region();
      return EXIT_FAILURE;
    }
    ret = open_bam_file(ctx);
    if(!ret) {
      cout << "Fail to open bam file" << endl;
      return EXIT_FAILURE;
    }

    ret = open_cpg_file(ctx);
    if(!ret) {
      cout << "Fail to open cpg file" << endl;
      return EXIT_FAILURE;
    }
    get_cpg_pos(ctx);
    map<string, int> res_map;
    res_map = itor_sam(ctx);

    ofstream out_stream("out.hap");

    for(auto r_map: res_map) {
      char direction = '+';
      if (r_map.first[r_map.first.size() - 1] == sam_read::DIRECTION_MINUS) {
        direction = '-';
      }
      string out_string = r_map.first.substr(0, r_map.first.size() - 1);
      out_string = out_string + '\t' + to_string(r_map.second) + '\t' + direction;
      cout << out_string << endl;
      out_stream << out_string << endl;

    }

    out_stream.close();
  }

  return EXIT_SUCCESS;
  }

