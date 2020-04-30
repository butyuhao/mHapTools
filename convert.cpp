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
  this->ctx = &ctx;

  len = ctx.aln->core.l_qseq; //length of the read.

  chr = ctx.hdr_bam->target_name[ctx.aln->core.tid] ; //checked
  start = ctx.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
  end = start + len - 1;


  qual = bam_get_qual(ctx.aln); //quality string
  qname = bam_get_qname(ctx.aln);
  flag = ctx.aln->core.flag; //
  //char *rname;
  mapq = ctx.aln->core.qual ; //mapping quality
  cigar = bam_get_cigar(ctx.aln);
//  for(int i=0; i < ctx.aln->core.n_cigar;++i){
//    int icigar = cigar[i];
//    printf("%d%d\n",bam_cigar_op(icigar),bam_cigar_oplen(icigar));
//  }
  //*rnext;
  //pnext;
  //tlen;

  //get the seq pointer and put it into a vector
  seq_p = bam_get_seq(ctx.aln);
  for(int i=0; i < len;++i){
    seq.push_back(_base[bam_seqi(seq_p, i)]);
  }
  if(strcmp(ctx.aligner, "BISMARK") == 0) {
    _get_bismark_std();

    if(_get_XM_tag()) {
      _get_bismark_QC();
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
  for(int i = 0; i < seq.size(); i++) {
  }
  for(int i = 0; i < ctx->cpg_pos.size(); i++) {
    pos = ctx->cpg_pos[i];
    if(pos < start || pos > end) {
      continue;
    }
    if(WC == DIRECTION_PLUS) {
      r_pos = pos - start;
    } else {
      r_pos = pos - start + 1;
    }
    if(r_pos >= len) {
      continue;
    }
    cpg.push_back(pos);
    hap_seq += seq[r_pos];
    hap_qual.push_back(qual[r_pos]);
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
    if(WC == DIRECTION_PLUS) {
      if(nucleobases == 'C') {
        hap_met += '1';
      } else if (nucleobases == 'T') {
        hap_met += '0';
      } else {
        hap_met += nucleobases;
        QC = false;
      }
    } else if (WC == DIRECTION_MINUS) {
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
    HT = HT_s(chr, _cpg[0], _cpg[_cpg.size() - 1], _hap_met, 1, WC);
  }
}

void sam_read::_get_bismark_std() {
  if(ctx->aln->core.flag == 99 || ctx->aln->core.flag == 147) {
    WC = DIRECTION_PLUS;
  } else if (ctx->aln->core.flag == 83 || ctx->aln->core.flag == 163) {
    WC = DIRECTION_MINUS;
  }
}

bool sam_read::_get_XM_tag() {

  for(int i = 0; i < ctx->aln->l_data - 1; i++) {
    //get XM tag. if success, return true, else, return false
    if(ctx->aln->data[i] == 'X' && ctx->aln->data[i + 1] == 'M') {
      XM_tag = ctx->aln->data + i + 3;
      //cout << XM_tag;
      return true;
    }
  }
  return false;
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

void sam_read::_get_bismark_QC() {
  for(int i = 0; XM_tag[i] != '\0'; i++) {
    if(XM_tag[i] == 'X' || XM_tag[i] == 'H' || XM_tag[i] == 'U') {
      QC = false;
    }
  }
  QC = true;
}

uint8_t* sam_read::_my_bam_get_seq(context &ctx) {
  bool pre_split_flag = 0;
  bool split_flag = 0;
  int split_cnt = 0;
  cout<<hex;
  for(int i = 0; i < ctx.aln->l_data; i++) {
    if(ctx.aln->data[i] == '\0') {
      cout << "X";
    } else {
      cout << ctx.aln->data[i];
    }
  }
  return NULL;
}

void itor_sam(context &ctx) {

  while(sam_read1(ctx.fp_bam, ctx.fp_bam->bam_header, ctx.aln) > 0){
    if(ctx.aln->core.pos + 1 == 48954) {
      sam_read sam_r = sam_read();
      sam_r.init(ctx);
      sam_r.haplo_type();
      sam_r.seq.clear();
    }

  }
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
    itor_sam(ctx);
  }

  return EXIT_SUCCESS;
  }

