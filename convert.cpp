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
  int i_start = ctx.i_start - shift;
  int i_end = ctx.i_end + shift;
  int i = 0;
  int pos = 0;
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
            pos = atoi(string(cpg_line_sv.substr(0, i)).c_str());
            //该cpg起始位点在用户要求的位点范围内
            if(pos >= i_start && pos <= i_end) {
              ctx.cpg_pos.push_back(pos);
            }
          }
        }
      }
    }

  }
}

bool get_cpg_pos(context &ctx) {
  //读取cpg文件，如果cpg的开始位点在用户指定的范围内，将其放到cpg_pos中
  //???这边可以只限定开始位点吗？还是需要开始和结束都在范围内？
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

  start = ctx.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
  len = ctx.aln->core.l_qseq; //length of the read.
  end = start + len - 1;
  chr = ctx.hdr_bam->target_name[ctx.aln->core.tid] ; //contig name (chromosome)

  quality = bam_get_seq(ctx.aln); //quality string
  map_quality = ctx.aln->core.qual ; //mapping quality


  for(int i=0; i< len ; i++){
    seq.push_back(seq_nt16_str[bam_seqi(quality,i)]); //gets nucleotide id and converts them into IUPAC id.
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
    //??？这个位置要确认一下确实可以运行

  } else if(strcmp(ctx.aligner, "MAQ") == 0) {

  } else {
    cout << "Only BSMAP, BISMARK and MAQ are supported." << endl;
    return false;
  }
  return true;
}

bool sam_read::haplo_type() {
  uint32_t r_pos;
  vector<int> cpg;
  string hap_seq = "";

  for(int i = 0; i < ctx->cpg_pos.size(); i++) {
    int pos = ctx->cpg_pos[i];
    if(pos < start || pos > end) {
      continue;
    }
    if(WC == DIRECTION_PLUS) {
      r_pos = pos - start;
    } else {
      r_pos = pos - start + 1;
    }
    if(r_pos > len) {
      continue;
    }
    cpg.push_back(pos);
  hap_seq += seq[r_pos];
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
      XM_tag = ctx->aln->data + i + 2;
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

void itor_sam(context &ctx) {

  while(sam_read1(ctx.fp_bam, ctx.hdr_bam, ctx.aln) > 0){
    sam_read sam_r = sam_read();
    sam_r.init(ctx);
    cout << sam_r.quality << endl;
    sam_r.haplo_type();
    sam_r.seq.clear();
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
    bool ret;
    ret = ctx.parse_region();
    if(!ret) {
      //fail to parse the query
      ctx.print_region();
      return EXIT_FAILURE;
    }

    ret = open_bam_file(ctx);
    if(!ret) {
      return EXIT_FAILURE;
    }

    ret = open_cpg_file(ctx);
    if(!ret) {
      return EXIT_FAILURE;
    }

    itor_sam(ctx);
  }

  return EXIT_SUCCESS;
  }

