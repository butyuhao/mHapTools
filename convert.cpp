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
        return EXIT_FAILURE;
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
        return EXIT_FAILURE;
      }
    }
  }
  if(i_end != -1) {
    is_get_end = true;
  }
  if(is_get_end && is_get_start && is_get_chr) {
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
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
  int  ret;
  ret = hts_getline(ctx.fp_cpg, KS_SEP_LINE, &ctx.fp_cpg->line);
  parse_cpg_line(ctx);
  cout << ctx.fp_cpg->line.s << endl;
  cout << ctx.fp_cpg->line.l << endl;
  cout << ctx.fp_cpg->line.m << endl;
  while(ret >= 0) {
    ret = hts_getline(ctx.fp_cpg, KS_SEP_LINE, &ctx.fp_cpg->line);
    parse_cpg_line(ctx);
  }
}

void parse_sam(context &ctx) {
  int ret;
  ret = get_cpg_pos(ctx);
}

inline context::~context() {
  //to_do明确一下哪些指针需要被关掉。
  if(fp_bam) {
    hts_close(fp_bam);
  }
  if(fp_cpg) {
    hts_close(fp_cpg);
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
    hts_close(ctx.fp_bam);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

bool open_cpg_file(context &ctx) {
  ctx.fp_cpg = hts_open(ctx.cpg_path, "r");
  if(ctx.fp_cpg == NULL) {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int main_convert(int argc, char *argv[]) {
  context ctx = context();
  ctx.init_ctx();

  int long_index;

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while(opt != -1) {
    switch(opt) {
      case 'i':
        ctx.bam_path = optarg;
        break;

      case 'a':
        cout << optarg << endl;
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
    if(ret == EXIT_FAILURE) {
      //fail to parse the query
      ctx.print_region();

      return EXIT_FAILURE;
    }
    ret = open_bam_file(ctx);
    if(ret == EXIT_FAILURE) {

      return EXIT_FAILURE;
    }
    ret = open_cpg_file(ctx);
    if(ret == EXIT_FAILURE) {

      return EXIT_FAILURE;
    }
    parse_sam(ctx);
  }

  return EXIT_SUCCESS;
  }
  //
//
//  //header parse
//  char *tar = ctx.bamHdr->text ;
//  uint32_t *tarlen = ctx.bamHdr->target_len ;
//
//  cout << tar << endl;
//  // content parse
//  char *chrom = "1";
//  int locus = 2000;
//  int comp;
//
//  while(sam_read1(ctx.fp_bam,ctx.bamHdr,ctx.aln) > 0){
//
//    int32_t pos = ctx.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
//    char *chr = ctx.bamHdr->target_name[ctx.aln->core.tid] ; //contig name (chromosome)
//    uint32_t len = ctx.aln->core.l_qseq; //length of the read.
//
//    uint8_t *q = bam_get_seq(ctx.aln); //quality string
//    uint32_t q2 = ctx.aln->core.qual ; //mapping quality
//
//
//    char *qseq = (char *)malloc(len);
//
//    for(int i=0; i< len ; i++){
//      qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
//    }
//
//    cout << chr << " " << pos << " " << len << " " << " " << qseq << " " << " " << q << " " << q2 << endl;
//
//    if(strcmp(chrom, chr) == 0){
//
//      if(locus > pos+len){
//        cout << chr << pos << len << qseq << q2;
//      }
//    }
//  }

