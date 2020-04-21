//
// Created by Yuhao Dan on 2020/4/13.
//
#include <getopt.h>
#include "convert.h"
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/hts.h>

using namespace std;

//Configure the getopt

struct context {
  htsFile *fp_bam;
  htsFile *fp_cpg;
  bam_hdr_t *hdr_bam;/* -i bam文件头部的指针 */

  bam1_t *aln;

  char *fn_bam;
  char *fn_cpg;
  char *bam_path;   /* -i option */
  char *output_path;  /* -o option */
  char *aligner;      /* -a option */
  char *bed_file;     /* -b option */
  char *cpg_path;      /* -c option */
  char *region;       /* -r option */
} ctx;

bool query_region::parse() {
  string_view arg_sv = string_view (arg);
  bool is_get_chr = false;
  bool is_get_start = false;
  bool is_get_end = false;
  int i;
  for(i = 0; i < arg_sv.size(); i++) {
    if(arg_sv[i] == ':') {
      //get the chr name
      i_chr = arg_sv.substr(0, i);
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

void query_region::print() {
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

void parse_sam(context &ctx, query_region &q_region) {
  
}

inline void init_ctx() {
  //Initialize ctx
  ctx.fp_bam = NULL;
  ctx.fp_cpg = NULL;
  ctx.hdr_bam = NULL;/* -i bam文件头部的指针 */

  ctx.aln = NULL;

  ctx.fn_bam = NULL;
  ctx.fn_cpg = NULL;
  ctx.bam_path = NULL;   /* -i option */
  ctx.cpg_path = NULL;      /* -c option */
  ctx.output_path = NULL;  /* -o option */
  ctx.aligner = NULL;      /* -a option */
  ctx.bed_file = NULL;     /* -b option */

  ctx.region = NULL;       /* -r option */
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

void destroy(context &ctx) {
  if (ctx.fp_bam) {
    hts_close(ctx.fp_bam);
  }
  if(ctx.fn_cpg) {
    hts_close(ctx.fp_cpg);
  }
  if(ctx.aln) {
    bam_destroy1(ctx.aln);
  }
}

int main_convert(int argc, char *argv[]) {
  init_ctx();

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
    query_region q_region = query_region(ctx.region);
    ret = q_region.parse();
    if(ret == EXIT_FAILURE) {
      //fail to parse the query
      q_region.print();
      destroy(ctx);
      return EXIT_FAILURE;
    }
    ret = open_bam_file(ctx);
    if(ret == EXIT_FAILURE) {
      destroy(ctx);
      return EXIT_FAILURE;
    }
    ret = open_cpg_file(ctx);
    if(ret == EXIT_FAILURE) {
      destroy(ctx);
      return EXIT_FAILURE;
    }
    parse_sam(ctx, q_region);
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

