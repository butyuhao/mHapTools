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

using namespace std;

//Configure the getopt

struct context {
  samFile *fp_in;
  bam_hdr_t *header;
  bam1_t *aln;
  char *input_file_path;   /* -i option */
  char *output_file_path;  /* -o option */
  char *aligner;      /* -a option */
  char *bed_file;     /* -b option */
  char *cpg_file;      /* -c option */
  char *region;       /* -r option */
  tbx_t *tbx_idx;
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
    { "cpg_file", required_argument, NULL, 'c' },
    { "region", optional_argument, NULL, 'r' },
    { NULL, no_argument, NULL, 0 }
};

bool open_tabix_file(context &ctx) {
  ctx.tbx_idx = tbx_index_load(ctx.cpg_file);
  cout << ctx.tbx_idx->idx;

}

void parse_sam(context &ctx, query_region &q_region) {
  open_tabix_file(ctx);
}

inline void init_ctx() {
  //Initialize ctx
  ctx.fp_in = NULL;
  ctx.header = NULL;
  ctx.aln = NULL;
  ctx.input_file_path = NULL;
  ctx.aligner = NULL;
  ctx.bed_file = NULL;
  ctx.cpg_file = NULL;
  ctx.region = NULL;
  ctx.output_file_path = NULL;
  ctx.tbx_idx = NULL;
}

bool open_sam_file(context &ctx) {

  ctx.fp_in = hts_open(ctx.input_file_path, "r");
  ctx.header = sam_hdr_read(ctx.fp_in); //read header
  if(ctx.header == NULL) {
    sam_close(ctx.fp_in);
    return EXIT_FAILURE;
  }
  //ctx.aln = bam_init1(); //initialize an alignment
  return EXIT_SUCCESS;
}

int main_convert(int argc, char *argv[]) {
  init_ctx();

  int long_index;

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while(opt != -1) {
    switch(opt) {
      case 'i':
        ctx.input_file_path = optarg;
        break;

      case 'a':
        cout << optarg << endl;
        break;

      case 'b':
        cout << optarg << endl;
        break;

      case 'c':
        ctx.cpg_file = optarg;
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
      goto DTOR;
    }
    ret = open_sam_file(ctx);
    parse_sam(ctx, q_region);
  }

  DTOR:
  bam_destroy1(ctx.aln);
  sam_close(ctx.fp_in);

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
//  while(sam_read1(ctx.fp_in,ctx.bamHdr,ctx.aln) > 0){
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

