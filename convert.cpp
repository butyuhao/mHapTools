//
// Created by Yuhao Dan on 2020/4/13.
//
#include <getopt.h>
#include "convert.h"
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include "htslib-1.10.2/htslib/sam.h"

using namespace std;

//Configure the getopt

struct convert_args_s {
  samFile *fp_in;
  bam_hdr_t *bamHdr;
  bam1_t *aln;
  char *input_file;   /* -i option */
  char *aligner;      /* -a option */
  char *bed_file;     /* -b option */
  char *cpg_pos;      /* -c option */
  char *region;       /* -r option */
  char *output_file;  /* -o option */
} convert_args;

static const char *opt_string = "i:abc:ro";

static const struct option long_opts[] = {
    { "input", required_argument, NULL, 'i' },
    { "aligner", optional_argument, NULL, 'a' },
    { "bed_file", optional_argument, NULL, 'b' },
    { "cpg_pos", required_argument, NULL, 'c' },
    { "region", optional_argument, NULL, 'r' },
    { NULL, no_argument, NULL, 0 }
};

int main_convert(int argc, char *argv[]) {
  //Initialize convert_args
  convert_args.fp_in = NULL;
  convert_args.bamHdr = NULL;
  convert_args.aln = NULL;
  convert_args.input_file = NULL;
  convert_args.aligner = NULL;
  convert_args.bed_file = NULL;
  convert_args.cpg_pos = NULL;
  convert_args.region = NULL;
  convert_args.output_file = NULL;

  int long_index;

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while(opt != -1) {
    switch(opt) {
      case 'i':
        convert_args.fp_in = sam_open(optarg,"r");
        convert_args.bamHdr = sam_hdr_read(convert_args.fp_in); //read header
        convert_args.aln = bam_init1(); //initialize an alignment
        break;

      case 'a':
        cout << optarg << endl;
        break;

      case 'b':
        cout << optarg << endl;
        break;

      case 'c':
        cout << optarg << endl;
        break;

      case 'r':
        cout << optarg << endl;
        break;

      default:
        break;
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }
  //header parse
  char *tar = convert_args.bamHdr->text ;
  uint32_t *tarlen = convert_args.bamHdr->target_len ;

  cout << tar << endl;
  // content parse
  char *chrom = "1";
  int locus = 2000;
  int comp;

  while(sam_read1(convert_args.fp_in,convert_args.bamHdr,convert_args.aln) > 0){

    int32_t pos = convert_args.aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
    char *chr = convert_args.bamHdr->target_name[convert_args.aln->core.tid] ; //contig name (chromosome)
    uint32_t len = convert_args.aln->core.l_qseq; //length of the read.

    uint8_t *q = bam_get_seq(convert_args.aln); //quality string
    uint32_t q2 = convert_args.aln->core.qual ; //mapping quality


    char *qseq = (char *)malloc(len);

    for(int i=0; i< len ; i++){
      qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
    }

    cout << chr << " " << pos << " " << len << " " << " " << qseq << " " << " " << q << " " << q2 << endl;

    if(strcmp(chrom, chr) == 0){

      if(locus > pos+len){
        cout << chr << pos << len << qseq << q2;
      }
    }
  }

  bam_destroy1(convert_args.aln);
  sam_close(convert_args.fp_in);

  return EXIT_SUCCESS;
}