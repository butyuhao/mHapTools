//
// Created by Yuhao Dan on 2020/4/13.
//
#include <getopt.h>
#include "convert.h"
#include <unistd.h>
#include <getopt.h>

using namespace std;

//Configure the getopt

struct convert_args_s {
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
        cout << optarg << endl;
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
  return EXIT_SUCCESS;
}