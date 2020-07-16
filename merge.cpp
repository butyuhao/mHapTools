#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include "hap.cpp"

namespace std {
  int main_merge(int argc, char *argv[]) {
      int long_index;

  static const char *opt_string = "i:o:";

  static const struct option long_opts[] = {
      { "input", required_argument, NULL, 'i' },
      { "output", optional_argument, NULL, 'o' },
      { NULL, no_argument, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  while (opt != -1) {
    switch (opt) {
      case 'i': {
        cout << optarg << endl;
        if (argc < 4) {
          return 1;
        } else {
          cout << argv[optind] << endl;
        }
        
        break;
      }
      case 'o': {
        break;
      }
      default: {
        break;
      }
    }
    opt = getopt_long(argc, argv, opt_string, long_opts, &long_index);
  }
  return 0;
}
}//namespace std
