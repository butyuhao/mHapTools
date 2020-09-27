#include <iostream>
#include <string.h>
#include <time.h>
#include "./include/merge.h"
#include "./include/convert.h"
#include "./include/beta.h"
#include "./include/summary.h"
#include "./include/version.h"
#include "./include/utils.h"

using namespace std;

static void help() {
  cout << endl;
  cout << "mHapTools " << HAP_VERSION_TEXT <<  " (Tools for analysing methylated reads)" << endl;
  cout << endl;
  cout << "Usage:    mhaptools <command> [options]" << endl;
  cout << endl;
  cout << "Commands:" << endl;
  cout << "    convert     SAM/BAM --> mhap conversion" << endl;
  cout << "    merge       merge two mhap files" << endl;
  cout << "    beta        count methylated reads on CpG positions" << endl;
  cout << "    summary     get local or global Summarized Information" << endl;
}

int main(int argc, char *argv[]) {
  
  hts_set_log_level(HTS_LOG_WARNING);
  //hts_set_log_level(HTS_LOG_TRACE);

  if (argc < 2) {
    help();
    return 0;
  }

  int ret = 0;

  time_t start, stop;
  start=time(NULL);

  if (strcmp(argv[1], "convert") == 0) {
    ret = main_convert(argc - 1, argv + 1);

  } else if (strcmp(argv[1], "merge") == 0) {
    ret = main_merge(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "beta") == 0) {
    ret = main_beta(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "summary") == 0) {
    ret = main_summary(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "help") == 0) {
    help();
  } else if (strcmp(argv[1], "--help") == 0) {
    help();
  } else if (strcmp(argv[1], "-h") == 0) {
    help();
  } else {
    cout << "unrecognized command " <<  argv[1] << endl;
    return 0;
  }

  stop=time(NULL);

  cout << "Process finished." << endl;

  cout << "The duration is: "<< difftime(stop,start) <<" seconds." << endl;

  return ret;
}

