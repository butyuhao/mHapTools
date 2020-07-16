#include <iostream>
#include <string.h>
#include <time.h>
#include "merge.cpp"
#include "convert.cpp"

using namespace std;

int main(int argc, char *argv[]) {

  hapFile *fp = hap_open("/Users/butyuhao/Documents/GitHub/haptools/out.hap", "r");
  
  hap_line_t h_line_t = {" ", 0, 0, " ", 0, HAP_DEFAULT_DIRECTION};
  
  while(hap_read(fp, &h_line_t) == 0) {
    h_line_t.print();
  }

  hts_set_log_level(HTS_LOG_TRACE);

  if (argc < 2) {
    cout << "See usage:" << endl;
    cout << "-i path of the bam file" << endl;
    cout << "-a aligner [BISMARK | BSMAP | UNKNOWN]" << endl;
    cout << "-b path of the bed file" << endl;
    cout << "-c path of the cpg file" << endl;
    cout << "-r region" << endl;
    cout << "-o output path" << endl;
    cout << "-i path of the bam file" << endl;
    exit(1);
  }
  int ret = 0;

  time_t start, stop;
  start=time(NULL);

  if (strcmp(argv[1], "convert") == 0) { 
    ret = main_convert(argc - 1, argv + 1); 

  } else if (strcmp(argv[1], "merge") == 0) {
    ret = main_merge(argc - 1, argv + 1);
  }

  stop=time(NULL);

  hts_log_info("The duration is: %f seconds", difftime(stop,start));

  cout << "Process finished" << endl;

  return ret;

}