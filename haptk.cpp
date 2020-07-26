#include <iostream>
#include <string.h>
#include <time.h>
#include "./include/merge.h"
#include "./include/convert.h"
#include "./include/beta.h"


using namespace std;

int main(int argc, char *argv[]) {

  hts_set_log_level(HTS_LOG_ERROR);

  if (argc < 2) {
    help:
    cout << "Command: convert" << endl;
    cout << "Convert SAM/BAM file to hap file." << endl;
    cout << "    Usage:   haptools convert  [options] -i test.bam -c hg38_CpG.gz [-a BSMAP] [-o out.hap]" << endl;
    cout << "Options:" << endl;
    cout << "required:" << endl;
    cout << "    -c  <str>   CpG position file, with the corresponding .tbi file in the same folder. Prefix of the index file’s name should be kept same as the corresponding CpG file." << endl;
    cout << "    -i  <str>   input file, SAM/BAM format, should be sorted by samtools. with the corresponding .bai file in the same folder. Prefix of the index file’s name should be the same as the corresponding SAM/BAM file." << endl;

    cout << "    optional:" << endl;
    cout << "    -a  <str>   aligner, BSMAP/BISMARK/UNKNOWN, UNKNOWN (default)." << endl;
    cout << "    -o  <str>   output filename, out.hap is the default output filename." << endl;
    cout << "    -r  <str>   query region, e.g. chr1:2000-20000, default is to query the whole genome." << endl;
    cout << "    -b  <str>   bed file of query regions." << endl;
    cout << "    -h  <str>   help." << endl;

    cout << "Command: merge" << endl;
    cout << "Convert SAM/BAM file to hap file." << endl;
    cout << "Usage:   haptools merge  [options] -i in1.hap in2.hap -c hg38_CpG.gz [-o out.hap]" << endl;
    cout << "    required:" << endl;
    cout << "    -i  <str>   path of two hap files." << endl;
    cout << "    -c  <str>   CpG position file." << endl;

    cout << "    optional:" << endl;
    cout << "    -o  <str>   output filename, out.hap is the default output filename." << endl;
    exit(1);
  }
  int ret = 0;

  time_t start, stop;
  start=time(NULL);

  if (strcmp(argv[1], "convert") == 0) { 
    ret = main_convert(argc - 1, argv + 1); 

  } else if (strcmp(argv[1], "merge") == 0) {
    ret = main_merge(argc - 1, argv + 1);
  }  else if (strcmp(argv[1], "beta") == 0) {
    ret = main_beta(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "help") == 0) {
    goto help;
  } else if (strcmp(argv[1], "--help") == 0) {
    goto help;
  } else if (strcmp(argv[1], "-h") == 0) {
    goto help;
  }

  stop=time(NULL);

  cout << "Process finished." << endl;

  hts_log_info("The duration is: %f seconds", difftime(stop,start));

  cout << "The duration is: "<< difftime(stop,start) <<" seconds." << endl;

  return ret;

}