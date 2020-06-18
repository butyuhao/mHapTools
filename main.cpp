#include <iostream>
#include <string.h>
#include "convert.cpp"
#include <time.h>
#include <chrono>
using namespace std;


int lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos)//二分查找求下界
//找到cpg_pos这个数字在v中的下标

{
  int low = 0, high = v.size() - 1;
  while(low < high)
  {
    int mid = low + (high - low)/2;
    if(v[mid] >= cpg_pos) high = mid;
    else low = mid + 1;
  }
  return low;
}

int main(int argc, char *argv[]) {
//
//  hts_set_log_level(HTS_LOG_INFO);
//
//  if (argc < 2) {
//    cout << "See usage" << endl;
//  }
//  int ret = 0;
//
//  time_t first, second;
//  first=time(NULL);
//
//  if (strcmp(argv[1], "convert") == 0)    ret = main_convert(argc - 1, argv + 1);
//
//  second=time(NULL);
//
//  hts_log_info("The duration is: %f seconds", difftime(second,first));
//
//  return ret;
int test_n = 0;


htsFile *b = hts_open("/Users/butyuhao/Documents/GitHub/haptools/data/CDX188_COLFR0468_tumor_S1_Rmdup.bam", "r");

htsFile *f = hts_open("/Users/butyuhao/Documents/GitHub/haptools/data/hg19_CpG.gz", "r");
kstring_t cpg_line = {0,0,NULL};
unordered_map<int, vector<hts_pos_t>> cpg_pos_map;
unordered_map<int, vector<hts_pos_t>>::iterator cpg_pos_map_itor;
while (hts_getline(f, KS_SEP_LINE, &cpg_line) > 0) {

  char *p ,*q;
  int tid;
  hts_pos_t cpg_start;
  p = q = cpg_line.s;
  while(*q && *q != '\t') {
    q++;
  }
  *q = '\0';
  tid = bam_name2id(b->bam_header, p);
  *q = '\t';
  p = q + 1;
  q = p;
  while(*q && *q != '\t') {
    q++;
  }
  *q = '\0';
  cpg_start = atoll(p);

  cpg_pos_map_itor = cpg_pos_map.find(tid);

  if (cpg_pos_map_itor == cpg_pos_map.end()) {
    vector<hts_pos_t> v;
    v.push_back(cpg_start);
    cpg_pos_map[tid] = v;
  } else {
    cpg_pos_map[tid].push_back(cpg_start);
  }
  int pos = lower_bound(cpg_pos_map[tid], cpg_start);
  test_n++;

}

}


