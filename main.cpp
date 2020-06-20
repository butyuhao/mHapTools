#include <iostream>
#include <string.h>
#include "convert.cpp"
#include <time.h>
using namespace std;

int main(int argc, char *argv[]) {

  hts_set_log_level(HTS_LOG_INFO);

  if (argc < 2) {
    cout << "See usage" << endl;
  }
  int ret = 0;

  time_t first, second;
  first=time(NULL);

  if (strcmp(argv[1], "convert") == 0)    ret = main_convert(argc - 1, argv + 1);

  second=time(NULL);

  hts_log_info("The duration is: %f seconds", difftime(second,first));

  cout << "Process finished" << endl;

  return ret;

}