#include <iostream>
#include <string.h>
#include "convert.cpp"

using namespace std;

int main(int argc, char *argv[]) {

  hts_set_log_level(HTS_LOG_WARNING);

  if (argc < 2) {
    cout << "See usage" << endl;
  }
  int ret = 0;
  if (strcmp(argv[1], "convert") == 0)    ret = main_convert(argc - 1, argv + 1);

  return ret;


}