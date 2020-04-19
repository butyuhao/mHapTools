//
// Created by Yuhao Dan on 2020/4/13.
//

#ifndef BAM2HAP__CONVERT_H_
#define BAM2HAP__CONVERT_H_
#include <vector>
#include <string_view>
using namespace std;

class query_region {
  //用于保存与处理一条查询的染色体、开始位点、结束位点。
 public:
  query_region(char* &_arg)
      : arg (_arg) {};
  bool parse();
  void print();
 private:
  char *arg;
  string_view i_chr;
  int i_start = -1;
  int i_end = -1;
};

#endif //BAM2HAP__CONVERT_H_
