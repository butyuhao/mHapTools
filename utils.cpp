//
// Created by Yuhao Dan on 2020/9/13.
//

#include "./include/utils.h"

namespace std {
  bool is_suffix(string str, string suffix) {
    if (str.size() < suffix.size()) {
      return false;
    }
    if (str.size() == 0 || suffix.size() == 0) {
      return false;
    }
    for (int i = 1; i <= suffix.size(); i++) {
      if (str[str.size() - i] != suffix[suffix.size() - i]) {
        return false;
      }
    }
    return true;
  }
}
