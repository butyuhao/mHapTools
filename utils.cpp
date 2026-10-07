//
// Created by Yuhao Dan on 2020/9/13.
//

#include <vector>
#include "./include/utils.h"
#include "./htslib-1.10.2/htslib/hts.h"

namespace std {
  // Index of the first CpG position >= cpg_pos (binary search); used by merge
  // and beta. Moved here unchanged from the former convert.cpp.
  int _lower_bound(vector<hts_pos_t> &v, hts_pos_t &cpg_pos) {
    int low = 0, high = v.size() - 1;
    while (low < high) {
      int mid = low + (high - low) / 2;
      if (v[mid] >= cpg_pos) high = mid;
      else low = mid + 1;
    }
    return low;
  }

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
