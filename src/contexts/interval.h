#ifndef INTERVAL_H
#define INTERVAL_H

#include "context.h"

#include <vector>

class Interval : public Context {
 public:
  Interval(const unsigned int& bit_context, const std::vector<int>& map,
      unsigned int num_bits);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned int& byte_;
  std::vector<int> map_;
  unsigned long long mask_;
  int shift_;
};

#endif
