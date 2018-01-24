#ifndef INTERVAL_HASH_H
#define INTERVAL_HASH_H

#include "context.h"

#include <vector>

class IntervalHash : public Context {
 public:
  IntervalHash(const unsigned int& bit_context, const std::vector<int>& map,
      unsigned int num_bits, unsigned int order, unsigned int hash_size);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned int& byte_;
  std::vector<int> map_;
  unsigned long long mask_;
  unsigned int hash_size_, interval_, shift_;
};

#endif
