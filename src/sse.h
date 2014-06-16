#ifndef SSE_H
#define SSE_H

#include "mixer/logistic.h"

#include <vector>
#include <array>

class SSE {
 public:
  SSE(const Logistic& logistic, const unsigned long long& byte_context,
      const unsigned int& bit_context, unsigned int num_buckets, float delta,
      unsigned long long table_size);
  float Process(float input);
  void Perceive(int bit);

 private:
  const unsigned int& bit_context_;
  const unsigned long long& byte_context_;
  unsigned int num_buckets_;
  unsigned int bucket_;
  float divisor_;
  std::vector<std::vector<std::array<float, 256>>> predictions_;

};

#endif
