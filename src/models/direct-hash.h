#ifndef DIRECT_HASH_H
#define DIRECT_HASH_H

#include "model.h"

#include <vector>
#include <array>

class DirectHash : public Model {
 public:
  DirectHash(const unsigned long long& byte_context,
      const unsigned int& bit_context, int limit, float delta, int size);
  const std::valarray<float>& Predict();
  void Perceive(int bit);
  void ByteUpdate();

 private:
  const unsigned long long& byte_context_;
  const unsigned int& bit_context_;
  unsigned long long index_;
  int limit_;
  float delta_, divisor_;
  std::vector<std::array<float, 256>> predictions_;
  std::vector<std::array<unsigned char, 256>> counts_;
  std::vector<unsigned long long> checksums_;
};

#endif
