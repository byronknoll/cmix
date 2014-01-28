#ifndef BIT_BUFFER_H
#define BIT_BUFFER_H

#include "model.h"
#include "../states/state.h"

#include <vector>
#include <array>

class BitBuffer : public Model {
 public:
  BitBuffer(int size, float delta, int limit);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {}

 private:
  std::vector<bool> bit_buffer_;
  unsigned int bit_pos_;
  float delta_;
  unsigned int limit_;
  float divisor_;
  std::array<float, 2> predictions_;
  std::array<unsigned int, 2> counts_;
};

#endif
