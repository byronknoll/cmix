#ifndef BYTE_RUN_H
#define BYTE_RUN_H

#include "model.h"

#include <vector>
#include <array>

class ByteRun : public Model {
 public:
  ByteRun(const unsigned long long& byte_context,
      const unsigned int& bit_context, float delta,
      unsigned long long map_size);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate();

 private:
  const unsigned long long& byte_context_;
  const unsigned int& bit_context_;
  unsigned char byte_prediction_;
  unsigned char run_length_;
  unsigned int map_index_;
  unsigned char bit_pos_;
  float divisor_;
  std::vector<unsigned char> map_;
  std::vector<unsigned char> counts_;
  std::array<float, 256> predictions_;
};

#endif
