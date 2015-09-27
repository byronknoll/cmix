#ifndef BYTE_MIXER_H
#define BYTE_MIXER_H

#include "logistic.h"
#include "../models/byte-model.h"

#include <vector>
#include <valarray>

class ByteMixer : public ByteModel {
 public:
  ByteMixer(const unsigned int& bit_context, const std::valarray<float>& inputs,
      const Logistic& logistic, float learning_rate);
  void Mix();
  void ByteUpdate();

 private:
  const unsigned int& byte_;
  const std::valarray<float>& inputs_;
  const Logistic& logistic_;
  float learning_rate_;
  std::vector<std::valarray<float>> weights_;
};

#endif
