#ifndef BYTE_MIXER_H
#define BYTE_MIXER_H

#include "../models/byte-model.h"
#include "logistic.h"

#include <valarray>

class ByteMixer : public ByteModel {
 public:
  ByteMixer(const Logistic& logistic, int input_neurons, int hidden_neurons,
      const unsigned int& bit_context, float learning_rate);
  void SetInput(int index, float val);
  void Train();
  void ByteUpdate();

 private:
  const unsigned int& byte_;
  const Logistic& logistic_;
  float learning_rate_;
  unsigned long long steps_;
  std::valarray<std::valarray<std::valarray<float>>> weights_;
  std::valarray<std::valarray<float>> states_, errors_;
  std::valarray<float> outputs_;
};

#endif
