#include "byte-mixer.h"

#include <numeric>

ByteMixer::ByteMixer(const unsigned int& bit_context,
    const std::valarray<float>& inputs, const Logistic& logistic,
    float learning_rate) : byte_(bit_context), inputs_(inputs),
    logistic_(logistic), learning_rate_(learning_rate),
    weights_(256, std::valarray<float>(0.0, inputs.size())) {}

void ByteMixer::Mix() {
  for (int i = 0; i < 256; ++i) {
    probs_[i] = logistic_.Squash(std::inner_product(&inputs_[0],
        &inputs_[inputs_.size()], &weights_[i][0], 0.0));
  }
  ByteModel::ByteUpdate();
}

void ByteMixer::ByteUpdate() {
  for (unsigned int i = 0; i < 256; ++i) {
    float bit = 0;
    if (i == byte_) {
      bit = 1;
    }
    float update = learning_rate_ * (bit - probs_[i]);
    weights_[i] += update * inputs_;
  }
}

