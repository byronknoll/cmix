#include "byte-model.h"

ByteModel::ByteModel() : top_(255), mid_(0), bot_(0), probs_(1.0 / 256, 256) {}

float ByteModel::Predict() {
  float num = 0, denom = 0;
  mid_ = bot_ + ((top_ - bot_) / 2);
  for (int i = bot_; i <= top_; ++i) {
    denom += probs_[i];
    if (i > mid_) num += probs_[i];
  }
  if (denom == 0) return 0.5;
  return num / denom;
}

const std::valarray<float>& ByteModel::BytePredict() {
  return probs_;
}

void ByteModel::Perceive(int bit) {
  if (bit) {
    bot_ = mid_ + 1;
  } else {
    top_ = mid_;
  }
}

void ByteModel::ByteUpdate() {
  top_ = 255;
  bot_ = 0;
}
