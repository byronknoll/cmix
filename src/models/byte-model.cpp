#include "byte-model.h"

#include <numeric>

ByteModel::ByteModel() : top_(255), mid_(0), bot_(0), probs_(1.0 / 256, 256) {}

float ByteModel::Predict() {
  mid_ = bot_ + ((top_ - bot_) / 2);
  float num = std::accumulate(&probs_[mid_ + 1], &probs_[top_ + 1], 0.0f);
  float denom = std::accumulate(&probs_[bot_], &probs_[mid_ + 1], num);
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
