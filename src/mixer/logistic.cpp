#include "logistic.h"

#include <math.h>

Logistic::Logistic(int stretch_size, int squash_size) :
    stretch_size_(stretch_size), squash_size_(squash_size),
    stretch_table_(stretch_size, 0), squash_table_(squash_size, 0) {
  for (int i = 0; i < stretch_size_; ++i) {
    stretch_table_[i] = SlowStretch((i + 0.5) / stretch_size_);
  }
  for (int i = 0; i < squash_size_; ++i) {
    squash_table_[i] = SlowSquash(100 * ((i + 0.5) / squash_size_) - 50);
  }
}

float Logistic::Stretch(float p) const {
  return stretch_table_[p * stretch_size_];
}

float Logistic::Squash(float p) const {
  int index = ((p + 50) / 100) * squash_size_;
  if (index >= squash_size_) index = squash_size_ - 1;
  else if (index < 0) index = 0;
  return squash_table_[index];
}

float Logistic::SlowStretch(float p) {
  return log(p / (1 - p));
}

float Logistic::SlowSquash(float p) {
  return 1 / (1 + exp(-p));
}
