#include "logistic.h"

#include <math.h>

Logistic::Logistic(int stretch_size) : stretch_size_(stretch_size),
    stretch_table_(stretch_size, 0) {
  for (int i = 0; i < stretch_size_; ++i) {
    stretch_table_[i] = SlowStretch((i + 0.5) / stretch_size_);
  }
}

float Logistic::Stretch(float p) const {
  return stretch_table_[p * stretch_size_];
}

float Logistic::Squash(float p) const {
  return 1 / (1 + exp(-p));
}

float Logistic::SlowStretch(float p) {
  return log(p / (1 - p));
}
