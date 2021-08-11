#include "sigmoid.h"

#include <math.h>

Sigmoid::Sigmoid(int logit_size) : logit_size_(logit_size),
    logit_table_(logit_size, 0) {
  for (int i = 0; i < logit_size_; ++i) {
    logit_table_[i] = SlowLogit((i + 0.5f) / logit_size_);
  }
}

float Sigmoid::Logit(float p) const {
  int index = p * logit_size_;
  if (index >= logit_size_) index = logit_size_ - 1;
  else if (index < 0) index = 0;
  return logit_table_[index];
}

float Sigmoid::Logistic(float p) {
  return 1 / (1 + exp(-p));
}

float Sigmoid::SlowLogit(float p) {
  return log(p / (1 - p));
}
