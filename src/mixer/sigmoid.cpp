#include "sigmoid.h"

#include <math.h>

Sigmoid::Sigmoid(int logit_size) : logit_size_(logit_size),
    logit_table_(logit_size, 0) {
  for (int i = 0; i < logit_size_; ++i) {
    logit_table_[i] = SlowLogit((i + 0.5) / logit_size_);
  }
}

float Sigmoid::Logit(float p) const {
  return logit_table_[p * logit_size_];
}

float Sigmoid::Logistic(float p) {
  return 1 / (1 + exp(-p));
}

float Sigmoid::SlowLogit(float p) {
  return log(p / (1 - p));
}
