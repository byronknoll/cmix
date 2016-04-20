#ifndef LOGISTIC_H
#define LOGISTIC_H

#include <vector>

class Logistic {
 public:
  Logistic(int stretch_size, int squash_size);
  float Stretch(float p) const;
  float Squash(float p) const;

 private:
  float SlowStretch(float p);
  float SlowSquash(float p);
  int stretch_size_, squash_size_;
  std::vector<float> stretch_table_, squash_table_;
};

#endif
