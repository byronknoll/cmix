#ifndef DMC_H
#define DMC_H

#include "model.h"

#include <vector>

struct Node {
  float count[2];
  unsigned int next[2];
};

class DMC : public Model {
 public:
  DMC(float delta, unsigned int max_size);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};

 private:
  void Reset();

  std::vector<Node> nodes_;
  unsigned int cur_;
  float delta_;
  unsigned int max_size_;
  int threshold_, big_threshold_;
};

#endif
