#ifndef DMC_H
#define DMC_H

#include "model.h"

#include <vector>
#include <array>

struct Node {
  std::array<float, 2> count;
  std::array<unsigned int, 2> next;
};

class DMC : public Model {
 public:
  DMC(float delta, unsigned int max_size);
  const std::valarray<float>& Predict();
  void Perceive(int bit);
  void ByteUpdate() {};

 private:
  void Reset();

  std::vector<Node> nodes_;
  unsigned int cur_, max_size_;
  float delta_;
  int threshold_, big_threshold_;
};

#endif
