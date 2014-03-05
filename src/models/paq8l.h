#ifndef PAQ8L_H
#define PAQ8L_H

#include "model.h"

class PAQ8L : public Model {
 public:
  PAQ8L(int memory);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};
};

#endif
