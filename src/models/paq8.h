#ifndef PAQ8_H
#define PAQ8_H

#include "model.h"

class PAQ8 : public Model {
 public:
  PAQ8();
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};
};

#endif
