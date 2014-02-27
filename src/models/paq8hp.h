#ifndef PAQ8HP_H
#define PAQ8HP_H

#include "model.h"

class PAQ8HP : public Model {
 public:
  PAQ8HP();
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};
};

#endif
