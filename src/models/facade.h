#ifndef FACADE_H
#define FACADE_H

#include "model.h"

class Facade : public Model {
 public:
  Facade(const float& output) : output_(output) {};
  float Predict() { return output_; };
  void Perceive(int bit) {};
  void ByteUpdate() {};

 private:
  const float& output_;
};

#endif
