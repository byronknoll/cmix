#ifndef MODEL_H
#define MODEL_H

#include <valarray>

class Model {
 public:
  Model() : outputs_(0.5, 1) {}
  virtual ~Model() {}
  virtual const std::valarray<float>& Predict() {return outputs_;}
  virtual unsigned int NumOutputs() {return 1;}
  virtual void Perceive(int bit) {}
  virtual void ByteUpdate() {}

 protected:
  std::valarray<float> outputs_;
};

#endif
