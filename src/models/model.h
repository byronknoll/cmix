#ifndef MODEL_H
#define MODEL_H

class Model {
 public:
  virtual ~Model() {}
  virtual float Predict() {return 0.5;}
  virtual void Perceive(int bit) {}
  virtual void ByteUpdate() {}
};

#endif
