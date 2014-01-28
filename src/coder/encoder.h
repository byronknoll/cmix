#ifndef ENCODER_H
#define ENCODER_H

#include <fstream>

#include "../predictor.h"

class Encoder {
 public:
  Encoder(std::ofstream* os, unsigned long long file_size);
  void Encode(int bit);
  void Flush();

 private:
  void WriteByte(unsigned int byte);

  std::ofstream* os_;
  unsigned int x1_, x2_;
  Predictor p_;
};

#endif
