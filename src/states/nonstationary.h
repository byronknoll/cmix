#ifndef NONSTATIONARY_H
#define NONSTATIONARY_H

#include "state.h"

#include <array>

class Nonstationary : public State {
 public:
  Nonstationary();
  float InitProbability(int state) const;
  int Next(int state, int bit) const;

 private:
  std::array<std::array<unsigned char, 2>, 256> table_;
};

#endif

