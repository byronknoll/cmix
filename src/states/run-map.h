#ifndef RUN_MAP_H
#define RUN_MAP_H

#include "state.h"

#include <array>

class RunMap : public State {
 public:
  RunMap();
  float InitProbability(int state) const;
  int Next(int state, int bit) const;

 private:
  std::array<unsigned char, 512> table_;
};

#endif
