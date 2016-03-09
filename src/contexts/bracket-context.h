#ifndef BRACKET_CONTEXT_H
#define BRACKET_CONTEXT_H

#include "context.h"

#include <unordered_map>
#include <vector>

class BracketContext : public Context {
 public:
  BracketContext(const unsigned int& bit_context);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned int& byte_;
  std::unordered_map<unsigned char, unsigned char> brackets_;
  int max_;
  std::vector<int> active_, distance_;
};

#endif
