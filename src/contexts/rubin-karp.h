#ifndef RUBIN_KARP_H
#define RUBIN_KARP_H

#include "context.h"
#include <vector>

class RubinKarp : public Context {
 public:
  RubinKarp(const unsigned int& bit_context, unsigned int order);
  void Update();
  bool IsEqual(Context* c);

  const unsigned int& byte_;
  unsigned int order_;
  unsigned long long pow_;
  unsigned int base_;
  std::vector<unsigned char> rolling_;
  unsigned int pos_;
};

#endif
