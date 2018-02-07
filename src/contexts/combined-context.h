#ifndef COMBINED_CONTEXT_H
#define COMBINED_CONTEXT_H

#include "context.h"

class CombinedContext : public Context {
 public:
  CombinedContext(const unsigned long long& context1, const unsigned long long&
      context2, unsigned long long context1_size, unsigned long long
      context2_size);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned long long& context1_;
  const unsigned long long& context2_;
  unsigned int shift_;
};

#endif
