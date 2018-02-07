#ifndef BIT_CONTEXT_H
#define BIT_CONTEXT_H

#include "context.h"

class BitContext : public Context {
 public:
  BitContext(const unsigned long long& bit_context, const unsigned long long&
      byte_context, unsigned long long byte_context_size);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned long long& bit_context_;
  const unsigned long long& byte_context_;
};

#endif
