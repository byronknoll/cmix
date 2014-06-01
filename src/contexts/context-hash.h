#ifndef CONTEXT_HASH_H
#define CONTEXT_HASH_H

#include "context.h"

class ContextHash : public Context {
 public:
  ContextHash(const unsigned int& bit_context, unsigned int order,
      unsigned int hash_size);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned int& byte_;
  unsigned int hash_size_;
};

#endif
