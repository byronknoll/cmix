#ifndef INDIRECT_HASH_H
#define INDIRECT_HASH_H

#include "context.h"

#include <vector>

class IndirectHash : public Context {
 public:
  IndirectHash(const unsigned int& bit_context, unsigned int order1,
      unsigned int hash_size1, unsigned int order2, unsigned int hash_size2);
  void Update();
  bool IsEqual(Context* c);

 private:
  const unsigned int& byte_;
  unsigned long long context1_;
  unsigned int hash_size1_, hash_size2_, size1_;
  std::vector<unsigned long long> hashes_;
};

#endif
