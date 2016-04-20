#ifndef SPARSE_H
#define SPARSE_H

#include "context.h"

#include <vector>

class Sparse : public Context {
 public:
  Sparse(const std::vector<unsigned long long>& recent_contexts,
      const std::vector<unsigned int>& orders);
  void Update();
  bool IsEqual(Context* c);

 private:
  const std::vector<unsigned long long>& recent_contexts_;
  std::vector<unsigned int> orders_, factors_;
};

#endif
