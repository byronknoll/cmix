#ifndef MANAGER_H
#define MANAGER_H

#include "states/nonstationary.h"
#include "states/run-map.h"
#include "contexts/context.h"

#include <vector>
#include <memory>

class Manager {
 public:
  Manager(unsigned long long file_size);
  const Context& AddContext(std::unique_ptr<Context> context);
  void Perceive(int bit);

  unsigned int bit_context_;
  unsigned int zero_context_;
  unsigned long long history_pos_;
  std::vector<unsigned char> history_;
  std::vector<unsigned long long> words_;
  std::vector<unsigned long long> recent_bytes_;
  std::vector<unsigned int> recent_bytes2_;
  std::vector<std::unique_ptr<Context>> contexts_;
  RunMap run_map_;
  Nonstationary nonstationary_;

 private:
  void UpdateHistory();
  void UpdateWords();
  void UpdateRecentBytes();
};

#endif

