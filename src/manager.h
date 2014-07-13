#ifndef MANAGER_H
#define MANAGER_H

#include "states/nonstationary.h"
#include "states/run-map.h"
#include "contexts/context.h"
#include "contexts/bit-context.h"

#include <vector>
#include <memory>

class Manager {
 public:
  Manager();
  const Context& AddContext(std::unique_ptr<Context> context);
  const BitContext& AddBitContext(std::unique_ptr<BitContext> bit_context);
  void Perceive(int bit);

  unsigned int bit_context_;
  unsigned long long long_bit_context_;
  unsigned long long zero_context_;
  unsigned long long history_pos_;
  unsigned long long line_break_;
  unsigned long long longest_match_;
  std::vector<unsigned char> history_;
  std::vector<unsigned long long> words_;
  std::vector<unsigned long long> recent_bytes_;
  std::vector<std::unique_ptr<Context>> contexts_;
  std::vector<std::unique_ptr<BitContext>> bit_contexts_;
  RunMap run_map_;
  Nonstationary nonstationary_;

 private:
  void UpdateHistory();
  void UpdateWords();
  void UpdateRecentBytes();
};

#endif

