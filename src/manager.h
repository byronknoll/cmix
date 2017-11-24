#ifndef MANAGER_H
#define MANAGER_H

#include "states/nonstationary.h"
#include "states/run-map.h"
#include "contexts/context.h"
#include "contexts/bit-context.h"

#include <vector>
#include <memory>

class Manager {
 friend class Predictor;
 public:
  Manager();
  const Context& AddContext(std::unique_ptr<Context> context);
  const BitContext& AddBitContext(std::unique_ptr<BitContext> bit_context);
  void UpdateContexts(int bit);

 private:
  void UpdateHistory();
  void UpdateWords();
  void UpdateRecentBytes();

  unsigned int bit_context_;
  unsigned long long long_bit_context_, zero_context_, history_pos_,
      line_break_, longest_match_;
  std::vector<unsigned char> history_, shared_map_;
  std::vector<unsigned long long> words_, recent_bytes_;
  std::vector<std::unique_ptr<Context>> contexts_;
  std::vector<std::unique_ptr<BitContext>> bit_contexts_;
  RunMap run_map_;
  Nonstationary nonstationary_;
};

#endif

