#ifndef CONTEXT_MANAGER_H
#define CONTEXT_MANAGER_H

#include "states/nonstationary.h"
#include "states/run-map.h"
#include "contexts/context.h"
#include "contexts/bit-context.h"

#include <vector>
#include <memory>

struct ContextManager {
  ContextManager();
  const Context& AddContext(std::unique_ptr<Context> context);
  const BitContext& AddBitContext(std::unique_ptr<BitContext> bit_context);
  void UpdateContexts(int bit);
  void UpdateHistory();
  void UpdateWords();
  void UpdateRecentBytes();
  void UpdateWRTContext();

  unsigned int bit_context_ = 1, wrt_state_ = 0;
  unsigned long long long_bit_context_ = 1, zero_context_ = 0, history_pos_ = 0,
      line_break_ = 0, longest_match_ = 0, auxiliary_context_ = 0,
      wrt_context_ = 0;
  std::vector<unsigned char> history_, shared_map_;
  std::vector<unsigned long long> words_, recent_bytes_;
  std::vector<std::unique_ptr<Context>> contexts_;
  std::vector<std::unique_ptr<BitContext>> bit_contexts_;
  RunMap run_map_;
  Nonstationary nonstationary_;
};

#endif
