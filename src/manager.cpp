#include "manager.h"

Manager::Manager() : bit_context_(1), long_bit_context_(1), zero_context_(0),
    history_pos_(0), line_break_(0), longest_match_(0), history_(100000000, 0),
    words_(8, 0), recent_bytes_(8, 0) {}

const Context& Manager::AddContext(std::unique_ptr<Context> context) {
  for (const auto& old : contexts_) {
    if (old->IsEqual(context.get())) return *old;
  }
  contexts_.push_back(std::move(context));
  return *(contexts_[contexts_.size() - 1]);
}

const BitContext& Manager::AddBitContext(std::unique_ptr<BitContext>
    bit_context) {
  bit_contexts_.push_back(std::move(bit_context));
  return *(bit_contexts_[bit_contexts_.size() - 1]);
}

void Manager::UpdateHistory() {
  history_[history_pos_] = bit_context_;
  ++history_pos_;
  if (history_pos_ == history_.size()) history_pos_ = 0;
}

void Manager::UpdateWords() {
  unsigned char c = bit_context_;
  if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c >= 0x80) {
    words_[7] = words_[7] * 997*16 + c;
  } else {
    words_[7] = 0;
  }
  if (c >= 'A' && c <= 'Z') c += 'a' - 'A';
  if ((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == 8 || c == 6 ||
      c >= 0x80) {
    words_[0] = words_[0] * 997*16 + c;
    words_[0] &= 0xfffffff;
    words_[1] = words_[1] * 263*32 + c;
  } else {
    for (int i = 6; i >= 2; --i) {
      words_[i] = words_[i-1];
    }
    words_[1] = 0;
  }
}

void Manager::UpdateRecentBytes() {
  for (int i = 7; i >= 1; --i) {
    recent_bytes_[i] = recent_bytes_[i-1];
  }
  recent_bytes_[0] = bit_context_;
}

void Manager::Perceive(int bit) {
  bit_context_ += bit_context_ + bit;
  long_bit_context_ = bit_context_;
  if (bit_context_ >= 256) {
    bit_context_ -= 256;
    long_bit_context_ = 1;
    longest_match_ = 0;

    if (bit_context_ == '\n') {
      line_break_ = 0;
    } else if (line_break_ < 99) {
      ++line_break_;
    }

    UpdateHistory();
    UpdateWords();
    UpdateRecentBytes();
    for (const auto& context : contexts_) {
      context->Update();
    }
  }
  for (const auto& context : bit_contexts_) {
    context->Update();
  }
}
