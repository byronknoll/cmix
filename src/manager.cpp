#include "manager.h"

Manager::Manager() : bit_context_(1), zero_context_(0), history_pos_(0),
    line_break_(0), recent_bits_pos_(0), history_(100000000, 0), words_(8, 0),
    recent_bytes_(8, 0), recent_bytes2_(4, 0), pic_context_(6, 0),
    recent_bits_(216 * 8 * 7, false) {}

const Context& Manager::AddContext(std::unique_ptr<Context> context) {
  for (const auto& old : contexts_) {
    if (old->IsEqual(context.get())) return *old;
  }
  contexts_.push_back(std::move(context));
  return *(contexts_[contexts_.size() - 1]);
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

bool Manager::Bit(int index) {
  index = recent_bits_pos_ - index;
  if (index < 0) index += recent_bits_.size();
  return recent_bits_[index];
}

void Manager::UpdatePicContext(int bit) {
  recent_bits_[recent_bits_pos_] = bit;
  ++recent_bits_pos_;
  if (recent_bits_pos_ == recent_bits_.size()) recent_bits_pos_ = 0;

  pic_context_[0] = (Bit(216 * 8 * 2 + 1) << 7) +
      (Bit(216 * 8 * 2 - 1) << 6) +
      (Bit(216 * 8 * 2) << 5) +
      (Bit(216 * 8 - 2) << 4) +
      (Bit(216 * 8 + 1) << 3) +
      (Bit(216 * 8) << 2) +
      (Bit(216 * 8 - 1) << 1) + 
      Bit(1);

  pic_context_[1] = (Bit(216 * 8 - 1) << 7) +
      (Bit(4) << 6) +
      (Bit(3) << 5) +
      (Bit(2) << 4) +
      (Bit(216 * 8 * 4) << 3) +
      (Bit(216 * 8 * 3) << 2) +
      (Bit(216 * 8 * 2) << 1) +
      Bit(216 * 8);

  pic_context_[2] = (Bit(3) << 7) +
      (Bit(2) << 6) +
      (Bit(1) << 5) +
      (Bit(216 * 8 + 3) << 4) +
      (Bit(216 * 8 + 2) << 3) +
      (Bit(216 * 8 + 1) << 2) +
      (Bit(216 * 8) << 1) +
      Bit(216 * 8 - 1);

  pic_context_[3] = (Bit(216 * 8 + 1) << 1) + Bit(216 * 8 - 1);

  pic_context_[4] = (Bit(216 * 8 * 6) << 5) +
      (Bit(216 * 8 * 5) << 4) +
      (Bit(216 * 8 * 4) << 3) +
      (Bit(216 * 8 * 3) << 2) +
      (Bit(216 * 8 * 2) << 1) +
      Bit(216 * 8);

  pic_context_[5] = (Bit(8) << 7) + (Bit(7) << 6) + (Bit(6) << 5) +
      (Bit(5) << 4) + (Bit(4) << 3) + (Bit(3) << 2) + (Bit(2) << 1) + Bit(1);
}

void Manager::UpdateRecentBytes() {
  for (int i = 7; i >= 1; --i) {
    recent_bytes_[i] = recent_bytes_[i-1];
  }
  recent_bytes_[0] = bit_context_;
  for (int i = 3; i >= 1; --i) {
    recent_bytes2_[i] = recent_bytes2_[i-1];
  }
  recent_bytes2_[0] = bit_context_;
}

void Manager::Perceive(int bit) {
  UpdatePicContext(bit);
  bit_context_ += bit_context_ + bit;
  if (bit_context_ >= 256) {
    bit_context_ -= 256;

    if (bit_context_ == '\n') {
      line_break_ = 0;
    } else if (line_break_ < 255) {
      ++line_break_;
    }

    UpdateHistory();
    UpdateWords();
    UpdateRecentBytes();
    for (const auto& context : contexts_) {
      context->Update();
    }
  }
}
