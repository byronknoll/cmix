#include "decoder.h"

Decoder::Decoder(std::ifstream* is, Predictor* p) : is_(is), x1_(0),
    x2_(0xffffffff), x_(0), p_(p) {
  for (int i = 0; i < 4; ++i) {
    x_ = (x_ << 8) + (ReadByte() & 0xff);
  }
}

int Decoder::ReadByte() {
  int byte = (unsigned char)(is_->get());
  if (!is_->good()) return 0;
  return byte;
}

unsigned int Decoder::Discretize(float p) {
  return 1 + 65534 * p;
}

int Decoder::Decode() {
  const unsigned int p = Discretize(p_->Predict());
  const unsigned int xmid = x1_ + ((x2_ - x1_) >> 16) * p +
      (((x2_ - x1_) & 0xffff) * p >> 16);
  int bit = 0;
  if (x_ <= xmid) {
    bit = 1;
    x2_ = xmid;
  } else {
    x1_ = xmid + 1;
  }
  p_->Perceive(bit);

  while (((x1_^x2_) & 0xff000000) == 0) {
    x1_ <<= 8;
    x2_ = (x2_ << 8) + 255;
    x_ = (x_ << 8) + ReadByte();
  }
  return bit;
}
