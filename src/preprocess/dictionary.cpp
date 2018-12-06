#include "dictionary.h"

#include <stdio.h>
#include <unordered_set>
#include <string>

namespace preprocessor {

namespace {

const unsigned char kCapitalized = 0x40;
const unsigned char kUppercase = 0x07;
const unsigned char kEndUpper = 0x06;
const unsigned char kEscape = 0x0C;

void EncodeByte(unsigned char c, FILE* output) {
  if (c == kEndUpper || c == kEscape || c == kUppercase ||
      c == kCapitalized || c >= 0x80) {
    putc(kEscape, output);
  }
  putc(c, output);
}

void EncodeBytes(unsigned int bytes, FILE* output) {
  putc((unsigned char)bytes&0xFF, output);
  if (bytes & 0xFF00) {
    putc((unsigned char)((bytes&0xFF00)>>8), output);
  } else {
    return;
  }
  if (bytes & 0xFF0000) {
    putc((unsigned char)((bytes&0xFF0000)>>16), output);
  }
}

}

Dictionary::Dictionary(FILE* dictionary, bool encode, bool decode) {
  fseek(dictionary, 0L, SEEK_END);
  unsigned long long len = ftell(dictionary);
  fseek(dictionary, 0L, SEEK_SET);
  std::string line;
  int line_count = 0;
  const int kBoundary1 = 80, kBoundary2 = kBoundary1 + 3840,
      kBoundary3 = kBoundary2 + 40960;
  for (unsigned pos = 0; pos < len; ++pos) {
    unsigned char c = getc(dictionary);
    if (c >= 'a' && c <= 'z') line += c;
    else if (!line.empty()) {
      unsigned int bytes;
      if (line_count < kBoundary1) {
        bytes = 0x80 + line_count;
      } else if (line_count < kBoundary2) {
        bytes = 0xD0 + ((line_count-kBoundary1) / 80);
        bytes += (0x80 + ((line_count-kBoundary1) % 80)) << 8;
      } else if (line_count < kBoundary3) {
        bytes = 0xF0 + (((line_count-kBoundary2) / 80) / 32);
        bytes += (0xD0 + (((line_count-kBoundary2) / 80) % 32)) << 8;
        bytes += (0x80 + ((line_count-kBoundary2) % 80)) << 16;
      }
      if (encode) byte_map_[line] = bytes;
      if (decode) reverse_map_[bytes] = line;
      ++line_count;
      line.clear();
    }
  }
}

void Dictionary::Encode(FILE* input, int len, FILE* output) {
  std::string word;
  int num_upper = 0, num_lower = 0;
  for (int pos = 0; pos < len; ++pos) {
    unsigned char c = getc(input);
    bool advance = false;
    if (word.size() > 50) {
      advance = true;
    } else if (c >= 'a' && c <= 'z') {
      if (num_upper > 1) {
        advance = true;
      } else {
        ++num_lower;
        word += c;
      }
    } else if (c >= 'A' && c <= 'Z') {
      if (num_lower > 0) {
        advance = true;
      } else {
        ++num_upper;
        word += (c - 'A' + 'a');
      }
    } else {
      advance = true;
    }
    if (pos == len - 1 && !advance) {
      EncodeWord(word, num_upper, false, output);
    }
    if (advance) {
      if (word.empty()) {
        EncodeByte(c, output);
      } else {
        bool next_lower = (c >= 'a' && c <= 'z');
        EncodeWord(word, num_upper, next_lower, output);
        num_lower = 0;
        num_upper = 0;
        word.clear();
        if (next_lower) {
          ++num_lower;
          word += c;
        } else if (c >= 'A' && c <= 'Z') {
          ++num_upper;
          word += (c - 'A' + 'a');
        } else {
          EncodeByte(c, output);
        }
        if (pos == len - 1 && !word.empty()) {
          EncodeWord(word, num_upper, false, output);
        }
      }
    }
  }
}

void Dictionary::EncodeWord(const std::string& word, int num_upper,
    bool next_lower, FILE* output) {
  if (num_upper > 1) putc(kUppercase, output);
  else if (num_upper == 1) putc(kCapitalized, output);
  if (byte_map_.find(word) != byte_map_.end()) {
    EncodeBytes(byte_map_[word], output);
  } else if (!EncodeSuffix(word, output) && !EncodePrefix(word, output)) {
    for (unsigned int i = 0; i < word.size(); ++i) {
      putc((unsigned char)word[i], output);
    }
  }
  if (num_upper > 1 && next_lower) {
    putc(kEndUpper, output);
  }
}

bool Dictionary::EncodePrefix(const std::string& word, FILE* output) {
  if (word.size() <= 7) return false;
  std::string prefix = word.substr(0, word.size() - 1);
  while (prefix.size() >= 7) {
    if (byte_map_.find(prefix) != byte_map_.end()) {
      EncodeBytes(byte_map_[prefix], output);
      for (unsigned int i = prefix.size(); i < word.size(); ++i) {
        putc((unsigned char)word[i], output);
      }
      return true;
    }
    prefix.erase(prefix.size() - 1, 1);
  }
  return false;
}

bool Dictionary::EncodeSuffix(const std::string& word, FILE* output) {
  if (word.size() <= 7) return false;
  std::string suffix = word.substr(1, word.size() - 1);
  while (suffix.size() >= 7) {
    if (byte_map_.find(suffix) != byte_map_.end()) {
      for (unsigned int i = 0; i < word.size() - suffix.size(); ++i) {
        putc((unsigned char)word[i], output);
      }
      EncodeBytes(byte_map_[suffix], output);
      return true;
    }
    suffix.erase(0, 1);
  }
  return false;
}

unsigned char Dictionary::Decode(FILE* input) {
  while (output_buffer_.empty()) {
    AddToBuffer(input);
  }
  unsigned char next = output_buffer_.front();
  output_buffer_.pop_front();
  return next;
}

void Dictionary::AddToBuffer(FILE* input) {
  unsigned char c = getc(input);
  if (!(c == kEndUpper || c == kEscape || c == kUppercase ||
          c == kCapitalized || c >= 0x80)) {
    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))) {
      decode_upper_ = false;
    }
    if (decode_capital_ || decode_upper_) {
      c = (c - 'a') + 'A';
    }
    if (decode_capital_) decode_capital_ = false;
    output_buffer_.push_back(c);
    return;
  }
  if (c == kEscape) {
    decode_upper_ = false;
    output_buffer_.push_back(getc(input));
    return;
  }
  if (c == kUppercase) {
    decode_upper_ = true;
    return;
  }
  if (c == kCapitalized) {
    decode_capital_ = true;
    return;
  }
  if (c == kEndUpper) {
    decode_upper_ = false;
    return;
  }
  if (c >= 0x80) {
    unsigned int bytes = c;
    if (c > 0xCF) {
      c = getc(input);
      bytes += c << 8;
      if (c > 0xCF) {
        c = getc(input);
        bytes += c << 16;
      }
    }
    std::string word = reverse_map_[bytes];
    for (unsigned int i = 0; i < word.size(); ++i) {
      if (i == 0 && decode_capital_) {
        word[i] = (word[i] - 'a') + 'A';
        decode_capital_ = false;
      }
      if (decode_upper_) {
        word[i] = (word[i] - 'a') + 'A';
      }
      output_buffer_.push_back(word[i]);
    }
    return;
  }
  output_buffer_.push_back(c);
}

}
