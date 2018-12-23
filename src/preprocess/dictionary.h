#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <stdio.h>
#include <unordered_map>
#include <string>
#include <array>
#include <list>

namespace preprocessor {

class Dictionary {
 public:
  Dictionary(FILE* dictionary, bool encode, bool decode);
  void Encode(FILE* input, int len, FILE* output);
  unsigned char Decode(FILE* input);

 private:
  void EncodeWord(const std::string& word, int num_upper, bool next_lower,
      FILE* output);
  bool EncodeSubstring(const std::string& word, FILE* output);
  void AddToBuffer(FILE* input);

  std::unordered_map<std::string, unsigned int> byte_map_;
  std::unordered_map<unsigned int, std::string> reverse_map_;
  std::list<unsigned char> output_buffer_;
  bool decode_upper_ = false, decode_capital_ = false;
  unsigned int longest_word_ = 0;
};

}

#endif
