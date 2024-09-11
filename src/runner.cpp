#include <fstream>
#include <ctime>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string.h>

#include "preprocess/preprocessor.h"
#include "coder/encoder.h"
#include "coder/decoder.h"
#include "predictor.h"

namespace {
  const int kMinVocabFileSize = 10000;
}

char* dictionary_path = NULL;

int Help() {
  printf("cmix version 21\n");
  printf("Compress:\n");
  printf("    with dictionary:    cmix -c [dictionary] [input] [output]\n");
  printf("    without dictionary: cmix -c [input] [output]\n");
  printf("    force text-mode:    cmix -t [dictionary] [input] [output]\n");
  printf("    no preprocessing:   cmix -n [input] [output]\n");
  printf("    only preprocessing: cmix -s [dictionary] [input] [output]\n");
  printf("                        cmix -s [input] [output]\n");
  printf("Decompress:\n");
  printf("    with dictionary:    cmix -d [dictionary] [input] [output]\n");
  printf("    without dictionary: cmix -d [input] [output]\n");
  return -1;
}

void WriteHeader(unsigned long long length, const std::vector<bool>& vocab,
    bool dictionary_used, std::ofstream* os) {
  for (int i = 4; i >= 0; --i) {
    char c = length >> (8*i);
    if (i == 4) {
      c &= 0x7F;
      if (dictionary_used) c |= 0x80;
    }
    os->put(c);
  }
  if (length < kMinVocabFileSize) return;
  for (int i = 0; i < 32; ++i) {
    unsigned char c = 0;
    for (int j = 0; j < 8; ++j) {
      if (vocab[i * 8 + j]) c += 1<<j;
    }
    os->put(c);
  }
}

void WriteStorageHeader(FILE* out, bool dictionary_used) {
  for (int i = 4; i >= 0; --i) {
    char c = 0;
    if (i == 4 && dictionary_used) c = 0x80;
    putc(c, out);
  }
}

void ReadHeader(std::ifstream* is, unsigned long long* length,
    bool* dictionary_used, std::vector<bool>* vocab) {
  *length = 0;
  for (int i = 0; i <= 4; ++i) {
    *length <<= 8;
    unsigned char c = is->get();
    if (i == 0) {
      if (c&0x80) *dictionary_used = true;
      else *dictionary_used = false;
      c &= 0x7F;
    }
    *length += c;
  }
  if (*length == 0) return;
  if (*length < kMinVocabFileSize) {
    std::fill(vocab->begin(), vocab->end(), true);
    return;
  }
  for (int i = 0; i < 32; ++i) {
    unsigned char c = is->get();
    for (int j = 0; j < 8; ++j) {
      if (c & (1<<j)) (*vocab)[i * 8 + j] = true;
    }
  }
}

void ExtractVocab(unsigned long long input_bytes, std::ifstream* is,
    std::vector<bool>* vocab) {
  for (unsigned long long pos = 0; pos < input_bytes; ++pos) {
    unsigned char c = is->get();
    (*vocab)[c] = true;
  }
}

void ClearOutput() {
  fprintf(stderr, "\r                     \r");
  fflush(stderr);
}

void Compress(unsigned long long input_bytes, std::ifstream* is,
    std::ofstream* os, unsigned long long* output_bytes, Predictor* p) {
  Encoder e(os, p);
  unsigned long long percent = 1 + (input_bytes / 10000);
  ClearOutput();
  for (unsigned long long pos = 0; pos < input_bytes; ++pos) {
    char c = is->get();
    for (int j = 7; j >= 0; --j) {
      e.Encode((c>>j)&1);
    }
    if (pos % percent == 0) {
      double frac = 100.0 * pos / input_bytes;
      fprintf(stderr, "\rprogress: %.2f%%", frac);
      fflush(stderr);
    }
  }
  e.Flush();
  *output_bytes = os->tellp();
}

void Decompress(unsigned long long output_length, std::ifstream* is,
                std::ofstream* os, Predictor* p) {
  Decoder d(is, p);
  unsigned long long percent = 1 + (output_length / 10000);
  ClearOutput();
  for(unsigned long long pos = 0; pos < output_length; ++pos) {
    int byte = 1;
    while (byte < 256) {
      byte += byte + d.Decode();
    }
    os->put(byte);
    if (pos % percent == 0) {
      double frac = 100.0 * pos / output_length;
      fprintf(stderr, "\rprogress: %.2f%%", frac);
      fflush(stderr);
    }
  }
}

bool Store(const std::string& input_path, const std::string& temp_path,
    const std::string& output_path, FILE* dictionary,
    unsigned long long* input_bytes, unsigned long long* output_bytes) {
  FILE* data_in = fopen(input_path.c_str(), "rb");
  if (!data_in) return false;
  FILE* data_out = fopen(output_path.c_str(), "wb");
  if (!data_out) return false;
  fseek(data_in, 0L, SEEK_END);
  *input_bytes = ftell(data_in);
  fseek(data_in, 0L, SEEK_SET);
  WriteStorageHeader(data_out, dictionary != NULL);
  fprintf(stderr, "\rpreprocessing...");
  fflush(stderr);
  preprocessor::Encode(data_in, data_out, false, *input_bytes, temp_path,
      dictionary);
  fseek(data_out, 0L, SEEK_END);
  *output_bytes = ftell(data_out);
  fclose(data_in);
  fclose(data_out);
  return true;
}

bool RunCompression(bool enable_preprocess, bool text_mode,
    const std::string& input_path, const std::string& temp_path,
    const std::string& output_path, FILE* dictionary,
    unsigned long long* input_bytes, unsigned long long* output_bytes) {
  FILE* data_in = fopen(input_path.c_str(), "rb");
  if (!data_in) return false;
  FILE* temp_out = fopen(temp_path.c_str(), "wb");
  if (!temp_out) return false;

  fseek(data_in, 0L, SEEK_END);
  *input_bytes = ftell(data_in);
  fseek(data_in, 0L, SEEK_SET);

  if (enable_preprocess) {
    fprintf(stderr, "\rpreprocessing...");
    fflush(stderr);
    preprocessor::Encode(data_in, temp_out, text_mode, *input_bytes, temp_path,
        dictionary);
  } else {
    preprocessor::NoPreprocess(data_in, temp_out, *input_bytes);
  }
  fclose(data_in);
  fclose(temp_out);

  std::ifstream temp_in(temp_path, std::ios::in | std::ios::binary);
  if (!temp_in.is_open()) return false;

  std::ofstream data_out(output_path, std::ios::out | std::ios::binary);
  if (!data_out.is_open()) return false;

  temp_in.seekg(0, std::ios::end);
  unsigned long long temp_bytes = temp_in.tellg();
  temp_in.seekg(0, std::ios::beg);

  std::vector<bool> vocab(256, false);
  if (temp_bytes < kMinVocabFileSize) {
    std::fill(vocab.begin(), vocab.end(), true);
  } else {
    ExtractVocab(temp_bytes, &temp_in, &vocab);
    temp_in.seekg(0, std::ios::beg);
  }

  WriteHeader(temp_bytes, vocab, dictionary != NULL, &data_out);
  Predictor p(vocab);
  if (enable_preprocess) preprocessor::Pretrain(&p, dictionary);
  Compress(temp_bytes, &temp_in, &data_out, output_bytes, &p);
  temp_in.close();
  data_out.close();
  remove(temp_path.c_str());
  return true;
}

bool RunDecompression(const std::string& input_path,
    const std::string& temp_path, const std::string& output_path,
    FILE* dictionary, unsigned long long* input_bytes,
    unsigned long long* output_bytes) {
  std::ifstream data_in(input_path, std::ios::in | std::ios::binary);
  if (!data_in.is_open()) return false;

  data_in.seekg(0, std::ios::end);
  *input_bytes = data_in.tellg();
  data_in.seekg(0, std::ios::beg);
  std::vector<bool> vocab(256, false);
  bool dictionary_used;
  ReadHeader(&data_in, output_bytes, &dictionary_used, &vocab);
  if (!dictionary_used && dictionary != NULL) return false;
  if (dictionary_used && dictionary == NULL) return false;

  if (*output_bytes == 0) {  // undo store
    data_in.close();
    FILE* in = fopen(input_path.c_str(), "rb");
    if (!in) return false;
    FILE* data_out = fopen(output_path.c_str(), "wb");
    if (!data_out) return false;
    fseek(in, 5L, SEEK_SET);
    fprintf(stderr, "\rdecoding...");
    fflush(stderr);
    preprocessor::Decode(in, data_out, dictionary);
    fseek(data_out, 0L, SEEK_END);
    *output_bytes = ftell(data_out);
    fclose(in);
    fclose(data_out);
    return true;
  }
  Predictor p(vocab);
  if (dictionary_used) preprocessor::Pretrain(&p, dictionary);

  std::ofstream temp_out(temp_path, std::ios::out | std::ios::binary);
  if (!temp_out.is_open()) return false;

  Decompress(*output_bytes, &data_in, &temp_out, &p);
  data_in.close();
  temp_out.close();

  FILE* temp_in = fopen(temp_path.c_str(), "rb");
  if (!temp_in) return false;
  FILE* data_out = fopen(output_path.c_str(), "wb");
  if (!data_out) return false;

  preprocessor::Decode(temp_in, data_out, dictionary);
  fseek(data_out, 0L, SEEK_END);
  *output_bytes = ftell(data_out);
  fclose(temp_in);
  fclose(data_out);
  remove(temp_path.c_str());
  return true;
}

int main(int argc, char* argv[]) {
  if (argc < 4 || argc > 5 || strlen(argv[1]) != 2 || argv[1][0] != '-' ||
      (argv[1][1] != 'c' && argv[1][1] != 'd' && argv[1][1] != 's' &&
      argv[1][1] != 'n' && argv[1][1] != 't')) {
    return Help();
  }

  clock_t start = clock();

  bool enable_preprocess = true;
  bool text_mode = false;
  if (argv[1][1] == 'n') enable_preprocess = false;
  std::string input_path = argv[2];
  std::string output_path = argv[3];
  FILE* dictionary = NULL;
  if (argc == 5) {
    if (argv[1][1] == 'n') return Help();
    if (argv[1][1] == 't') text_mode = true;
    dictionary = fopen(argv[2], "rb");
    if (!dictionary) return Help();
    dictionary_path = argv[2];
    input_path = argv[3];
    output_path = argv[4];
  } else {
    if (argv[1][1] == 't') return Help();
  }

  std::string temp_path = output_path + ".cmix.temp";

  unsigned long long input_bytes = 0, output_bytes = 0;

  if (argv[1][1] == 's') {
    if (!Store(input_path, temp_path, output_path, dictionary, &input_bytes,
        &output_bytes)) {
      return Help();
    }
  } else if (argv[1][1] == 'c' || argv[1][1] == 'n' || argv[1][1] == 't') {
    if (!RunCompression(enable_preprocess, text_mode, input_path, temp_path,
        output_path, dictionary, &input_bytes, &output_bytes)) {
      return Help();
    }
  } else {
    if (!RunDecompression(input_path, temp_path, output_path, dictionary,
        &input_bytes, &output_bytes)) {
      return Help();
    }
  }

  printf("\r%lld bytes -> %lld bytes in %1.2f s.\n",
      input_bytes, output_bytes,
      ((double)clock() - start) / CLOCKS_PER_SEC);

  if (argv[1][1] == 'c') {
    double cross_entropy = output_bytes;
    cross_entropy /= input_bytes;
    cross_entropy *= 8;
    printf("cross entropy: %.3f\n", cross_entropy);
  }

  return 0;
}
