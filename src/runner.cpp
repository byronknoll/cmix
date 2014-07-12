#include <fstream>
#include <ctime>
#include <stdio.h>

#include "preprocess/preprocessor.h"
#include "coder/encoder.h"
#include "coder/decoder.h"

void WriteHeader(unsigned long long length, std::ofstream* os) {
  for (int i = 4; i >= 0; --i) {
    char c = length >> (8*i);
    os->put(c);
  }
}

void ReadHeader(std::ifstream* is, unsigned long long* length) {
  *length = 0;
  for (int i = 0; i < 5; ++i) {
    *length <<= 8;
    *length += (unsigned char)(is->get());
  }
}

void Compress(unsigned long long input_bytes, std::ifstream* is,
    std::ofstream* os, unsigned long long* output_bytes) {
  Encoder e(os);
  unsigned long long percent = 1 + (input_bytes / 100);
  for (unsigned long long pos = 0; pos < input_bytes; ++pos) {
    char c = is->get();
    for (int j = 7; j >= 0; --j) {
      e.Encode((c>>j)&1);
    }
    if (pos % percent == 0) {
      printf("\r%lld%%", pos / percent);
      fflush(stdout);
    }
  }
  e.Flush();
  *output_bytes = os->tellp();
}

void Decompress(unsigned long long output_length, std::ifstream* is,
                std::ofstream* os) {
  Decoder d(is);
  unsigned long long percent = 1 + (output_length / 100);
  for(unsigned long long pos = 0; pos < output_length; ++pos) {
    int byte = 1;
    while (byte < 256) {
      byte += byte + d.Decode();
    }
    os->put(byte);
    if (pos % percent == 0) {
      printf("\r%lld%%", pos / percent);
      fflush(stdout);
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc != 4 || argv[1][0] != '-' ||
      (argv[1][1] != 'c' && argv[1][1] != 'd')) {
    printf("To compress:   cmix -c input output\n");
    printf("To decompress: cmix -d input output\n");
    return -1;
  }

  clock_t start = clock();

  bool enable_preprocess = true;

  std::string temp_path = argv[3];
  if (enable_preprocess) temp_path += ".cmix.temp";

  unsigned long long input_bytes, output_bytes;

  if (argv[1][1]=='c') {
    if (enable_preprocess) {
      FILE* data_in = fopen(argv[2], "rb");
      FILE* temp_out = fopen(temp_path.c_str(), "wb");
      if (!data_in || !temp_out) return -1;

      fseek(data_in, 0L, SEEK_END);
      input_bytes = ftell(data_in);
      fseek(data_in, 0L, SEEK_SET);

      preprocessor::encode(data_in, temp_out, input_bytes);
      fclose(data_in);
      fclose(temp_out);
    } else {
      temp_path = argv[2];
    }

    std::ifstream temp_in(temp_path, std::ios::in | std::ios::binary);
    std::ofstream data_out(argv[3], std::ios::out | std::ios::binary);
    if (!temp_in.is_open() || !data_out.is_open()) return -1;

    temp_in.seekg(0, std::ios::end);
    unsigned long long temp_bytes = temp_in.tellg();
    if (!enable_preprocess) input_bytes = temp_bytes;
    temp_in.seekg(0, std::ios::beg);

    WriteHeader(temp_bytes, &data_out);
    Compress(temp_bytes, &temp_in, &data_out, &output_bytes);
    temp_in.close();
    data_out.close();
  } else {
    std::ifstream data_in(argv[2], std::ios::in | std::ios::binary);
    std::ofstream temp_out(temp_path, std::ios::out | std::ios::binary);
    if (!data_in.is_open() || !temp_out.is_open()) return -1;

    data_in.seekg(0, std::ios::end);
    input_bytes = data_in.tellg();
    data_in.seekg(0, std::ios::beg);

    ReadHeader(&data_in, &output_bytes);
    Decompress(output_bytes, &data_in, &temp_out);
    data_in.close();
    temp_out.close();

    if (enable_preprocess) {
      FILE* temp_in = fopen(temp_path.c_str(), "rb");
      FILE* data_out = fopen(argv[3], "wb");
      if (!temp_in || !data_out) return -1;

      preprocessor::decode(temp_in, data_out);
      fseek(data_out, 0L, SEEK_END);
      output_bytes = ftell(data_out);
      fclose(temp_in);
      fclose(data_out);
    }
  }
  if (enable_preprocess) remove(temp_path.c_str());

  printf("\r%lld bytes -> %lld bytes in %1.2f s.\n",
      input_bytes, output_bytes,
      ((double)clock() - start) / CLOCKS_PER_SEC);
  return 0;
}
