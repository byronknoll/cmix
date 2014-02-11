#include <fstream>
#include <ctime>

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
  Encoder e(os, input_bytes);
  int percent = 0;
  for (unsigned long long pos = 0; pos < input_bytes; ++pos) {
    char c = is->get();
    for (int j = 7; j >= 0; --j) {
      e.Encode((c>>j)&1);
    }
    if ((100.0 * pos) / input_bytes > percent) {
      printf("%d%%\r", percent);
      fflush(stdout);
      ++percent;
    }
  }
  e.Flush();
  *output_bytes = os->tellp();
}

void Decompress(unsigned long long output_length, std::ifstream* is,
                std::ofstream* os) {
  Decoder d(is, output_length);
  int percent = 0;
  for(unsigned long long pos = 0; pos < output_length; ++pos) {
    int byte = 1;
    while (byte < 256) {
      byte += byte + d.Decode();
    }
    os->put(byte);
    if ((100.0 * pos) / output_length > percent) {
      printf("%d%%\r", percent);
      fflush(stdout);
      ++percent;
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

  std::ifstream is(argv[2], std::ios::in | std::ios::binary);
  std::ofstream os(argv[3], std::ios::out | std::ios::binary);
  if (!is.is_open() || !os.is_open()) {
    printf("Error opening file.\n");
    abort();
  }

  is.seekg(0, std::ios::end);
  unsigned long long input_bytes = is.tellg();
  is.seekg(0, std::ios::beg);
  unsigned long long output_bytes;

  if (argv[1][1]=='c') {
    WriteHeader(input_bytes, &os);
    Compress(input_bytes, &is, &os, &output_bytes);
  } else {
    ReadHeader(&is, &output_bytes);
    Decompress(output_bytes, &is, &os);
  }

  is.close();
  os.close();

  printf("%lld bytes -> %lld bytes in %1.2f s.\n",
      input_bytes, output_bytes,
      ((double)clock()-start)/CLOCKS_PER_SEC);
  return 0;
}
