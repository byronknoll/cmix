// This preprocessor is adapted from paq8l and paq8hp12any.

#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string.h>

#define preprocFlag 1220
#define OPTION_UTF8							1
#define OPTION_USE_NGRAMS					2
#define OPTION_CAPITAL_CONVERSION			4
#define OPTION_WORD_SURROROUNDING_MODELING	8
#define OPTION_SPACE_AFTER_EOL				16
#define OPTION_EOL_CODING					32
#define OPTION_NORMAL_TEXT_FILTER			64
#define OPTION_USE_DICTIONARY				128
#define OPTION_RECORD_INTERLEAVING			256
#define OPTION_DNA_QUARTER_BYTE				512
#define OPTION_TRY_SHORTER_WORD				1024
#define OPTION_TO_LOWER_AFTER_PUNCTUATION	2048
#define OPTION_SPACELESS_WORDS				4096
#define OPTION_ADD_SYMBOLS_0_5				8192
#define OPTION_ADD_SYMBOLS_14_31			16384
#define OPTION_ADD_SYMBOLS_A_Z				32768
#define OPTION_ADD_SYMBOLS_MISC				65536
#define OPTION_SPACE_AFTER_CC_FLAG			131072
#define IF_OPTION(option) ((preprocFlag & option)!=0)

#include "textfilter.cpp"

namespace preprocessor {

typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

typedef enum {DEFAULT, JPEG, EXE, TEXT} Filetype;

inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}

/////////////////////////// Filters /////////////////////////////////
//
// Before compression, data is encoded in blocks with the following format:
//
//   <type> <size> <encoded-data>
//
// Type is 1 byte (type Filetype): DEFAULT=0, JPEG, EXE
// Size is 4 bytes in big-endian format.
// Encoded-data decodes to <size> bytes.  The encoded size might be
// different.  Encoded data is designed to be more compressible.
//
//   void encode(FILE* in, FILE* out, int n);
//
// Reads n bytes of in (open in "rb" mode) and encodes one or
// more blocks to temporary file out (open in "wb+" mode).
// The file pointer of in is advanced n bytes.  The file pointer of
// out is positioned after the last byte written.
//
//   en.setFile(FILE* out);
//   int decode(Encoder& en);
//
// Decodes and returns one byte.  Input is from en.decompress(), which
// reads from out if in COMPRESS mode.  During compression, n calls
// to decode() must exactly match n bytes of in, or else it is compressed
// as type 0 without encoding.
//
//   Filetype detect(FILE* in, int n, Filetype type);
//
// Reads n bytes of in, and detects when the type changes to
// something else.  If it does, then the file pointer is repositioned
// to the start of the change and the new type is returned.  If the type
// does not change, then it repositions the file pointer n bytes ahead
// and returns the old type.
//
// For each type X there are the following 2 functions:
//
//   void encode_X(FILE* in, FILE* out, int n, ...);
//
// encodes n bytes from in to out.
//
//   int decode_X(Encoder& en);
//
// decodes one byte from en and returns it.  decode() and decode_X()
// maintain state information using static variables.

bool IsAscii(int byte) {
  if (byte >= 9 && byte <= 13) return true;
  if (byte >= 32 && byte <= 126) return true;
  return false;
}

// Detect data type
Filetype detect(FILE* in, int n, Filetype type) {
  U32 buf1=0, buf0=0;  // last 8 bytes
  long start=ftell(in);

  // For EXE detection
  std::vector<int> abspos(256, 0),  // CALL/JMP abs. addr. low byte -> last offset
    relpos(256, 0);    // CALL/JMP relative addr. low byte -> last offset
  int e8e9count=0;  // number of consecutive CALL/JMPs
  int e8e9pos=0;    // offset of first CALL or JMP instruction
  int e8e9last=0;   // offset of most recent CALL or JMP

  // For JPEG detection
  int soi=0, sof=0, sos=0;  // position where found

  // For TEXT detection
  int ascii_start = -1;
  int ascii_run = 0;

  for (int i=0; i<n; ++i) {
    int c=getc(in);
    if (c==EOF) return (Filetype)(-1);
    buf1=buf1<<8|buf0>>24;
    buf0=buf0<<8|c;

    // Detect JPEG by code SOI APPx (FF D8 FF Ex) followed by
    // SOF0 (FF C0 xx xx 08) and SOS (FF DA) within a reasonable distance.
    // Detect end by any code other than RST0-RST7 (FF D9-D7) or
    // a byte stuff (FF 00).

    if (i>=3 && (buf0&0xfffffff0)==0xffd8ffe0) soi=i;
    if (soi && i-soi<0x10000 && (buf1&0xff)==0xff
        && (buf0&0xff0000ff)==0xc0000008)
      sof=i;
    if (soi && sof && sof>soi && i-soi<0x10000 && i-sof<0x1000
        && (buf0&0xffff)==0xffda) {
      sos=i;
      if (type!=JPEG) return fseek(in, start+soi-3, SEEK_SET), JPEG;
    }
    if (type==JPEG && sos && i>sos && (buf0&0xff00)==0xff00
        && (buf0&0xff)!=0 && (buf0&0xf8)!=0xd0)
      return DEFAULT;

    // Detect EXE if the low order byte (little-endian) XX is more
    // recently seen (and within 4K) if a relative to absolute address
    // conversion is done in the context CALL/JMP (E8/E9) XX xx xx 00/FF
    // 4 times in a row.  Detect end of EXE at the last
    // place this happens when it does not happen for 64KB.

    if ((buf1&0xfe)==0xe8 && ((buf0+1)&0xfe)==0) {
      int r=buf0>>24;  // relative address low 8 bits
      int a=((buf0>>24)+i)&0xff;  // absolute address low 8 bits
      int rdist=i-relpos[r];
      int adist=i-abspos[a];
      if (adist<rdist && adist<0x1000 && abspos[a]>5) {
        e8e9last=i;
        ++e8e9count;
        if (e8e9pos==0 || e8e9pos>abspos[a]) e8e9pos=abspos[a];
      }
      else e8e9count=0;
      if (type!=EXE && e8e9count>=4 && e8e9pos>5)
        return fseek(in, start+e8e9pos-5, SEEK_SET), EXE;
      abspos[a]=i;
      relpos[r]=i;
    }
    if (type==EXE && i-e8e9last>0x1000)
      return fseek(in, start+e8e9last, SEEK_SET), DEFAULT;

    // Detect TEXT
    if (type == DEFAULT) {
      if (IsAscii(c)) {
        if (ascii_start == -1) {
          ascii_start = i;
          ascii_run = 0;
        }
        ++ascii_run;
        if (ascii_run > 500) {
          return fseek(in, start + ascii_start, SEEK_SET), TEXT;
        }
      } else {
        ascii_start = -1;
      }
    } else if (type == TEXT) {
      if (IsAscii(c)) {
        ascii_run -= 2;
        if (ascii_run < 0) ascii_run = 0;
      } else {
        ascii_run += 3;
        if (ascii_run > 300) {
          return fseek(in, ftell(in) - 100, SEEK_SET), DEFAULT;
        }
      }
    }
  }
  return type;
}

// Default encoding as self
void encode_default(FILE* in, FILE* out, int len) {
  while (len--) putc(getc(in), out);
}

int decode_default(FILE* in) {
  return getc(in);
}

// JPEG encode as self.  The purpose is to shield jpegs from exe transform.
void encode_jpeg(FILE* in, FILE* out, int len) {
  while (len--) putc(getc(in), out);
}

int decode_jpeg(FILE* in) {
  return getc(in);
}

// EXE transform: <encoded-size> <begin> <block>...
// Encoded-size is 4 bytes, MSB first.
// begin is the offset of the start of the input file, 4 bytes, MSB first.
// Each block applies the e8e9 transform to strings falling entirely
// within the block starting from the end and working backwards.
// The 5 byte pattern is E8/E9 xx xx xx 00/FF (x86 CALL/JMP xxxxxxxx)
// where xxxxxxxx is a relative address LSB first.  The address is
// converted to an absolute address by adding the offset mod 2^25
// (in range +-2^24).

void encode_exe(FILE* in, FILE* out, int len, int begin) {
  const int BLOCK=0x10000;
  std::vector<U8> blk(BLOCK);
  fprintf(out, "%c%c%c%c", len>>24, len>>16, len>>8, len); // size, MSB first
  fprintf(out, "%c%c%c%c", begin>>24, begin>>16, begin>>8, begin); 

  // Transform
  for (int offset=0; offset<len; offset+=BLOCK) {
    int size=min(len-offset, BLOCK);
    int bytesRead=fread(&blk[0], 1, size, in);
    if (bytesRead!=size) abort();
    for (int i=bytesRead-1; i>=4; --i) {
      if ((blk[i-4]==0xe8||blk[i-4]==0xe9) && (blk[i]==0||blk[i]==0xff)) {
        int a=(blk[i-3]|blk[i-2]<<8|blk[i-1]<<16|blk[i]<<24)+offset+begin+i+1;
        a<<=7;
        a>>=7;
        blk[i]=a>>24;
        blk[i-1]=a>>16;
        blk[i-2]=a>>8;
        blk[i-3]=a;
      }
    }
    fwrite(&blk[0], 1, bytesRead, out);
  }
}

int decode_exe(FILE* in) {
  const int BLOCK=0x10000;  // block size
  static int offset=0, q=0;  // decode state: file offset, queue size
  static int size=0;  // where to stop coding
  static int begin=0;  // offset in file
  static U8 c[5];  // queue of last 5 bytes, c[0] at front

  // Read size from first 4 bytes, MSB first
  while (offset==size && q==0) {
    offset=0;
    size=getc(in)<<24;
    size|=getc(in)<<16;
    size|=getc(in)<<8;
    size|=getc(in);
    begin=getc(in)<<24;
    begin|=getc(in)<<16;
    begin|=getc(in)<<8;
    begin|=getc(in);
  }

  // Fill queue
  while (offset<size && q<5) {
    memmove(c+1, c, 4);
    c[0]=getc(in);
    ++q;
    ++offset;
  }

  // E8E9 transform: E8/E9 xx xx xx 00/FF -> subtract location from x
  if (q==5 && (c[4]==0xe8||c[4]==0xe9) && (c[0]==0||c[0]==0xff)
      && (((offset-1)^(offset-5))&-BLOCK)==0) { // not crossing block boundary
    int a=(c[3]|c[2]<<8|c[1]<<16|c[0]<<24)-offset-begin;
    a<<=7;
    a>>=7;
    c[3]=a;
    c[2]=a>>8;
    c[1]=a>>16;
    c[0]=a>>24;
  }

  // return oldest byte in queue
  //assert(q>0 && q<=5);
  return c[--q];
}

void encode_text(FILE* in, FILE* out, int len) {
  FILE* temp_input = tmpfile();
  if (!temp_input) abort();

  for (int i = 0; i < len; ++i) {
    putc(getc(in), temp_input);
  }
  rewind(temp_input);

  FILE* temp_output = tmpfile();
  if (!temp_output) abort();

  WRT wrt;
  wrt.defaultSettings(0, NULL);
  wrt.WRT_start_encoding(temp_input, temp_output, len, false);

  int size = ftell(temp_output);
  if (size > len - 50) {
    for (int i = 0; i < 4; ++i) {
      putc(0, out);
    }
    rewind(temp_input);
    for (int i = 0; i < len; ++i) {
      putc(getc(temp_input), out);
    }
  } else {
    rewind(temp_output);
    for (int i = 0; i < size; ++i) {
      putc(getc(temp_output), out);
    }
  }

  fclose(temp_input);
  fclose(temp_output);
}

FILE* wrt_temp;
WRT* wrt_decoder = NULL;
bool wrt_enabled = true;

void reset_text_decoder(FILE* in) {
  if (wrt_temp) fclose(wrt_temp);
  wrt_temp = tmpfile();
  if (!wrt_temp) abort();

  unsigned int size = 0;
  for (int i = 4; i != 0; --i) {
    int c = getc(in);
    size = size * 256 + c;
    putc(c, wrt_temp);
  }
  if (size == 0) {
    wrt_enabled = false;
    return;
  }
  wrt_enabled = true;
  size -= 8;

  for (unsigned int i = 0; i < size; ++i) putc(getc(in), wrt_temp);
  rewind(wrt_temp);

  if (wrt_decoder != NULL) delete wrt_decoder;
  wrt_decoder = new WRT();
  wrt_decoder->defaultSettings(0, NULL);
  wrt_decoder->WRT_prepare_decoding();
}

int decode_text(FILE* in) {
  if (!wrt_enabled) return getc(in);
  return wrt_decoder->WRT_decode_char(wrt_temp, NULL, 0);
}

// Split n bytes into blocks by type.  For each block, output
// <type> <size> and call encode_X to convert to type X.
void encode(FILE* in, FILE* out, int n) {
  Filetype type=DEFAULT;
  long begin=ftell(in);
  while (n>0) {
    Filetype nextType=detect(in, n, type);
    long end=ftell(in);
    fseek(in, begin, SEEK_SET);
    int len=int(end-begin);
    if (len>0) {
      fprintf(out, "%c%c%c%c%c", type, len>>24, len>>16, len>>8, len);
      //printf("type: %d\tlength: %d\n", type, len);
      switch(type) {
        case JPEG: encode_jpeg(in, out, len); break;
        case EXE:  encode_exe(in, out, len, begin); break;
        case TEXT: encode_text(in, out, len); break;
        default:   encode_default(in, out, len); break;
      }
    }
    n-=len;
    type=nextType;
    begin=end;
  }
}

// Decode <type> <len> <data>...
int decode2(FILE* in) {
  static Filetype type=DEFAULT;
  static int len=0;
  while (len==0) {
    int c = getc(in);
    if (c == EOF) return -1;
    type=(Filetype)c;
    len=getc(in)<<24;
    len|=getc(in)<<16;
    len|=getc(in)<<8;
    len|=getc(in);
    if (len<0) len=1;
    if (type == TEXT) reset_text_decoder(in);
  }
  --len;
  switch (type) {
    case JPEG: return decode_jpeg(in);
    case EXE:  return decode_exe(in);
    case TEXT: return decode_text(in);
    default:   return decode_default(in);
  }
}

void decode(FILE* in, FILE* out) {
  while (true) {
    int result = decode2(in);
    if (result == -1) return;
    putc(result, out);
  }
}

}  // namespace preprocessor
