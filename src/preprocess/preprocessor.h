#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include <stdio.h>

namespace preprocessor {

typedef enum {DEFAULT, HDR, JPEG, EXE, TEXT, IMAGE1, IMAGE4, IMAGE8, IMAGE8GRAY,
    IMAGE24, IMAGE32} Filetype;

inline bool HasInfo(Filetype ft) { return ft==TEXT || ft==IMAGE1 || ft==IMAGE4
    || ft==IMAGE8 || ft==IMAGE8GRAY || ft==IMAGE24 || ft==IMAGE32; }

void Encode(FILE* in, FILE* out, int n, std::string temp_path,
    FILE* dictionary);

void NoPreprocess(FILE* in, FILE* out, int n);

void Decode(FILE* in, FILE* out, std::string temp_path, FILE* dictionary);

}

#endif
