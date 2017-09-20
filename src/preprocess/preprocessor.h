#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include "../predictor.h"

#include <stdio.h>

typedef enum {DEFAULT, HDR, JPEG, EXE, TEXT, IMAGE1, IMAGE4, IMAGE8, IMAGE8GRAY, IMAGE24, IMAGE32} Filetype;

inline bool hasInfo(Filetype ft) { return ft==TEXT || ft==IMAGE1 || ft==IMAGE4 || ft==IMAGE8 || ft==IMAGE8GRAY || ft==IMAGE24 || ft==IMAGE32; }

namespace preprocessor {

void pretrain(Predictor* p, const std::vector<bool>& vocab, FILE* dictionary);

void encode(FILE* in, FILE* out, int n, std::string temp_path,
    FILE* dictionary);

void no_preprocess(FILE* in, FILE* out, int n);

void decode(FILE* in, FILE* out, std::string temp_path, FILE* dictionary);

}

#endif
