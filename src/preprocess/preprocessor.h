#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include "../predictor.h"

#include <stdio.h>

namespace preprocessor {

typedef enum {DEFAULT, JPEG, EXE, TEXT, BMP} Filetype;

void pretrain(Predictor* p, const std::vector<bool>& vocab, FILE* dictionary);

void encode(FILE* in, FILE* out, int n, std::string temp_path,
    FILE* dictionary);

void no_preprocess(FILE* in, FILE* out, int n);

void decode(FILE* in, FILE* out, std::string temp_path, FILE* dictionary);

}

#endif
