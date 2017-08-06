#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include "../predictor.h"

#include <stdio.h>

namespace preprocessor {

void pretrain(Predictor* p, const std::vector<bool>& vocab, FILE* dictionary);

void encode(FILE* in, FILE* out, int n, std::string temp_path,
    FILE* dictionary);

void decode(FILE* in, FILE* out, std::string temp_path, FILE* dictionary);

}

#endif
