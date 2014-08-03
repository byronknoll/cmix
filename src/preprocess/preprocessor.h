#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include "../predictor.h"

#include <stdio.h>

namespace preprocessor {

void pretrain(Predictor* p, FILE* dictionary);

void encode(FILE* in, FILE* out, int n, FILE* dictionary);

void decode(FILE* in, FILE* out, FILE* dictionary);

}  // namespace preprocessor

#endif
