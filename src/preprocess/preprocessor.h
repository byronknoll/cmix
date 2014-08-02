#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include "../predictor.h"

#include <stdio.h>

namespace preprocessor {

void pretrain(Predictor* p);

void encode(FILE* in, FILE* out, int n);

void decode(FILE* in, FILE* out);

}  // namespace preprocessor

#endif
