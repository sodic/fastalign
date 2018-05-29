//
// Created by filip on 28.05.18..
//

#ifndef ZAVRSNI_TYPES_H
#define ZAVRSNI_TYPES_H

#endif //ZAVRSNI_TYPES_H


typedef uint64_t minhash_t;

typedef struct {
    minhash_t hash;
    uint32_t index;
    int32_t strand;
} minimizer;