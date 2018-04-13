//
// Created by filip on 13.04.18..
//

#ifndef ZAVRSNI_MAPPER_H
#define ZAVRSNI_MAPPER_H

#include <FASTARead/FASTARead.h>
#include <winnowing/winnowing.h>

namespace mapper {

    typedef struct{
        int32_t start;
        int32_t end;
    } region;

    void find_candidates(const char *readA,
                         uint32_t lengthA,
                         std::unordered_map<winnowing::minhash_t, std::vector<int32_t>> lookup_map,
                         uint32_t s,
                         double tau);
}

#endif //ZAVRSNI_MAPPER_H
