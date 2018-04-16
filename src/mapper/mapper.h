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

    typedef struct{
        uint32_t position;
        winnowing::minhash_t jaccard;
    } estimate;

    void find_candidates(std::vector<winnowing::minimizer>& A_minimizers,
                         uint32_t lengthA,
                         std::unordered_map<winnowing::minhash_t, std::vector<int32_t>> &lookup_table,
                         uint32_t s,
                         double tau,
                         std::vector<region> &candidates);

    void compute_estimates(std::vector<winnowing::minimizer> &ref_minimizers,
                           std::vector<winnowing::minimizer> &query_minimizers,
                           uint32_t lengthA,
                           std::vector<region> &candidates,
                           uint32_t s,
                           double tau,
                           std::vector<estimate> estimates);
}

#endif //ZAVRSNI_MAPPER_H
