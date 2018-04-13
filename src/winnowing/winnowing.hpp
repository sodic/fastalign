//
// Created by filip on 13.04.18..
//

#ifndef ZAVRSNI_WINNOWING_H
#define ZAVRSNI_WINNOWING_H

#include <vector>
#include <unordered_map>

namespace winnowing{

    typedef uint64_t minhash_t;

    typedef struct {
        minhash_t hash;
        int32_t index;
    } minimizer;


    void index_sequence(const char* sequence,
                            uint32_t sequence_l,
                            uint32_t w,
                            uint32_t k,
                            std::vector<winnowing::minimizer>& minimizers,
                            std::unordered_map<minhash_t, std::vector<std::int32_t >>& lookup_table);

}

#endif //ZAVRSNI_WINNOWING_H