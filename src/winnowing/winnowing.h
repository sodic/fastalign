//
// Created by filip on 13.04.18..
//

#ifndef ZAVRSNI_WINNOWING_H
#define ZAVRSNI_WINNOWING_H

#include <vector>
#include <unordered_map>

namespace winnowing {

    const uint32_t DEFAULT_W = 5;
    const uint32_t DEFAULT_K = 15;

    typedef uint64_t minhash_t;

    typedef struct {
        minhash_t hash;
        uint32_t index;
    } minimizer;

    void compute_minimizers(const char *seq,
                            uint32_t seq_l,
                            int32_t w,
                            uint32_t k,
                            std::vector<minimizer> &minimizers);


    void index_sequence(const char *sequence,
                        uint32_t sequence_l,
                        uint32_t w,
                        uint32_t k,
                        std::vector<winnowing::minimizer> &minimizers,
                        std::unordered_map<minhash_t, std::vector<std::int32_t >> &lookup_table);

}

#endif //ZAVRSNI_WINNOWING_H