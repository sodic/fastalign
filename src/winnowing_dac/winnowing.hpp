//
// Created by dario on 02.11.17..
//

#ifndef PROJEKT_MINIMIZER_H
#define PROJEKT_MINIMIZER_H

#include <vector>
#include <string>
#include <unordered_set>
#include <map>
#include <unordered_map>


namespace winnowing {
    typedef uint64_t minhash_t;

    struct minimizer {
        minhash_t hash;
        int index;
        int strand;

        minimizer(minhash_t _h, int _position, bool _rc)
                : hash(_h), index(_position), strand(_rc) {}

        inline bool operator<(const minimizer &other) const {
            return hash < other.hash;
        }

    };

    // vrati minimizere redom kako se oni nalaze u stringu
    void compute_minimizers(const char *target,
                            uint32_t target_length,
                            int w,
                            int k,
                            std::vector<minimizer> &minimizers);

    void index_sequence(const char *sequence,
                        uint32_t seq_length,
                        uint32_t w,
                        uint32_t k,
                        std::vector<minimizer> &minimizers,
                        std::unordered_map<minhash_t, std::vector<uint32_t >> &lookup_table);

} // namespace winnowing

#endif //PROJEKT_MINIMIZER_H
