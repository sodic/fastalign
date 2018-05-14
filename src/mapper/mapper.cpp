//
// Created by filip on 13.04.18..
//

#include <cstdint>
#include <vector>
#include <FASTARead/FASTARead.h>
#include <winnowing/winnowing.hpp>
#include<cmath>
#include <algorithm>
#include <map>
#include <assert.h>
#include "mapper.hpp"


namespace mapper {

    void find_candidates(std::vector<winnowing::minimizer>& query_minimizers,
                         uint32_t query_length,
                         std::unordered_map<winnowing::minhash_t, std::vector<uint32_t>> &lookup_table,
                         uint32_t s,
                         double tau,
                         std::vector<region> &candidates) {
        uint32_t m = (uint32_t) ceil(tau * s);

        std::vector<uint32_t> positions;         // the array L

        for (auto &minimizer : query_minimizers) {
            auto matches = lookup_table.find(minimizer.hash);
            if (matches == lookup_table.end()) {
                continue;
            }
            for (uint32_t pos : matches->second) {
                positions.push_back(pos);
            }
        }

        std::sort(positions.begin(), positions.end());

        int total_positions = (int) positions.size();
        for (int i = 0; i <= total_positions - m; ++i) {
            int j = i + m - 1;

            // check if consecutive hits are close enough
            if (positions[j] - positions[i] < query_length) {

                region candidate = {
                        start: std::max(0u, positions[j] - query_length + 1),
                        end: positions[i]
                };

                //check overlap with the previous candidate
                auto prev_it = candidates.end();
                prev_it--;

                if (!candidates.empty() && prev_it->end >= candidate.start) {
                    // if there is an overlap, extend the previous candidate
                    prev_it->end = std::max(candidate.end, prev_it->end);
                } else {
                    candidates.push_back(candidate);
                }
            }

        }
    }


    void init_L(std::map<winnowing::minhash_t, matchInfo> &L,
                std::vector<winnowing::minimizer> &minimizers){
        for(auto &minimizer : minimizers){
            L[minimizer.hash] = (matchInfo) {
                    mutual: 0,
                    strand: minimizer.strand
            };
        }
    }

    double solve_jaccard(std::map<winnowing::minhash_t, matchInfo> &L, uint32_t s) {
        assert(s > 0);

        uint32_t shared_sketch = 0;
        auto it = L.begin();
        uint32_t collected = 0;
        while (it != L.end() && collected < s) {
            shared_sketch += it->second.mutual;
            it++;
            collected++;
        }
        return 1.0 * shared_sketch / s;
    }


    void remove_from_map(uint32_t start,
                         uint32_t end,
                         std::map<winnowing::minhash_t, matchInfo> &L,
                         std::vector<winnowing::minimizer> &minimizers){
        int i = 0;
        while(minimizers[i].index < start){
            i++;
        }

        while(minimizers[i].index < end){
            L.erase(minimizers[i].hash);
        }
    }


    void insert_into_map(uint32_t start,
                         uint32_t end,
                         std::map<winnowing::minhash_t, matchInfo> &L,
                         std::vector<winnowing::minimizer> &minimizers){
        int i = 0;
        while(minimizers[i].index < start){
            i++;
        }

        while(minimizers[i].index < end){
            if (L.find(minimizers[i].hash) != L.end()) {
                L[minimizers[i].hash].mutual = 1;
                L[minimizers[i].hash].strand *= minimizers[i].strand;
            } else {
                L[minimizers[i].hash] = matchInfo {
                        mutual: 0,
                        strand: minimizers[i].strand
                };
            }
            i++;
        }
    }

    bool consensus_strand(std::map<winnowing::minhash_t, matchInfo> L, uint32_t s) {
        assert(s > 0);

        int32_t sum = 0;
        uint32_t counter = 0;
        auto it = L.begin();
        while (it != L.end() && counter < s) {
            if (it->second.mutual) {
                sum += it->second.strand;
                counter++;
            }
            it++;
        }
        return sum >= 0;
    }



    void compute_estimates(std::vector<winnowing::minimizer> &ref_minimizers,
                           std::vector<winnowing::minimizer> &query_minimizers,
                           uint32_t query_length,
                           std::vector<region> &candidates,
                           uint32_t s,
                           double tau,
                           std::vector<estimate> &estimates){
        for (auto &candidate : candidates) {
            uint32_t  i = candidate.start;
            uint32_t  j = candidate.start + query_length;

            std::map<winnowing::minhash_t, matchInfo> L;
            init_L(L, query_minimizers);                // treba ih i za ostale opet postavit na nulu, ili? mislim da ne
            insert_into_map(i, j, L, ref_minimizers);

            winnowing::minhash_t J = solve_jaccard(L, s);
            if (J >= tau){
                estimates.push_back((estimate) {i, J});
            }

            while(i <= candidate.end){
                L.erase(i);
                remove_from_map(i, i+1, L, ref_minimizers);
                insert_into_map(j, j + 1, L, ref_minimizers);

                J = solve_jaccard(L, s);
                if (J >= tau){
                    estimates.push_back((estimate) {
                            position: i,
                            strand: consensus_strand(L, s),
                            jaccard: J
                    });
                }

                i++;
                j++;
            }
        }
    }


    int binary_search(std::vector<winnowing::minimizer> &minimizers,
                      uint32_t start,
                      uint32_t length,
                      uint32_t element){
        if(!length){
            return -1;
        }

        uint32_t middle = start + length/2;

        if (minimizers[middle].index == element){
            return middle;
        }

        if (minimizers[middle].index < element) {
            return binary_search(minimizers, middle + 1, length/2 - 1 + length%2, element);
        }

        return binary_search(minimizers, start, length/2, element);
    }


}
