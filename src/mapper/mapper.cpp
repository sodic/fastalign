//
// Created by filip on 13.04.18..
//

#include <cstdint>
#include <vector>
#include <FASTARead/FASTARead.h>
#include <winnowing/winnowing.h>
#include<cmath>
#include <algorithm>
#include <map>
#include "mapper.h"


namespace mapper {

    void find_candidates(std::vector<winnowing::minimizer>& A_minimizers,
                         uint32_t lengthA,
                         std::unordered_map<winnowing::minhash_t, std::vector<int32_t>> &lookup_table,
                         uint32_t s,
                         double tau,
                         std::vector<region> &candidates) {
        uint32_t m = (uint32_t) ceil(tau * s);

        std::vector<int32_t> positions;

        for (auto &minimizer : A_minimizers) {
            auto matches = lookup_table.find(minimizer.hash);
            if (matches == lookup_table.end()) {
                continue;
            }
            for (int32_t pos : matches->second) {
                positions.push_back(pos);
            }
        }

        std::sort(positions.begin(), positions.end());
        int total_positions = (int) positions.size();
        for (int i = 0; i < total_positions - m; ++i) {
            int j = i + m - 1;
            if (positions[j] - positions[i] < lengthA) {
                candidates.emplace_back((region) {positions[j] - lengthA + 1, positions[i]});
            }
        }

    }

    void init_L(std::map<winnowing::minhash_t, bool>& L,
                std::vector<winnowing::minimizer> &minimizers){
        for(auto &minimizer : minimizers){
            L[minimizer.hash] = false;
        }
    }

    winnowing::minhash_t solve_jaccard(std::map<winnowing::minhash_t, bool> L, uint32_t s){
        winnowing::minhash_t shared_sketch = 0
        for(uint32_t i = 0; i < s; i++){
            shared_sketch += L[i];
        }
        return shared_sketch/s;
    }

    void remove_from_map(uint32_t start,
                         uint32_t end,
                         std::map<winnowing::minhash_t, bool> &L,
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
                         std::map<winnowing::minhash_t, bool> &L,
                         std::vector<winnowing::minimizer> &minimizers){
        int i = 0;
        while(minimizers[i].index < start){
            i++;
        }

        while(minimizers[i].index < end){
            L[minimizers[i].hash] = true;
        }
    }



    void compute_estimates(std::vector<winnowing::minimizer> &ref_minimizers,
                           std::vector<winnowing::minimizer> &query_minimizers,
                           uint32_t lengthA,
                           std::vector<region> &candidates,
                           uint32_t s,
                           double tau,
                           std::vector<estimate> estimates){
        for (auto &candidate : candidates) {
            uint32_t  i = candidate.start;
            uint32_t  j = candidate.start + lengthA;
            std::map<winnowing::minhash_t, bool> L;
            init_L(L, query_minimizers);
            winnowing::minhash_t J = solve_jaccard(L, s);
            if (J >= tau){
                estimates.push_back((estimate) {position: i, jaccard: J});
            }
            while(i <= candidate.end){
                remove_from_map(i, i+1, L, ref_minimizers);
                insert_into_map(j, j + 1, L, ref_minimizers);
                J = solve_jaccard(L, s);
                if (J >= tau){
                    estimates.push_back((estimate) {position: i, jaccard: J});
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


    void get_minimizers_between(uint32_t start,
                                uint32_t end,
                                std::vector<winnowing::minimizer> &minimizers,
                                std::vector<winnowing::minhash_t > &result){
        int i = 0;
        while(minimizers[i].index < start){
            i++;
        }

        while(minimizers[i].index < end){
            result.push_back(minimizers[i].hash);
        }
    }
}