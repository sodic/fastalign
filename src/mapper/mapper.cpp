//
// Created by filip on 13.04.18..
//

#include <cstdint>
#include <vector>
#include <FASTARead/FASTARead.h>
#include <winnowing/winnowing.h>
#include<cmath>
#include <algorithm>
#include "mapper.h"


namespace mapper {

    void find_candidates(const char *sequenceA,
                         uint32_t lengthA,
                         std::unordered_map<winnowing::minhash_t, std::vector<int32_t>> lookup_table,
                         uint32_t s,
                         double tau) {
        uint32_t m = (uint32_t) ceil(tau * s);

        std::vector<winnowing::minimizer> minimizersA;
        std::vector<int32_t> positions;
        winnowing::compute_minimizers(sequenceA, lengthA, winnowing::DEFAULT_W, winnowing::DEFAULT_K, minimizersA);

        for (auto &minimizer : minimizersA) {
            auto matches = lookup_table.find(minimizer.hash);
            if (matches == lookup_table.end()) {
                continue;
            }
            for (int32_t pos : matches->second) {
                positions.push_back(pos);
            }
        }

        std::sort(positions.begin(), positions.end());
        std::vector<region> candidates;
        int total_positions = (int) positions.size();
        for (int i = 0; i < total_positions - m; ++i) {
            int j = i + m - 1;
            if (positions[j] - positions[i] < positions.size()) {
                candidates.emplace_back((region) {positions[j] - total_positions + 1, positions[i]});
            }
        }

    }
}