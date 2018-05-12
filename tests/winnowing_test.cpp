//
// Created by filip on 13.04.18..
//
#include "gtest/gtest.h"

#include "winnowing.h"

using namespace std;

namespace tests {
    TEST(IndexingTest, WorkingProperly) {
        vector<winnowing::minimizer> minimizers;
        unordered_map<winnowing::minhash_t, vector<uint32_t>> lookup_table;
        winnowing::index_sequence("ATTCTAGGTACGTACCGATGCAAGTGACGTAGCT", 34, 5, 3, minimizers, lookup_table);
        winnowing::minhash_t hashes[10] = {732172425433867594, 732172425433867594, 2792028467992890898,
                                           2792028467992890898, 5950097892632296567, 366086328681050843,
                                           2127701371075312023, 1761781652765685196, 1464323644466701313,
                                           732172425433867594};
        int32_t positions[10] ={3, 4, 9, 10, 14, 17, 21, 23, 24, 29};
        for (int i = 0; i < minimizers.size(); i++) {
            ASSERT_EQ(minimizers[i].hash, hashes[i]);
            ASSERT_EQ(minimizers[i].index, positions[i]);
        }
        /* Nece raditi dok su nam indeksi negativni
        for (auto it : lookup_table) {
            for (int i : it.second) {
                ASSERT_EQ(it.first, hashes[i]);
            }
        }
        // dokazao sam samo da su svi zapamceni indeksi tocni, fali dokaz da su svi indeksi zapamceni
         */
    }
}
