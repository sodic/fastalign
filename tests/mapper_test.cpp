//
// Created by filip on 14.05.18..
//

#include <mapper/mapper.cpp>
#include "gtest/gtest.h"

using namespace std;
using namespace mapper;

namespace winnowing {


    void init_sets(set<uint32_t> &s1, set<uint32_t> &s2) {
        srand(time(NULL));
        for (int i = 0; i < 100; ++i) {
            s1.insert(rand() % 1000);
            s2.insert(rand() % 1000);
        }
    }


    TEST(SolveJaccardTest, WorkingPoperly) {
        std::map<winnowing::minhash_t, matchInfo> L = {
                {1, matchInfo {mutual: 0, strand: 1}},
                {2, matchInfo {mutual: 1, strand: 1}},
                {3, matchInfo {mutual: 0, strand: 1}},
                {4, matchInfo {mutual: 0, strand: 1}},
                {5, matchInfo {mutual: 1, strand: 1}},
                {6, matchInfo {mutual: 0, strand: 1}}

        };
        uint32_t s = 4;
        ASSERT_EQ(solve_jaccard(L, 2), 0.5);
    }

    TEST(SolveJaccardTest, PerfectEstimate) {
        set<uint32_t> s1;
        set<uint32_t> s2;
        init_sets(s1, s2);

        set<uint32_t> union_set;
        set<uint32_t> intersection_set;

        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                         inserter(intersection_set, intersection_set.begin()));
        set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(union_set, union_set.begin()));

        std::map<winnowing::minhash_t, matchInfo> L;
        for (auto &e : union_set) {
            L[e] = mapper::matchInfo{0, 0};
        }

        for (auto &e : intersection_set) {
            L[e].mutual = 1;
        }

        ASSERT_EQ(1.0 * intersection_set.size() / union_set.size(), solve_jaccard(L, union_set.size()));
    }

    TEST(SolveJaccardTest, STooBig) {
        std::map<winnowing::minhash_t, matchInfo> L;
        srand(time(NULL));
        for (int i = 0; i < 100; i++) {
            L[i] = matchInfo {mutual: (uint32_t) rand() % 1, strand: rand() % 2 ? -1 : 1};
        }

        ASSERT_EQ(solve_jaccard(L, L.size() + 200), solve_jaccard(L, L.size()));
    }

    TEST(SolveJaccardTest, SIs0) {
        std::map<winnowing::minhash_t, matchInfo> L;
        for (int i = 0; i < 100; i++) {
            L[i] = matchInfo {mutual: 1, strand: 1};
        }

        ASSERT_DEATH(solve_jaccard(L, 0), "");
    }


    TEST(ConsensusStrandTest, WorkingPoperly) {
        std::map<winnowing::minhash_t, matchInfo> L = {
                {1,  matchInfo {mutual: 0, strand: -1}},
                {2,  matchInfo {mutual: 1, strand: -1}},
                {3,  matchInfo {mutual: 0, strand: 1}},
                {4,  matchInfo {mutual: 0, strand: 1}},
                {5,  matchInfo {mutual: 1, strand: -1}},
                {6,  matchInfo {mutual: 1, strand: 1}},
                {7,  matchInfo {mutual: 1, strand: 1}},
                {8,  matchInfo {mutual: 0, strand: -1}},
                {9,  matchInfo {mutual: 1, strand: 1}},
                {10, matchInfo {mutual: 1, strand: -1}},
                {11, matchInfo {mutual: 1, strand: -1}},

        };

        ASSERT_EQ(consensus_strand(L, 2), false);
        ASSERT_EQ(consensus_strand(L, 3), false);
        ASSERT_EQ(consensus_strand(L, 4), true);
        ASSERT_EQ(consensus_strand(L, 5), true);
        ASSERT_EQ(consensus_strand(L, 7), false);
    }

    TEST(ConsensusStrandTest, STooBig) {
        std::map<winnowing::minhash_t, matchInfo> L;
        srand(time(NULL));
        for (int i = 0; i < 100; i++) {
            L[i] = matchInfo {mutual: (uint32_t) rand() % 1, strand: rand() % 2 ? -1 : 1};
        }

        ASSERT_EQ(consensus_strand(L, L.size() + 10), consensus_strand(L, L.size()));
    }

    TEST(ConsensusStrandTest, SIs0) {
        std::map<winnowing::minhash_t, matchInfo> L;
        for (int i = 0; i < 100; i++) {
            L[i] = matchInfo {mutual: 1, strand: 1};
        }

        ASSERT_DEATH(consensus_strand(L, 0), "");
    }

    TEST(InsertIntoMap, WorkingProperly) {
        std::vector<winnowing::minimizer> minimizers = {
                winnowing::minimizer {hash: 1, index: 34, strand: 1},
                winnowing::minimizer {hash: 2, index: 56, strand: 1},
                winnowing::minimizer {hash: 3, index: 78, strand: 1},
                winnowing::minimizer {hash: 4, index: 93, strand: 1},
                winnowing::minimizer {hash: 5, index: 123, strand: 1},
                winnowing::minimizer {hash: 6, index: 156, strand: 1},
                winnowing::minimizer {hash: 7, index: 167, strand: 1},
                winnowing::minimizer {hash: 8, index: 169, strand: 1},
                winnowing::minimizer {hash: 8, index: 171, strand: 1},
                winnowing::minimizer {hash: 9, index: 185, strand: 1},
        };

        std::map<winnowing::minhash_t, matchInfo> L;
        insert_into_map(75, 170, L, minimizers);
        ASSERT_EQ(L.size(), 6);
        for (int i = 3; i <= 8; ++i) {
            ASSERT_NE(L.find(i), L.end());
            ASSERT_EQ(L.find(i)->second.mutual, 0);
        }

        insert_into_map(75, 170, L, minimizers);
        ASSERT_EQ(L.size(), 6);
        for (int i = 3; i <= 8; ++i) {
            ASSERT_EQ(L.find(i)->second.mutual, 1);
        }
    }

    TEST(RemoveFromMap, WorkingProperly) {
        std::vector<winnowing::minimizer> minimizers = {
                winnowing::minimizer {hash: 1, index: 34, strand: 1},
                winnowing::minimizer {hash: 2, index: 56, strand: 1},
                winnowing::minimizer {hash: 3, index: 78, strand: 1},
                winnowing::minimizer {hash: 4, index: 93, strand: 1},
                winnowing::minimizer {hash: 5, index: 123, strand: 1},
                winnowing::minimizer {hash: 6, index: 156, strand: 1},
                winnowing::minimizer {hash: 7, index: 167, strand: 1},
                winnowing::minimizer {hash: 8, index: 169, strand: 1},
                winnowing::minimizer {hash: 8, index: 171, strand: 1},
                winnowing::minimizer {hash: 9, index: 185, strand: 1},
        };

        std::map<winnowing::minhash_t, matchInfo> L;
        insert_into_map(75, 170, L, minimizers);
        ASSERT_EQ(L.size(), 6);
        for (int i = 3; i <= 8; ++i) {
            ASSERT_NE(L.find(i), L.end());
            ASSERT_EQ(L.find(i)->second.mutual, 0);
        }

        insert_into_map(75, 170, L, minimizers);
        ASSERT_EQ(L.size(), 6);
        for (int i = 3; i <= 8; ++i) {
            ASSERT_EQ(L.find(i)->second.mutual, 1);
        }
    }

    /*TEST(SolveJaccardTest, OKEstimate) {
        set<uint32_t> s1;
        set<uint32_t> s2;
        init_sets(s1, s2);

        set<uint32_t> union_set;
        set<uint32_t> intersection_set;

        set_intersection(s1.begin(), s2.end(), s2.begin(), s2.end(), inserter(intersection_set, intersection_set.begin()));
        set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(union_set, union_set.begin()));


        std::map<winnowing::minhash_t, matchInfo> L;
        for(auto &e : union_set){
            L[e] = mapper::matchInfo{0, 0};
        }

        for(auto &e : intersection_set){
            L[e].mutual = 1;
        }

        ASSERT_EQ(1.0*intersection_set.size()/union_set.size(), solve_jaccard(L, union_set.size()/1.3));
        }*/
}
