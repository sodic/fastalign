//
// Created by filip on 13.04.18..
//
#include <winnowing_dac/winnowing.hpp>
#include "gtest/gtest.h"


using namespace std;

namespace winnowing {
    TEST(minimizerTest, minimizerTest) {
        std::vector<winnowing::minimizer> result;
        winnowing::compute_minimizers("ACTGATGC", 8, 2, 3, result);
        std::vector<winnowing::minimizer> expected = {winnowing::minimizer(9, -1, 1,1), winnowing::minimizer(30, -1,1, 1),
                                                      winnowing::minimizer(14, -1,1, -1)};
        ASSERT_EQ((int) result.size(), (int) expected.size());
        for (int i = 0; i < result.size(); i++) {
            ASSERT_EQ(result[i].hash, expected[i].hash);
            ASSERT_EQ(result[i].strand, expected[i].strand);
        }

    }
}
