//
// Created by filip on 29.05.18..
//

#include <statistics/statistics.hpp>
#include "gtest/gtest.h"

namespace statistics {
    TEST(JacccardToETest, WorkingPoperly) {
        ASSERT_NEAR(0.101699, statistics::jaccard_to_e(0.43, 5), 0.000001);
        ASSERT_NEAR(0.019491, statistics::jaccard_to_e(0.56, 17), 0.000001);
        ASSERT_NEAR(0.098354, statistics::jaccard_to_e(0.23, 10), 0.000001);

        ASSERT_DOUBLE_EQ(0, statistics::jaccard_to_e(1, 23));
        ASSERT_DOUBLE_EQ(1, statistics::jaccard_to_e(0, 23));
    }

    TEST(EToJaccardTest, WorkingPoperly) {
        ASSERT_NEAR(0.177298, statistics::e_to_jaccard(0.1, 12), 0.000001);
        ASSERT_NEAR(0.047514, statistics::e_to_jaccard(0.24, 10), 0.000001);
        ASSERT_NEAR(0.006561, statistics::e_to_jaccard(0.62, 7), 0.000001);

        ASSERT_DOUBLE_EQ(1, statistics::jaccard_to_e(0, 1));
    }

    TEST(ELowerBoundTest, WorkingPoperly) {
        ASSERT_NEAR(0.0610763, statistics::e_lower_bound(0.1, 12, 5, 0.78), 0.000001);
        ASSERT_NEAR(0.186539, statistics::e_lower_bound(0.3, 10, 3, 0.8), 0.000001);
        ASSERT_NEAR(0.297063, statistics::e_lower_bound(0.4, 15, 7, 0.74), 0.000001);
    }

    TEST(EstimtePValueTest, WorkingPoperly) {
        ASSERT_NEAR(27.356028755725401, statistics::estimate_p_value(10, 16, 5000, 4699745), 0.0001);
    }

    TEST(RecommendedWindowSizeTest, WorkingPoperly) {
        ASSERT_EQ(111, statistics::calculate_window_size(16, 5000, 4699745));
        ASSERT_EQ(125, statistics::calculate_window_size(20, 5000, 3000000));
        ASSERT_EQ(100, statistics::calculate_window_size(21, 5000, 5000000));
    }
}
