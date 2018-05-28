//
// Created by filip on 13.04.18..
//
#include <Sweeper/Sweeper.h>
#include "gtest/gtest.h"

#include "winnowing.hpp"

using namespace std;

namespace sweeper {

    class MockMapping {
        double _score;
        uint32_t _start;
    public:
        MockMapping(uint32_t start, double score) : _start(start), _score(score) {}

        double score() {
            return _score;
        }

        uint32_t start() {
            return _start;
        }
    };

    TEST(SweeperTest, ProperSorting) {
        std::vector<MockMapping> mappings;
        mappings.emplace_back(MockMapping(2, 13.0));
        mappings.emplace_back(MockMapping(6, 20.0));
        mappings.emplace_back(MockMapping(4, 20.0));
        mappings.emplace_back(MockMapping(3, 15.0));
        mappings.emplace_back(MockMapping(6, 21.5));
        mappings.emplace_back(MockMapping(1, 17.0));
        mappings.emplace_back(MockMapping(6, 21.0));
        mappings.emplace_back(MockMapping(7, 25.0));
        mappings.emplace_back(MockMapping(0, 14.0));
        mappings.emplace_back(MockMapping(10, 17.0));
        mappings.emplace_back(MockMapping(7, 24.0));
        Sweeper<MockMapping> sweeper(mappings);
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            sweeper.insert(i);
        }
        std::vector<MockMapping> sorted = mappings;
        std::sort(sorted.begin(), sorted.end(), [](MockMapping &m1, MockMapping &m2) {
            if (m1.score() == m2.score()) {
                return m1.start() >= m2.start();
            }
            return m1.score() > m2.score();
        });
        for (MockMapping mapping : sorted) {
            uint32_t idx = sweeper.pop();
            MockMapping from_sweeper = mappings[idx];
            ASSERT_EQ(mapping.start(), from_sweeper.start());
            ASSERT_EQ(mapping.score(), from_sweeper.score());
        }
    }
}
