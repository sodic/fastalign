//
// Created by filip on 24.05.18..
//

#ifndef ZAVRSNI_SWEEPER_H
#define ZAVRSNI_SWEEPER_H


#include <set>
#include <vector>
#include <cstdint>
#include <tuple>

template<class Mapping>
class MappingComparator {
private:
    std::vector<Mapping> &_mappings;
public:
    explicit MappingComparator(std::vector<Mapping> &mappings) : _mappings(mappings) {}

    bool operator()(const uint32_t &i1, const uint32_t &i2) {
        Mapping m1 = this->_mappings[i1];
        double score1 = m1.score();
        double start1 = m1.start();

        Mapping m2 = this->_mappings[i2];
        double score2 = m2.score();
        double start2 = m2.start();

        return std::tie(score1, start1) > std::tie(score2, start2);
    }
};


template<class Mapping>
class Sweeper {
private:
    std::set<uint32_t, MappingComparator<Mapping>> _bst;
    std::vector<Mapping> &_mappings;
public:
    explicit Sweeper(std::vector<Mapping> &mappings) : _mappings(mappings),
                                                       _bst(MappingComparator<Mapping>(mappings)) {
    }

    void insert(uint32_t index) {
        this->_bst.insert(index);
    }

    void remove(uint32_t index) {
        this->_bst.erase(index);
    }

    uint32_t pop() {
        uint32_t index = *(this->_bst.begin());
        this->_bst.erase(index);
        return index;
    }

    void mark_good() {
        auto it = this->_bst.begin();
        uint32_t top_score = _mappings[*it].score();
        while (it != this->_bst.end() && _mappings[*it].score() >= top_score) {
            _mappings[*it].mark_good();
            it++;
        }
    }
};


#endif //ZAVRSNI_SWEEPER_H
