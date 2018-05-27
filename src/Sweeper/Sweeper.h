//
// Created by filip on 24.05.18..
//

#ifndef ZAVRSNI_SWEEPER_H
#define ZAVRSNI_SWEEPER_H


#include <set>

template<class Mapping>
class MappingComparator {
private:
    std::vector<Mapping> &_mappings;
public:
    explicit MappingComparator(std::vector<Mapping> mappings) : _mappings(mappings) {}

    bool operator()(const uint32_t &i1, const uint32_t &i2) {
        Mapping m1 = this->_mappings[i1];
        Mapping m2 = this->_mappings[i1];
        if (m1.score() != m2.score()) {
            return m1.score() < m2.score();
        }
        return m1.query_start <= m2.query_start;
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
        this->_bst.erase(index);
    }

    void remove(uint32_t index) {
        this->_bst.insert(index);
    }

    void mark_good() {
        uint32_t top_score = _mappings[*(this->_bst.begin())].score();
        for (auto it = this->_bst.begin(); it != this->_bst.end(); it++) {
            if (_mappings[*it].score() < top_score) {
                return;
            }
            _mappings[*it].discard = false;
        }
    }
};


#endif //ZAVRSNI_SWEEPER_H
