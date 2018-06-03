//
// Created by filip on 13.04.18..
//

#ifndef ZAVRSNI_MAPPER_H
#define ZAVRSNI_MAPPER_H

#include <FASTARead/FASTARead.h>
#include "winnowing_dac/winnowing.hpp"
#include "statistics/statistics.hpp"
#include "Sweeper/Sweeper.h"

namespace mapper {

    typedef struct {
        uint32_t query_seq_length; // treba li nam ovo
        uint32_t query_start;
        uint32_t query_end;
        uint32_t ref_start;
        uint32_t ref_end;
        double identity_estimate;
        bool strand;

        uint32_t split_id; // split reading combinations
        bool discard;

        void setInvalid() {
            this->discard = false;
        }

        double score() {
            return this->identity_estimate * (this->query_end - this->query_start + 1);
        }

        uint32_t start() {
            return this->query_start;
        }

        void mark_good() {
            this->discard = false;
        }

        void mark_redundant() {
            this->discard = true;
        }
    } Mapping;


    void compute_mappings(const char *R,
                          uint32_t r_length,
                          const char *Q,
                          uint32_t q_length,
                          std::vector<Mapping> &mappings);

}

#endif //ZAVRSNI_MAPPER_H
