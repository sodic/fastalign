//
// Created by filip on 13.04.18..
//

#ifndef ZAVRSNI_MAPPER_H
#define ZAVRSNI_MAPPER_H

#include <FASTARead/FASTARead.h>
#include <winnowing/winnowing.h>

namespace mapper {

    const uint32_t DEFAULT_SEG_LENGTH = 5000; // taken from MashMap, default l0
    const double DEFAULT_IDENTITY_THRESHOLD = 0.85; // hopefully taken from MashMap

    typedef struct {
        uint32_t query_seq_length; // treba li nam ovo
        uint32_t query_start;
        uint32_t query_end;
        uint32_t ref_start;
        uint32_t ref_end;
        uint32_t score;
        bool strand;
    } Alignment;

    typedef struct {
        uint32_t mutual;
        int32_t strand;
    } matchInfo;

    typedef struct{
        uint32_t start;
        uint32_t end;
    } region;

    typedef struct{
        uint32_t position;
        bool strand:;
        winnowing::minhash_t jaccard;
    } estimate;

    void find_candidates(std::vector<winnowing::minimizer>& query_minimizers,
                         uint32_t query_length,
                         std::unordered_map<winnowing::minhash_t, std::vector<uint32_t>> &lookup_table,
                         uint32_t s,
                         double tau,
                         std::vector<region> &candidates);

    void compute_estimates(std::vector<winnowing::minimizer> &ref_minimizers,
                           std::vector<winnowing::minimizer> &query_minimizers,
                           uint32_t query_length,
                           std::vector<region> &candidates,
                           uint32_t s,
                           double tau,
                           std::vector<estimate> &estimates);


    void map_fragment(const char *query_seq,
                      uint32_t query_length,
                      std::vector<winnowing::minimizer> &ref_minimizers,
                      std::unordered_map<winnowing::minhash_t, std::vector<uint32_t>> &lookup_table,
                      uint32_t s,
                      double tau,
                      std::vector<estimate> &estimates) {
        // find the minimizers for the query sequence
        std::vector<winnowing::minimizer> query_minimizers;
        winnowing::compute_minimizers(query_seq, query_length, query_minimizers);

        // first stage -> find candidates for jaccard estimations
        std::vector<region> candidates;
        find_candidates(query_minimizers, query_length, lookup_table, s, tau, candidates);

        // second stage -> estimate jaccard for alignments
        compute_estimates(ref_minimizers, query_minimizers, query_length, candidates, s, tau, estimates);
    }


    void process_fragment(const char *Q,
                          uint32_t seg_length,
                          uint32_t offset,
                          std::vector<winnowing::minimizer> &ref_minimizers,
                          std::unordered_map<winnowing::minhash_t, std::vector<std::uint32_t >> &lookup_table,
                          std::vector<Alignment> alignments,
                          uint32_t s,
                          double tau) {
        std::vector<mapper::estimate> mappings;
        mapper::map_fragment(Q + offset,
                             seg_length,
                             ref_minimizers,
                             lookup_table,
                             s,
                             tau,
                             mappings);

        for (mapper::estimate &e : mappings) {
            alignments.push_back(Alignment{
                    query_start : offset,
                    query_end : offset + seg_length - 1,
                    ref_start: e.position,
                    ref_end: e.position + seg_length
            });
        }
    }

    void find_alignments(const char *R,
                         uint32_t r_length,
                         const char *Q,
                         uint32_t q_length,
                         std::vector<Alignment> &alignments) {

        uint32_t s = 6;  // todo rijesi ovo
        double tau = 12;

        // indexing the reference
        std::vector<winnowing::minimizer> ref_minimizers;
        std::unordered_map<winnowing::minhash_t, std::vector<std::uint32_t >> lookup_table;
        winnowing::index_sequence(R, r_length, ref_minimizers, lookup_table);

        // map each of the l0/2 fragments of the query
        uint32_t seg_length = DEFAULT_SEG_LENGTH;
        uint32_t number_of_fragments = r_length / seg_length;
        for (uint32_t i = 0; i < number_of_fragments; ++i) {
            process_fragment(Q, seg_length, i * seg_length, ref_minimizers, lookup_table, alignments, s, tau);

        }

        // check if we have one last fragment left
        if (number_of_fragments > 0 && r_length % seg_length != 0) {
            process_fragment(Q, seg_length, q_length - seg_length, ref_minimizers, lookup_table, alignments, s, tau);
        }


    }

}

#endif //ZAVRSNI_MAPPER_H
