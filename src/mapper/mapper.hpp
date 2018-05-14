//
// Created by filip on 13.04.18..
//

#ifndef ZAVRSNI_MAPPER_H
#define ZAVRSNI_MAPPER_H

#include <FASTARead/FASTARead.h>
#include <winnowing/winnowing.hpp>

namespace mapper {

    const uint32_t DEFAULT_FRAGMENT_LENGTH = 5000; // taken from MashMap, default l0/2
    const double DEFAULT_IDENTITY_THRESHOLD = 0.85; // hopefully taken from MashMap

    double jaccard_to_e(float jaccard, int k) {
        if (jaccard == 0) {
            return 1.0; //jaccard estimate 0 -> 1.0 mash distance
        }

        if (jaccard == 1) {
            return 0.0; //jaccard estimate 1 -> 0.0 mash distance
        }

        return (-1.0 / k) * log(2.0 * jaccard / (1 + jaccard));
    }

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
    } Mapping;

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
        bool strand;
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
                          uint32_t fragment_length,
                          uint32_t offset,
                          std::vector<winnowing::minimizer> &ref_minimizers,
                          std::unordered_map<winnowing::minhash_t, std::vector<std::uint32_t >> &lookup_table,
                          std::vector<Mapping> &mappings,
                          uint32_t s,
                          double tau) {
        std::vector<mapper::estimate> estimates;
        mapper::map_fragment(Q + offset,
                             fragment_length,
                             ref_minimizers,
                             lookup_table,
                             s,
                             tau,
                             estimates);

        for (mapper::estimate &e : estimates) {
            Mapping m{};
            m.query_start = offset;
            m.query_end = offset + fragment_length - 1;
            m.ref_start = e.position;
            m.ref_end = e.position + fragment_length;
            m.identity_estimate = 100 * (1 - jaccard_to_e(e.jaccard, winnowing::DEFAULT_K));
            mappings.push_back(m);
        }
    }


    void merge_mappings(std::vector<Mapping> &mappings, uint32_t fragment_length) {
        if (mappings.size() < 2) {
            return;
        }

        // sort mappings by positions on the reference
        std::sort(mappings.begin(), mappings.end(), [](const Mapping a, const Mapping b) {
            return std::tie(a.ref_start, a.query_start) < std::tie(b.ref_start, b.ref_end);
        });

        // assign sequential IDs to all mappings
        for (uint32_t i = 0, num_of_mappings = mappings.size(); i < num_of_mappings; ++i) {
            mappings[i].split_id = i;
        }

        for (auto it = mappings.begin(); it != mappings.end();) {
            double current_frag_index = std::ceil(
                    it->query_start * 1.0 / fragment_length); // todo nisam siguran cemu ovo, mogu li se samo dijeliti
            for (auto it_chain = std::next(it); it_chain != mappings.end(); it_chain++) {

                // if fragments are further apart than l0/2, no chaining
                if (it_chain->ref_start - it->ref_end > 2 * fragment_length) {
                    break;
                }

                double chain_frag_index = std::ceil(
                        it_chain->query_start * 1.0 / fragment_length);  // todo ista stvar, vidi gore
                if (it->strand == it_chain->strand
                    && chain_frag_index == current_frag_index + (it->strand ? 1 : -1)) {
                    it_chain->split_id = it->split_id;
                }

            }
        }

        // sort mappings by split IDs
        std::sort(mappings.begin(), mappings.end(), [](const Mapping a, const Mapping b) {
            return a.split_id < b.split_id;
        });

        // each chain should be represented by only one mapping
        auto it = mappings.begin();
        while (it != mappings.end()) {
            auto it_chain_end = std::find_if(it, mappings.end(), [&it](Mapping m) {
                return m.split_id != it->split_id;
            });

            // go through the entire chain and collect data into the first mapping
            std::for_each(it, it_chain_end, [&it](Mapping &m) {
                it->query_start = std::min(it->query_start, m.query_start);
                it->ref_start = std::min(it->ref_start, m.ref_start);

                it->query_end = std::max(it->query_end, m.query_end);
                it->ref_end = std::max(it->ref_end, m.ref_end);
            });

            // the identity is set to be the average of all mappings in the chain
            it->identity_estimate = std::accumulate(it, it_chain_end, 0.0, [](double x, Mapping &m) {
                return x + m.identity_estimate;
            });
            it->identity_estimate = it->identity_estimate / std::distance(it, it_chain_end);

            // mark all other mappings of the chain for removal
            std::for_each(std::next(it), it_chain_end, [](Mapping &m) {
                m.discard = true;
            });

            // next mapping is the first not included in the chain
            it = it_chain_end;
        }

        // remove all unneeded mappings
        mappings.erase(std::remove_if(
                mappings.begin(), mappings.end(),
                [](Mapping &m) {
                    return m.discard;
                }), mappings.end());
    }

    void compute_mappings(const char *R,
                          uint32_t r_length,
                          const char *Q,
                          uint32_t q_length,
                          std::vector<Mapping> &mappings) {

        uint32_t s = 6;  // todo rijesi ovo
        double tau = 12;

        // indexing the reference
        std::vector<winnowing::minimizer> ref_minimizers;
        std::unordered_map<winnowing::minhash_t, std::vector<std::uint32_t >> lookup_table;
        winnowing::index_sequence(R, r_length, ref_minimizers, lookup_table);

        // map each of the l0/2 fragments of the query
        uint32_t fragment_length = DEFAULT_FRAGMENT_LENGTH;
        uint32_t number_of_fragments = r_length / fragment_length;
        for (uint32_t i = 0; i < number_of_fragments; ++i) {
            process_fragment(Q, fragment_length, i * fragment_length, ref_minimizers, lookup_table, mappings, s, tau);

        }

        // check if we have one last fragment left
        if (number_of_fragments > 0 && r_length % fragment_length != 0) {
            process_fragment(Q, fragment_length, q_length - fragment_length, ref_minimizers, lookup_table, mappings, s,
                             tau);
        }

        merge_mappings(mappings, fragment_length);
    }

}

#endif //ZAVRSNI_MAPPER_H
