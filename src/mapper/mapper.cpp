//
// Created by filip on 13.04.18..
//

#include <cstdint>
#include <vector>
#include "FASTARead/FASTARead.h"
#include<cmath>
#include <algorithm>
#include <map>
#include <assert.h>
#include <config/config.hpp>
#include <iostream>
#include "mapper.hpp"


namespace mapper {


    typedef struct {
        uint32_t mutual;
        int32_t strand;
    } matchInfo;

    typedef struct {
        uint32_t start;
        uint32_t end;
    } region;

    typedef struct {
        uint32_t position;
        bool strand;
        double identity;
    } estimate;

    void find_candidates(std::vector<winnowing::minimizer>& query_minimizers,
                         uint32_t query_length,
                         std::unordered_map<winnowing::minhash_t, std::vector<uint32_t>> &lookup_table,
                         uint32_t s,
                         double tau,
                         std::vector<region> &candidates) {
        int m = (uint32_t) ceil(0.001 * tau * s);

        int alt_m = statistics::min_hits_relaxed_estimate(s, config::constants::default_k);
        alt_m = alt_m < 1 ? 1 : alt_m;

        std::vector<uint32_t> positions;         // the array L

        for (auto &minimizer : query_minimizers) {
            auto matches = lookup_table.find(minimizer.hash);
            if (matches == lookup_table.end()) {
                continue;
            }
            for (uint32_t pos : matches->second) {
                positions.push_back(pos);
            }
        }

        std::sort(positions.begin(), positions.end());

        int total_positions = (int) positions.size();
        for (int i = 0; i <= total_positions - alt_m; ++i) {
            int j = i + m - 1;

            // check if consecutive hits are close enough
            if (positions[j] - positions[i] < query_length) {

                region candidate = {
                        start: std::max(0u, positions[j] - query_length + 1),
                        end: positions[i]
                };

                //check overlap with the previous candidate
                auto prev_it = candidates.end();
                prev_it--;

                if (!candidates.empty() && prev_it->end >= candidate.start) {
                    // if there is an overlap, extend the previous candidate
                    prev_it->end = std::max(candidate.end, prev_it->end);
                } else {
                    candidates.push_back(candidate);
                }
            }

        }
    }


    void init_L(std::map<winnowing::minhash_t, matchInfo> &L,
                std::vector<winnowing::minimizer> &minimizers){
        for(auto &minimizer : minimizers){
            L[minimizer.hash] = (matchInfo) {
                    mutual: 0,
                    strand: minimizer.strand
            };
        }
    }

    double solve_jaccard(std::map<winnowing::minhash_t, matchInfo> &L, uint32_t s) {
        assert(s > 0);

        uint32_t shared_sketch = 0;
        auto it = L.begin();
        uint32_t collected = 0;
        while (it != L.end() && collected < s) {
            shared_sketch += it->second.mutual;
            it++;
            collected++;
        }
        return 1.0 * shared_sketch / s;
    }


    void remove_from_map(uint32_t start,
                         uint32_t end,
                         std::map<winnowing::minhash_t, matchInfo> &L,
                         std::vector<winnowing::minimizer> &minimizers){
        if (minimizers.empty()) {
            return;
        }
        int i = 0;
        while (i < minimizers.size() && minimizers[i].index < start) {
            i++;
        }

        while (i < minimizers.size() && minimizers[i].index < end) {
            L.erase(minimizers[i].hash);
            i++;
        }
    }


    void insert_into_map(uint32_t start,
                         uint32_t end,
                         std::map<winnowing::minhash_t, matchInfo> &L,
                         std::vector<winnowing::minimizer> &minimizers){
        if (minimizers.empty()) {
            return;
        }
        int i = 0;
        while (i < minimizers.size() && minimizers[i].index < start) {
            i++;
        }

        while (i < minimizers.size() && minimizers[i].index < end) {
            if (L.find(minimizers[i].hash) != L.end()) {
                L[minimizers[i].hash].mutual = 1;
                L[minimizers[i].hash].strand *= minimizers[i].strand;
            } else {
                L[minimizers[i].hash] = matchInfo {
                        mutual: 0,
                        strand: minimizers[i].strand
                };
            }
            i++;
        }
    }

    bool consensus_strand(std::map<winnowing::minhash_t, matchInfo> L, uint32_t s) {
        assert(s > 0);

        int32_t sum = 0;
        uint32_t counter = 0;
        auto it = L.begin();
        while (it != L.end() && counter < s) {
            if (it->second.mutual) {
                sum += it->second.strand;
                counter++;
            }
            it++;
        }
        return sum >= 0;
    }



    void compute_estimates(std::vector<winnowing::minimizer> &ref_minimizers,
                           std::vector<winnowing::minimizer> &query_minimizers,
                           uint32_t query_length,
                           std::vector<region> &candidates,
                           uint32_t s,
                           double tau,
                           std::vector<estimate> &estimates){
        for (auto &candidate : candidates) {
            uint32_t  i = candidate.start;
            uint32_t  j = candidate.start + query_length;

            std::map<winnowing::minhash_t, matchInfo> L;
            init_L(L, query_minimizers);                // treba ih i za ostale opet postavit na nulu, ili? mislim da ne
            insert_into_map(i, j, L, ref_minimizers);

            double J = solve_jaccard(L, s);
            double e = statistics::jaccard_to_e(J, 16);

            //Compute lower bound to mash distance within 90% confidence interval
            double e_lower_bound = statistics::e_lower_bound(e, s, 16, config::constants::confidence_interval);

            double identity = 100 * (1 - e);
            double identity_upper_boud = 100 * (1 - e_lower_bound);
            if (identity_upper_boud >= tau) {
                estimates.push_back((estimate) {
                        position: i,
                        strand: consensus_strand(L, s),
                        identity:  identity
                });
            }

            while(i <= candidate.end){
                L.erase(i);
                remove_from_map(i, i+1, L, ref_minimizers);
                insert_into_map(j, j + 1, L, ref_minimizers);

                J = solve_jaccard(L, s);
                e = statistics::jaccard_to_e(J, 16);

                //Compute lower bound to mash distance within 90% confidence interval
                e_lower_bound = statistics::e_lower_bound(e, s, 16, config::constants::confidence_interval);

                identity = 100 * (1 - e);
                identity_upper_boud = 100 * (1 - e_lower_bound);
                if (identity_upper_boud >= tau) {
                    estimates.push_back((estimate) {
                            position: i,
                            strand: consensus_strand(L, s),
                            identity: identity
                    });
                }

                i++;
                j++;
            }
        }
    }


    enum Event {
        Begin, End
    };
    typedef std::tuple<uint32_t, Event, uint32_t> sweeppoint;

    void filter_on_query(std::vector<Mapping> &mappings) {
        //assume all mappings are redundant
        for (Mapping &m : mappings) {
            m.mark_redundant();
        }

        //preparing the 2n points for the sweep
        std::vector<sweeppoint> points;
        for (int i = 0; i < mappings.size(); ++i) {
            points.emplace_back(mappings[i].query_start, Begin, i);
            points.emplace_back(mappings[i].query_end + 1, End, i);
        }
        std::sort(points.begin(), points.end());

        Sweeper<mapper::Mapping> sweeper(mappings);

        auto it = points.begin();
        while (it != points.end()) {
            uint32_t current_pos = std::get<0>(*it);
            while (std::get<0>(*it) == current_pos) {

                if (std::get<1>(*it) == Event::Begin) {
                    sweeper.insert(std::get<2>(*it));
                } else {
                    sweeper.remove(std::get<2>(*it));
                }
                it++;
            }
            sweeper.mark_good();
        }

        mappings.erase(std::remove_if(mappings.begin(), mappings.end(), [](const Mapping &m) {
            return m.discard;
        }), mappings.end());
    }


    void map_fragment(const char *query_seq,
                      uint32_t query_length,
                      std::vector<winnowing::minimizer> &ref_minimizers,
                      std::unordered_map<winnowing::minhash_t, std::vector<uint32_t>> &lookup_table,
                      double tau,
                      std::vector<estimate> &estimates) {
        // find the minimizers for the query sequence
        std::vector<winnowing::minimizer> query_minimizers;
        winnowing::compute_minimizers(query_seq, query_length, 111, 16, query_minimizers);

        // sketch size is |Wh(a)|
        uint32_t s = query_minimizers.size();
        if (s == 0) {
            return;
        }

        // first stage -> find candidates for identity estimations
        std::vector<region> candidates;
        find_candidates(query_minimizers, query_length, lookup_table, s, tau, candidates);

        // second stage -> estimate identity for alignments
        compute_estimates(ref_minimizers, query_minimizers, query_length, candidates, s, tau, estimates);
    }


    void process_fragment(const char *Q,
                          uint32_t fragment_length,
                          uint32_t offset,
                          std::vector<winnowing::minimizer> &ref_minimizers,
                          std::unordered_map<winnowing::minhash_t, std::vector<std::uint32_t >> &lookup_table,
                          std::vector<Mapping> &mappings,
                          double tau) {
        std::vector<mapper::estimate> estimates;
        mapper::map_fragment(Q + offset,
                             fragment_length,
                             ref_minimizers,
                             lookup_table,
                             tau,
                             estimates);

        for (mapper::estimate &e : estimates) {
            Mapping m{};
            m.query_start = offset;
            m.query_end = offset + fragment_length - 1;
            m.ref_start = e.position;
            m.ref_end = e.position + fragment_length;
            m.identity_estimate = e.identity;
            m.strand = e.strand;
            mappings.push_back(m);
        }
    }


    void merge_mappings(std::vector<Mapping> &mappings, uint32_t fragment_length) {
        if (mappings.size() < 2) {
            return;
        }

        // sort mappings by positions on the reference
        std::sort(mappings.begin(), mappings.end(), [](const Mapping a, const Mapping b) {
            return std::tie(a.ref_start, a.query_start) < std::tie(b.ref_start, b.query_start);
        });

        // assign sequential IDs to all mappings
        for (uint32_t i = 0, num_of_mappings = mappings.size(); i < num_of_mappings; ++i) {
            mappings[i].split_id = i;
        }

        for (auto it = mappings.begin(); it != mappings.end(); it++) {
            auto current_frag_index = std::ceil(
                    it->query_start * 1.0 / fragment_length); // which fragment is this wrt. the whole query sequence
            for (auto it_chain = std::next(it); it_chain != mappings.end(); it_chain++) {

                // if fragments are further apart than l0/2, no chaining
                if (it_chain->ref_start - it->ref_end > 2 * fragment_length) {
                    break;
                }

                auto chain_frag_index = std::ceil(
                        it_chain->query_start * 1.0 / fragment_length);
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


    int binary_search(std::vector<winnowing::minimizer> &minimizers,
                      uint32_t start,
                      uint32_t length,
                      uint32_t element){
        if(!length){
            return -1;
        }

        uint32_t middle = start + length/2;

        if (minimizers[middle].index == element){
            return middle;
        }

        if (minimizers[middle].index < element) {
            return binary_search(minimizers, middle + 1, length/2 - 1 + length%2, element);
        }

        return binary_search(minimizers, start, length/2, element);
    }

    void compute_mappings(const char *R,
                          uint32_t r_length,
                          const char *Q,
                          uint32_t q_length,
                          std::vector<Mapping> &mappings) {

        double tau = config::constants::default_percentage_identity;

        // indexing the reference
        std::vector<winnowing::minimizer> ref_minimizers;
        std::unordered_map<winnowing::minhash_t, std::vector<std::uint32_t >> lookup_table;
        auto k = config::constants::default_k;
        auto w = statistics::calculate_window_size(k, config::constants::default_segment_length, r_length);
        winnowing::index_sequence(R, r_length, w, k, ref_minimizers, lookup_table);

        // map each of the l0/2 fragments of the query
        uint32_t fragment_length = config::constants::default_segment_length;
        uint32_t number_of_fragments = r_length / fragment_length;

        std::cout << "Fragment count: " << number_of_fragments << "(+/- 1)" << std::endl;

        for (uint32_t i = 0; i < number_of_fragments; ++i) {
            std::cout << i << std::endl;
            process_fragment(Q, fragment_length, i * fragment_length, ref_minimizers, lookup_table, mappings, tau);

        }

        // check if we have one last fragment left
        if (number_of_fragments > 0 && r_length % fragment_length != 0) {
            process_fragment(Q, fragment_length, q_length - fragment_length, ref_minimizers, lookup_table, mappings,
                             tau);
        }

        std::cout << "All fragments processed." << std::endl;

        merge_mappings(mappings, fragment_length);

        std::cout << "Merged all found mappings, preparing to filter." << std::endl;

        //linear sweep
        filter_on_query(mappings);
    }


}
