//
// Created by filip on 29.05.18..
//

#include <algorithm>
#include <config/config.hpp>
#include <boost/math/distributions/binomial.hpp>

namespace statistics {
    using namespace ::boost::math;

    double jaccard_to_e(double jaccard, int k) {
        if (jaccard == 0) {
            return 1.0; //identity estimate 0 -> 1.0 mash distance
        }

        if (jaccard == 1) {
            return 0.0; //identity estimate 1 -> 0.0 mash distance
        }

        return (-1.0 / k) * log(2.0 * jaccard / (1 + jaccard));
    }

    double e_to_jaccard(float e, int k) {
        return 1.0 / (2 * exp(e * k) - 1);
    }

    double e_lower_bound(double e,
                         unsigned int sketch_size,
                         unsigned int k,
                         double confidence_interval) {
        double q2 = (1.0 - confidence_interval) / 2;
        double x = quantile(complement(binomial(sketch_size, e_to_jaccard(e, k)), q2));
        double j = x / sketch_size;
        return jaccard_to_e(j, k);
    }

    unsigned int minimum_hits_estimate(unsigned int sketch_size,
                                       unsigned int k) {
        double e = 1 - config::constants::default_percentage_identity / 100.0;
        double j = e_to_jaccard(e, k);

        return ceil(sketch_size * j);
    }

    unsigned int min_hits_relaxed_estimate(unsigned int sketch_size,
                                           unsigned int k) {
        unsigned int upper_limit = minimum_hits_estimate(sketch_size, k);
        unsigned int min_hits_relaxed = upper_limit;
        for (unsigned i = upper_limit; i > 0; --i) {
            double j = 1.0 * i / sketch_size;
            double e = jaccard_to_e(j, k);
            double e_lower = e_lower_bound(e, sketch_size, k, config::constants::confidence_interval);

            double upper_boud_identity = 100.0 * (1.0 - e_lower);

            if (upper_boud_identity >= config::constants::default_percentage_identity) {
                min_hits_relaxed = i;
            } else {
                break;
            }
        }
        return min_hits_relaxed;
    }

    double estimate_p_value(unsigned int sketch_size,
                            unsigned int k,
                            unsigned int q_length,
                            unsigned long r_length) {
        double possible_kmers = pow(config::constants::alphabet_size, k);
        double p_random_kmer = 1.0 - pow((1.0 - 1.0 / possible_kmers), q_length);
        double njihov = 1. / (1. + possible_kmers / q_length);
        double J = njihov * njihov / (njihov * (2 - njihov));

        unsigned int x = min_hits_relaxed_estimate(sketch_size, k);

        double ret = r_length * (x ? cdf(complement(binomial(sketch_size, J), x - 1)) : 1.0);
        return ret;


    }

    unsigned int calculate_window_size(unsigned int k,
                                       unsigned int q_length,
                                       unsigned long ref_length) {
        std::vector<unsigned int> sketch_values = {1, 2, 5};
        for (unsigned int i = 10; i < q_length; i += 10) {
            sketch_values.push_back(i);
        }

        double p_value;
        unsigned int optimal_sketch_size;
        for (unsigned int &sketch_value : sketch_values) {
            p_value = estimate_p_value(sketch_value, k, q_length, ref_length);
            if (p_value <= config::constants::p_value_cutoff) {
                optimal_sketch_size = sketch_value;
                break;
            }
        }
        unsigned int w = 2.0 * q_length / optimal_sketch_size;
        return std::min(std::max(1u, w), q_length);

    }
}