//
// Created by filip on 29.05.18..
//

#ifndef ZAVRSNI_STATISTICS_H_H
#define ZAVRSNI_STATISTICS_H_H



namespace statistics {

    double jaccard_to_e(double jaccard, int k);

    double e_to_jaccard(float e, int k);

    double e_lower_bound(double e,
                         unsigned int sketch_size,
                         unsigned int k,
                         double confidence_interval);

    unsigned int minimum_hits_estimate(unsigned int sketch_size,
                                       unsigned int k);

    unsigned int min_hits_relaxed_estimate(unsigned int sketch_size,
                                           unsigned int k);

    double estimate_p_value(unsigned int sketch_size,
                            unsigned int k,
                            unsigned int q_length,
                            unsigned long r_length);

    unsigned int calculate_window_size(unsigned int k,
                                       unsigned int q_length,
                                       unsigned long ref_length);
}

#endif //ZAVRSNI_STATISTICS_H_H
