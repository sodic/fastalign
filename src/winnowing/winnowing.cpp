//
// Created by filip on 04.04.18..
//
#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include "winnowing.h"

using namespace std;

namespace winnowing{

    //Thomas Wang's integer hash function
    uint64_t invertible__hash(uint64_t x) {
        uint64_t key = x;

        key = (~key) + (key << 21); // key = (key << 21) - key - 1;
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return key;
    }

    uint64_t invertible_hash_inverse(uint64_t key) {
        uint64_t tmp;

        // Invert key = key + (key << 31)
        tmp = key - (key << 31);
        key = key - (tmp << 31);

        // Invert key = key ^ (key >> 28)
        tmp = key ^ key >> 28;
        key = key ^ tmp >> 28;

        // Invert key *= 21
        key *= 14933078535860113213u;

        // Invert key = key ^ (key >> 14)
        tmp = key ^ key >> 14;
        tmp = key ^ tmp >> 14;
        tmp = key ^ tmp >> 14;
        key = key ^ tmp >> 14;

        // Invert key *= 265
        key *= 15244667743933553977u;

        // Invert key = key ^ (key >> 24)
        tmp = key ^ key >> 24;
        key = key ^ tmp >> 24;

        // Invert key = (~key) + (key << 21)
        tmp = ~key;
        tmp = ~(key - (tmp << 21));
        tmp = ~(key - (tmp << 21));
        key = ~(key - (tmp << 21));

        return key;
    }

    /**
     * Returns the associated hash value of a given nucleotide
     * @param c nucleotide representation
     * @return nucleotide hash value
     */
    size_t value_of_nucleotide(char c) {
        switch (c) {
            case 'A':
                return 0;

            case 'C':
                return 1;

            case 'G':
                return 2;

            case 'T':
                return 3;

            default:
                return 0;
        }
    }


    /**
     * Returns the complement base
     * @param c
     * @return
     */
    char complement(char c)
    {
        switch (c) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            default:
                throw runtime_error("Invalid value.");
        }
    }

    /**
     * Naive implementation of the minimizer hashing function. The kmer hash is calculated based on the position of each
     * nucleotide in the sequence as well as its hash value.
     * @param seq  pointer to the start of the nucleotide char sequence
     * @param seq_l length of the kmer
     * @return 64-bit hash value
     */
    uint64_t hash(const char *seq, uint32_t seq_l)
    {
        uint64_t hash = 0;

        for (int i = 0; i < seq_l; i++) {
            hash += value_of_nucleotide(seq[i]) * pow(4, seq_l - i - 1);
        }

        return hash;
    }

    /**
     * Naive implementation of the minimizer hashing function. The kmer hash is calculated based on the position of each
     * nucleotide in the reversed sequence as well as its hash value.
     * @param seq  pointer to the start of the nucleotide char sequence
     * @param seq_l length of the kmer
     * @return 64-bit hash value
     */
    uint64_t hash_rev(const char *seq, uint32_t seq_l) {
        uint64_t hash = 0;

        for (int i = seq_l - 1; i >= 0; i--) {
            hash += value_of_nucleotide(complement(seq[i])) * pow(4, i);
        }

        return hash;
    }


    /**
     * Calculates the minimizer hash of a kmer based on the value of the kmer before it and the last current nucleotide
     * @param seq sequence
     * @param index start index of the kmer in the sequence
     * @param last_hash pointer to last recorded kmer hash valie
     * @param power power to which values should be raised
     * @param first_nucleotide_value poitner to the first value of the nucleotide
     * @return new hash value
     */
    uint64_t based_hash(const char *seq, int32_t index, uint64_t *last_hash, uint32_t power,
                        uint64_t *first_nucleotide_value)
    {
        uint64_t  hash = 0;
        uint64_t remove = *first_nucleotide_value * (1 << (power*2));
        *first_nucleotide_value = value_of_nucleotide(seq[index]);
        hash = (*last_hash - remove) * 4 + value_of_nucleotide(seq[index + power]);
        *last_hash = hash;
        return hash;
    }

    /**
    * Calculates the reverse minimizer hash of a kmer based on the value of the kmer before it and the last current nucleotide
    * @param seq sequence
    * @param index start index of the kmer in the sequence
    * @param last_hash pointer to last recorded kmer hash valie
    * @param power power to which values should be raised
    * @param first_nucleotide_value poitner to the first value of the nucleotide
    * @return new hash value
    */
    uint64_t based_hash_rev(const char *seq, int32_t index, uint64_t *last_hash, uint32_t power,
                            uint64_t *first_nucleotide_value)
    {
        uint64_t  hash = 0;
        uint64_t remove = *first_nucleotide_value ;
        *first_nucleotide_value = value_of_nucleotide(complement(seq[index]));
        hash = (*last_hash - remove) / 4 + (value_of_nucleotide(complement(seq[index + power])) * (1 << (power*2)));
        *last_hash = hash;
        return hash;
    }

    /**
 * Find all minimizers for a given strand and stores them as a minim struct in a vector.
 * @param seq sequence
 * @param seq_l sequence length
 * @param w window size
 * @param k kmer length
 * @param minimizers  reference to vector where all minimizers will be stored
 */
    void compute_minimizers
            (const char *seq,
             uint32_t seq_l,
             int32_t w,
             uint32_t k,
             std::vector <minimizer> &minimizers
            )
    {

        if(seq_l < w + k){
            return;
        }

        const uint32_t kmers_l = seq_l - k + 1;

        uint64_t *hash_buffer = new uint64_t[w];
        uint64_t *r_hash_buffer = new uint64_t[w];
        uint64_t first_nuc_val = value_of_nucleotide(seq[0]);
        uint64_t first_nuc_val_r = value_of_nucleotide(complement(seq[0]));
        uint32_t power = k - 1;
        uint64_t prev_hash = hash(seq, k);
        uint64_t prev_hash_r = hash_rev(seq, k);

        hash_buffer[0] = invertible__hash(prev_hash);
        r_hash_buffer[0] = invertible__hash(prev_hash_r);

        for (uint32_t i = 1; i < w; i++) {
            hash_buffer[i] = invertible__hash(based_hash(seq, i, &prev_hash, power, &first_nuc_val));
            r_hash_buffer[i] = invertible__hash(based_hash_rev(seq, i, &prev_hash_r, power, &first_nuc_val_r));   //HASH
        }

        uint32_t min_l_pred = seq_l - w - k + 2;
        uint64_t last_min_hash = UINT64_MAX;
        int64_t last_min_position = -1;

        for (uint32_t i = 0; i < min_l_pred; i++) {
            uint64_t u;
            uint64_t v;

            if(last_min_position != -1 && abs(last_min_position) >= i){
                u = hash_buffer[(i + w - 1) % w];
                v = r_hash_buffer[(i + w - 1) % w];

                if(u == v){
                    continue;
                }

                else if(u < v && u <= last_min_hash){
                    minimizers.emplace_back((minimizer) {u, (i + w - 1), 1});
                    last_min_position = i + w - 1;
                    last_min_hash = u;
                }

                else if(u > v && v <= last_min_hash) {
                    minimizers.emplace_back((minimizer) {v, (i + w - 1), -1});
                    last_min_position = (i + w - 1);
                    last_min_hash = v;
                }
            }
            else {
                uint64_t m = UINT64_MAX;

                uint32_t *min_positions = new uint32_t[w];
                uint16_t min_pos_size = 0;
                for (int j = 0; j < w; j++) {
                    u = hash_buffer[(i + j) % w];
                    v = r_hash_buffer[(i + j) % w];

                    if(u == v){
                        continue;
                    }

                    if(u < m || v <  m){
                        if(u < v){
                            min_positions[0] = i + j;
                            m = u;

                        } else {
                            min_positions[0] = (i + j);
                            m = v;
                        }

                        min_pos_size = 1;
                    }

                    else if(u == m){
                        min_positions[min_pos_size++] = i + j;
                    }

                    else if(v == m){
                        min_positions[min_pos_size++] = (i + j);
                    }
                }

                last_min_hash = m;
                last_min_position = min_positions[min_pos_size - 1];

                for(uint32_t j = 0; j < min_pos_size; j++) {
                    minimizers.emplace_back((minimizer) {m,  min_positions[j]}); //TREBA NEGATIVNO
                }


                delete[] min_positions;
            }

            int next_end = i + w;
            if (next_end < kmers_l) {
                hash_buffer[next_end % w] = invertible__hash(
                        based_hash(seq, next_end, &prev_hash, power, &first_nuc_val));//HASH
                r_hash_buffer[next_end % w] = invertible__hash(
                        based_hash_rev(seq, next_end, &prev_hash_r, power, &first_nuc_val_r));   //HASH
            }
        }

        delete[] hash_buffer;
        delete[] r_hash_buffer;
    }

    void compute_minimizers(const char *seq,
                            uint32_t seq_l,
                            std::vector<minimizer> &minimizers) {
        compute_minimizers(seq, seq_l, DEFAULT_W, DEFAULT_K, minimizers);
    }

    /**
 * Wrapper function for finding all minimizers of a sequence and embedding the resulting vector into a vector of vectors
 * which keeps track of all minimizer vectors for all sequnces.
 * @param sequence nucleotide sequence from which minimizers should be extracted
 * @param sequence_l sequence length
 * @param w window size
 * @param k kmers size
 * @param the vector where the minimizers should be stored
 */
    void index_sequence(const char* sequence,
                        uint32_t sequence_l,
                        uint32_t w,
                        uint32_t k,
                        std::vector<minimizer>& minimizers,
                        unordered_map<minhash_t, vector<uint32_t >> &lookup_table) {

        compute_minimizers(sequence, sequence_l, w, k, minimizers);
        for (auto &minimizer : minimizers) {
            lookup_table[minimizer.hash].push_back(minimizer.index);
        }
        minimizers.shrink_to_fit();
    }

    void index_sequence(const char *sequence,
                        uint32_t sequence_l,
                        std::vector<winnowing::minimizer> &minimizers,
                        std::unordered_map<minhash_t, std::vector<std::uint32_t >> &lookup_table) {
        index_sequence(sequence, sequence_l, DEFAULT_W, DEFAULT_K, minimizers, lookup_table);
    }


}