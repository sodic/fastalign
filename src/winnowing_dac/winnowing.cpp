//
// Created by dario on 02.11.17..
//

#include <queue>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "winnowing.hpp"

namespace {
    std::unordered_map<winnowing::minhash_t, int> hashCnt;
    std::map<char, int> baseValue = {{'A', 3},
                                     {'T', 2},
                                     {'G', 1},
                                     {'C', 0}};
    const int BASE = 4;
    std::map<char, char> rcMap = {{'A', 'T'},
                                  {'C', 'G'},
                                  {'T', 'A'},
                                  {'G', 'C'}};

    void push(const winnowing::minimizer &triple, std::deque<winnowing::minimizer> &dq) {
        while (!dq.empty() && triple < dq.back()) {
            dq.pop_back();
        }
        dq.push_back(triple);
    }

    void processState(std::deque<winnowing::minimizer> &dq,
                      std::vector<winnowing::minimizer> &result, int &lastPositionTaken) {
        assert(!dq.empty());
        winnowing::minimizer &front = dq.front();
        dq.pop_front();

        if (lastPositionTaken < front.index) {
            result.push_back(front);
            hashCnt[front.hash]++;
            lastPositionTaken = front.index;
        }
        while (!dq.empty() && dq.front().hash == front.hash) {
            front = dq.front();
            dq.pop_front();
            if (lastPositionTaken < front.index) {
                result.push_back(front);
                hashCnt[front.hash]++;
                lastPositionTaken = front.index;
            }
        }
        dq.push_front(front);
    }

    void pop(int position, std::deque<winnowing::minimizer> &dq) {
        while (!dq.empty() && dq.front().index == position)
            dq.pop_front();
    }


} // namespace

namespace winnowing {

    void compute_minimizers(const char *target,
                            const uint32_t target_length,
                            int w,
                            int k,
                            std::vector<minimizer> &minimizers) {

        if (target_length < k) {
            return; // ne postoji ni jedan kmer
        }

        if (target_length < k + w - 1) {//ne postoji ni jedan window od w kmera
            w = target_length - k +
                1; // smanji velicinu trazenog prozora na najvise sta mozes, da se nadje barem jedan winnowing
        }


        std::deque<minimizer> dq;

        minhash_t lastPower = 1;
        for (int i = 0; i < k - 1; i++)
            lastPower *= BASE;

        minhash_t tmpHash = 0;
        minhash_t tmpRcHash = 0;

        minhash_t tmpPot = 1;
        for (int i = 0; i < k; i++) {
            tmpHash *= BASE;
            tmpHash += baseValue[target[i]];
            tmpRcHash += tmpPot * baseValue[rcMap[target[i]]];
            tmpPot *= BASE;
        }

        // queue s maksimumom algoritam
        for (int i = 0; i < w; i++) {
            minimizer mp1 = minimizer(tmpHash, i, 1);
            minimizer mp2 = minimizer(tmpRcHash, i, -1);
            push(mp1, dq);
            push(mp2, dq);
            tmpHash -= lastPower * baseValue[target[i]];
            tmpHash *= BASE;
            tmpHash += baseValue[target[i + k]];

            tmpRcHash -= baseValue[rcMap[target[i]]];
            tmpRcHash /= BASE;
            tmpRcHash += lastPower * baseValue[rcMap[target[i + k]]];
        }


        int lastPositionTaken = -1;

        processState(dq, minimizers, lastPositionTaken);

        for (int i = w; i < target_length - k + 1; i++) {
            pop(i - w, dq);
            minimizer mp1 = minimizer(tmpHash, i, 1);
            minimizer mp2 = minimizer(tmpRcHash, i, -1);
            push(mp1, dq);
            push(mp2, dq);
            processState(dq, minimizers, lastPositionTaken);

            tmpHash -= lastPower * baseValue[target[i]];
            tmpHash *= BASE;
            tmpHash += baseValue[target[i + k]];

            tmpRcHash -= baseValue[rcMap[target[i]]];
            tmpRcHash /= BASE;
            tmpRcHash += lastPower * baseValue[rcMap[target[i + k]]];
        }

    }


    void index_sequence(const char *sequence,
                        const uint32_t seq_length,
                        uint32_t w,
                        uint32_t k,
                        std::vector<minimizer> &minimizers,
                        std::unordered_map<minhash_t, std::vector<uint32_t >> &lookup_table) {

        compute_minimizers(sequence, seq_length, w, k, minimizers);
        for (auto &minimizer : minimizers) {
            lookup_table[minimizer.hash].push_back(minimizer.index);
        }
        minimizers.shrink_to_fit();
    }


} // namespace winnowing
