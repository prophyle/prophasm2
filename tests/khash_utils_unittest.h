#pragma once
#include "../src/khash_utils.h"

#include "gtest/gtest.h"

namespace {
    TEST(KHASH_UTILS, KMersToVec) {
        struct TestCase {
            std::vector<kmer_t> kMers;
        };
        std::vector<TestCase> tests = {
                // {TCC, CTA, ACT, CCT}
                {{0b110101, 0b011100, 0b000111, 0b010111}},
                // {TCC, ACT, CCT}
                {{0b110101, 0b000111, 0b010111}},
        };

        for (auto t: tests) {
            auto kMers = kh_init_S64();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64(kMers, kMer, &ret);

            auto got = kMersToVec(kMers);
            sort(t.kMers.begin(), t.kMers.end());
            sort(got.begin(), got.end());

            EXPECT_EQ(t.kMers, got);
        }
    }

    TEST(KHASH_UTILS, Intersection) {
        struct TestCase {
            std::vector<std::vector<kmer_t>> kMerSets;
            int k;
            bool complements;
            std::vector<kmer_t> wantResult;
        };
        std::vector<TestCase> tests = {
                // {{TCC, CTA, ACT, CCT}, {TCC, ACT, CCT}}
                {{{0b110101, 0b011100, 0b000111, 0b010111}, {0b110101, 0b000111, 0b010111}}, 3, false, {0b110101, 0b000111, 0b010111}},
                // {{TCC, CTA, ACT, CCT}, {TCC, ACT, CCT}, {TCC, CTA, ACT, TCT}}
                {{{0b110101, 0b011100, 0b000111, 0b010111}, {0b110101, 0b000111, 0b010111}, {0b110101, 0b011100, 0b000111, 0b110111}}, 3, false, {0b110101, 0b000111}},
                // {{CC}, {GG}}
                {{{0b0101}, {0b1010}}, 2, true, {0b0101}},
        };

        for (auto t : tests) {
            std::vector<kh_S64_t *> input (t.kMerSets.size());
            for (size_t i = 0; i < t.kMerSets.size(); ++i) {
                input[i] = kh_init_S64();
                int ret;
                for (auto &&kMer : t.kMerSets[i]) kh_put_S64(input[i], kMer, &ret);
            }

            auto gotResult = getIntersection(input, t.k, t.complements);
            auto gotResultVec = kMersToVec(gotResult);
            sort(t.wantResult.begin(), t.wantResult.end());
            sort(gotResultVec.begin(), gotResultVec.end());

            EXPECT_EQ(t.wantResult, gotResultVec);
        }
    }

    TEST(KHASH_UTILS, DifferenceInPlace) {
        struct TestCase {
            std::vector<std::vector<kmer_t>> kMers;
            std::vector<kmer_t> toRemove;
            int k;
            bool complements;
            std::vector<std::vector<kmer_t>> wantResult;
        };

        std::vector<TestCase> tests = {
                // {{TCC, CTA, ACT, CCT}, {TCC, ACT, CCT}}
                {{{0b110101, 0b011100, 0b000111, 0b010111}, {0b110101, 0b000111, 0b010111}},
                 {0b110101}, 3, false,
                 {{0b011100, 0b000111, 0b010111}, {0b000111, 0b010111}}},
                // {{TCC, CTA, ACT, CCT}, {TCC, CTA, CCT}, {TCC, CTA, ACT, TCT}}
                {{{0b110101, 0b011100, 0b000111, 0b010111},
                  {0b110101, 0b011100, 0b010111}, {0b110101, 0b011100, 0b000111, 0b110111}},
                 {0b110101, 0b000111}, 3, false,
                 {{0b011100, 0b010111}, {0b011100, 0b010111}, {0b110111, 0b011100}}},
                // {{CC}, {GG}}
                {{{0b0101}, {0b1010}}, {0b0101}, 2, true, {{}, {}}},
        };

        for (auto t : tests) {
            std::vector<kh_S64_t *> input (t.kMers.size());
            for (size_t i = 0; i < t.kMers.size(); ++i) {
                input[i] = kh_init_S64();
                int ret;
                for (auto &&kMer : t.kMers[i]) kh_put_S64(input[i], kMer, &ret);
            }
            auto intersection = kh_init_S64();
            int ret;
            for (auto &&kMer : t.toRemove) kh_put_S64(intersection, kMer, &ret);

            differenceInPlace(input, intersection, t.k, t.complements);

            std::vector<std::vector<kmer_t>> gotResult (t.kMers.size());
            for (size_t i = 0; i < t.kMers.size(); ++i) {
                gotResult[i] = kMersToVec(input[i]);
                sort(t.wantResult[i].begin(), t.wantResult[i].end());
                sort(gotResult[i].begin(), gotResult[i].end());
                EXPECT_EQ(t.wantResult[i], gotResult[i]);
            }
        }
    }
}
