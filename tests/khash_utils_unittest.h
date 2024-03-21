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
            auto kMers = kh_init_S64M();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64M(kMers, kMer, &ret);

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
                {{{0b0101}, {0b0101}}, 2, true, {0b0101}},
        };

        for (auto t : tests) {
            std::vector<kh_S64M_t *> input (t.kMerSets.size());
            for (size_t i = 0; i < t.kMerSets.size(); ++i) {
                input[i] = kh_init_S64M();
                int ret;
                for (auto &&kMer : t.kMerSets[i]) kh_put_S64M(input[i], kMer, &ret);
            }

            kh_S64M_t * gotResult = kh_init_S64M();
            getIntersection(gotResult, input, t.k, t.complements);
            auto gotResultVec = kMersToVec(gotResult);
            sort(t.wantResult.begin(), t.wantResult.end());
            sort(gotResultVec.begin(), gotResultVec.end());

            EXPECT_EQ(t.wantResult, gotResultVec);
        }
    }

    TEST(KHASH_UTILS, DifferenceInPlace) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            std::vector<kmer_t> toRemove;
            int k;
            bool complements;
            std::vector<kmer_t> wantResult;
        };

        std::vector<TestCase> tests = {
                // {TCC, CTA, ACT, CCT}
                {{0b110101, 0b011100, 0b000111, 0b010111},
                 {0b110101}, 3, false,
                 {0b011100, 0b000111, 0b010111}},
                // {TCC, CTA, CCT}
                {{0b110101, 0b011100, 0b010111},
                 {0b110101, 0b000111}, 3, false,
                 {0b011100, 0b010111},},
                // {GG}
                { {0b0101}, {0b0101}, 2, true, {}},
        };

        for (auto t : tests) {
            kh_S64M_t * input = kh_init_S64M();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64M(input, kMer, &ret);
            auto intersection = kh_init_S64M();
            for (auto &&kMer : t.toRemove) kh_put_S64M(intersection, kMer, &ret);

            differenceInPlace(input, intersection, t.k, t.complements);

            std::vector<kmer_t> gotResult = kMersToVec(input);
            sort(t.wantResult.begin(), t.wantResult.end());
            sort(gotResult.begin(), gotResult.end());
            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(KHASH_UTILS, AbundanceIntegration) {
        struct TestCase {
            std::vector<std::pair<byte, kmer_t>> kMers;
            int k;
            bool complements;
            std::vector<std::pair<byte, std::vector<kmer_t>>> wantResults;
        };

        std::vector<TestCase> tests = {
                {
                    // 1x TAT, 2x CCC, 4x TTT, 1x AAA
                    {{1, 0b110011}, {2, 0b010101}, {4, 0b111111}, {1, 0b000000}},
                    3, false,
                    {{4, {0b111111}}, {5, {}}, {2, {0b010101, 0b111111}}}
                },
        };

        for (auto t : tests) {
            // Turn off optimizations for MINIMUM_ABUNDANCE being 1.
            MINIMUM_ABUNDANCE = 255;
            kh_S64M_t * input = kh_init_S64M();
            for (auto &&[abundance, kMer] : t.kMers) for (byte i = 0; i < abundance; ++i)
                insertKMer(input, kMer, t.k, t.complements);

            for (auto &&[abundance, wantResult] : t.wantResults) {
                MINIMUM_ABUNDANCE = abundance;
                for (auto &&[_, kMer] : t.kMers) {
                    bool contains = containsKMer(input, kMer, t.k, t.complements);
                    bool should_contain = wantResult.end() != find(wantResult.begin(), wantResult.end(), kMer);
                    EXPECT_EQ(should_contain, contains);
                }
            }
        }
    }
}
