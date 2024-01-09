#pragma once
#include "../src/prophasm.h"

#include "gtest/gtest.h"

namespace {
    TEST(Prophasm, RightExtension) {
        struct TestCase {
            kmer_t last;
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            kmer_t wantExt;
            kmer_t  wantNext;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, CTA, ACT, CCT}; A; CTA
                {0b000111, {0b110101, 0b011100, 0b000111, 0b010111}, 3, false, 0b00, 0b011100},
                // ACT; {TCC, ACT, CCT}
                {0b000111, {0b110101, 0b000111, 0b010111}, 3, false,  kmer_t(-1), kmer_t(-1)},
        };

        for (auto t: tests) {
            auto kMers = kh_init_S64();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64(kMers, kMer, &ret);

            auto got = RightExtension(t.last, kMers, t.k, 1, t.complements);
            auto gotExt = got.first;
            auto gotNext = got.second;

            EXPECT_EQ(t.wantNext, gotNext);
            EXPECT_EQ(t.wantExt, gotExt);
        }
    }

    TEST(Prophasm, LeftExtension) {
        struct TestCase {
            kmer_t first;
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            kmer_t wantExt;
            kmer_t  wantNext;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, ACT, CCT}
                {0b000111, {0b110101, 0b000111, 0b010111}, 3, false,   kmer_t(-1), kmer_t(-1)},
                // TAC; {TCC, CTA, ACT, CCT}; C; CTA
                {0b110001, {0b110101, 0b011100, 0b000111, 0b010111}, 3, false, 0b01, 0b011100},
        };

        for (auto t: tests) {
            auto kMers = kh_init_S64();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64(kMers, kMer, &ret);

            auto got = LeftExtension(t.first, kMers, t.k, 1, t.complements);
            auto gotExt = got.first;
            auto gotNext = got.second;
            EXPECT_EQ(t.wantNext, gotNext);
            EXPECT_EQ(t.wantExt, gotExt);
        }
    }


    TEST(Prophasm, NextSimplitig) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            int simplitigID;
            std::string wantResult;
            // Sorted k-mers.
            std::vector<kmer_t> wantKMers;
        };
        std::vector<TestCase> tests = {
                // {ACAA, AACA}
                {{0b00010000,  0b00000100}, 4, false, 0,
                 ">0\nAACAA\n", {},},
                // {ACAA, ATTT, TGTT} complements: {AAAT, TTGT, AACA}
                {{0b00010000, 0b00111111, 0b11101111}, 4, true, 42,
                 ">42\nAACAA\n", {0b00111111},},
        };

        for (auto &&t: tests) {
            std::stringstream of;
            auto kMers = kh_init_S64();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64(kMers, kMer, &ret);

            NextSimplitig(kMers, t.kMers.front(), of, t.k, t.complements, t.simplitigID);
            auto remainingKmers = kMersToVec(kMers);
            std::sort(remainingKmers.begin(), remainingKmers.end());

            EXPECT_EQ(t.wantKMers, remainingKmers);
            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(Prophasm, Prophasm) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                // {GCT, TAA, AAA}
                {{0b100111, 0b110000, 0b000000}, 3, false, ">0\nTAAA\n>1\nGCT\n",},
                // {TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG, TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA}
                {{0b11111101111111111111111111111111111111111111111111111111111110,
                  0b11110111111111111111111111111111111111111111111111111111111000},
                 31, false, ">0\nTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA\n",},
        };

        for (auto t: tests) {
            std::stringstream of;
            auto kMers = kh_init_S64();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64(kMers, kMer, &ret);

            ComputeSimplitigs(kMers, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
