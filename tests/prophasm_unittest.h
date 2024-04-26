#pragma once
#include "../src/prophasm.h"

#include "gtest/gtest.h"

namespace {
    TEST(Prophasm, RightExtension) {
        struct TestCase {
            kmer_t last;
            kmer_t complement; // Irrelevant if complements is set to false.
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            uint32_t wantExt;
            kmer_t wantNext;
            kmer_t wantComplement;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, CTA, ACT, CCT}; A; CTA
                {0b000111, 0b001011, {0b110101, 0b011100, 0b000111, 0b010111}, 3, false, 0b00, 0b011100, 0b110010},
                // ACT; {TCC, ACT, CCT}
                {0b000111, 0b001011, {0b110101, 0b000111, 0b010111}, 3, false,  uint32_t(-1), 0b000111, 0b001011},
        };

        for (auto t: tests) {
            auto kMers = kh_init_S64M();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64M(kMers, kMer, &ret);

            auto got = RightExtension(t.last, t.complement, kMers, t.k, t.complements);

            EXPECT_EQ(t.wantNext, t.last);
            EXPECT_EQ(t.wantComplement, t.complement);
            EXPECT_EQ(t.wantExt, got);
        }
    }

    TEST(Prophasm, LeftExtension) {
        struct TestCase {
            kmer_t first;
            kmer_t complement; // Irrelevant if complements is set to false.
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            uint32_t wantExt;
            kmer_t wantNext;
            kmer_t wantComplement;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, ACT, CCT}
                {0b000111, 0b001011, {0b110101, 0b000111, 0b010111}, 3, false,  uint32_t(-1), 0b000111, 0b001011},
                // TAC; {TCC, CTA, ACT, CCT}; C; CTA
                {0b110001, 0b101100, {0b110101, 0b011100, 0b000111, 0b010111}, 3, false, 0b01, 0b011100, 0b110010},
        };

        for (auto t: tests) {
            auto kMers = kh_init_S64M();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64M(kMers, kMer, &ret);

            auto got = LeftExtension(t.first, t.complement, kMers, t.k, t.complements);
            EXPECT_EQ(t.wantNext, t.first);
            EXPECT_EQ(t.wantComplement, t.complement);
            EXPECT_EQ(t.wantExt, got);
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
                // {ACAA, AAAT, AACA} complements: {TTGT, ATTT, TGTT}
                {{0b00010000, 0b00000011, 0b00000100}, 4, true, 42,
                 ">42\nAACAA\n", {0b00000011},},
        };

        for (auto &&t: tests) {
            std::stringstream of;
            auto kMers = kh_init_S64M();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64M(kMers, kMer, &ret);

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
                // {TTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG, TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA}
                {{0b1111111101111111111111111111111111111111111111111111111111111110,
                  0b1111110111111111111111111111111111111111111111111111111111111000},
                 32, false, ">0\nTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA\n",},
        };

        for (auto t: tests) {
            std::stringstream of;
            auto kMers = kh_init_S64M();
            int ret;
            for (auto &&kMer : t.kMers) kh_put_S64M(kMers, kMer, &ret);

            ComputeSimplitigs(kMers, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
