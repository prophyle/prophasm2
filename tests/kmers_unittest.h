#pragma once
#include "../src/kmers.h"

#include "gtest/gtest.h"

namespace {
    TEST(KMers, BitSuffix) {
        EXPECT_EQ(0b011110, BitSuffix(0b1100011110, 3));
        EXPECT_EQ(0b1110, BitSuffix(0b1110, 2));
        EXPECT_EQ(0b0, BitSuffix(0b1110, 0));
        EXPECT_EQ(0b111111'11111111'11111111'11111111'11111111'11111110LL, BitSuffix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL,  23));
        EXPECT_EQ(0b111111'11111110LL, BitSuffix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 7));
    }

    TEST(KMers, BitPrefix) {
        EXPECT_EQ(0b110001, BitPrefix(0b1100011110, 5, 3));
        EXPECT_EQ(0b1110, BitPrefix(0b1110, 2, 2));
        EXPECT_EQ(0b0, BitPrefix(0b1110, 2, 0));
        EXPECT_EQ(0b111111'01111111'11111111'11111111'11111111'11111111LL, BitPrefix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, 23));
        EXPECT_EQ(0b111111'01111111LL, BitPrefix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, 7));
    }

    TEST(KMers, NucleotideToInt) {
        struct TestCase {
            char nucleotide;
            int wantResult;
        };
        std::vector<TestCase> tests = {
                {'A', 0},
                {'C', 1},
                {'G', 2},
                {'T', 3},
                {'B', -1},
        };

        for (auto t : tests) {
            int gotResult = NucleotideToInt(t.nucleotide);
            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(KMers, NucleotideAtIndex) {
        EXPECT_EQ('G', NucleotideAtIndex(0b111001, 3, 1));
        EXPECT_EQ('G', NucleotideAtIndex(0b11100111, 4, 1));
        EXPECT_EQ('C', NucleotideAtIndex(0b11100001, 4, 3));
        EXPECT_EQ('T', NucleotideAtIndex(0b1101100001, 5, 0));
        EXPECT_EQ('T', NucleotideAtIndex(0b11, 1, 0));
    }

    TEST(KMers, NumberToKMer) {
        struct TestCase {
            kmer_t encoded;
            int d;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, 2, "GC"},
                {0b1011LL, 3, "AGT"},
                {0b111LL, 1, "T"},
                {0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, "TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"},
        };

        for (auto &&t: tests) {
            std::string gotResult = NumberToKMer(t.encoded, t.d);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(KMers, MaskForK) {
        struct TestCase {
            int k;
            kmer_t wantResult;
        };
        std::vector<TestCase> tests = {
                {2, 0b1111LL},
                {3, 0b111111LL},
                {1, 0b11LL},
                {32, 0b11111111'11111111'11111111'11111111'11111111'11111111'11111111'11111111LL },
        };

        for (auto t: tests) {
            kmer_t gotResult = MaskForK(t.k);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(KMers, ReverseComplement) {
        struct TestCase {
            kmer_t input;
            int k;
            kmer_t wantResult;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, 2, 0b1001LL},
                {0b101111LL, 3, 0b000001LL},
                {0b11LL, 1, 0b00LL},
                {0b11111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 32, 0b010000'00000000'00000000'00000000'00000000'00000000'00000000'1000000000LL },
        };

        for (auto t: tests) {
            kmer_t gotResult = ReverseComplement(t.input, t.k);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(KMers, CanonicalKMer) {
        struct TestCase {
            kmer_t input;
            int k;
            kmer_t wantResult;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, 2, 0b1001LL},
                {0b101111LL, 3, 0b000001LL},
                {0b11LL, 1, 0b00LL},
                {0b11111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 32, 0b010000'00000000'00000000'00000000'00000000'00000000'00000000'1000000000LL },
        };

        for (auto t: tests) {
            kmer_t gotResult = CanonicalKMer(t.input, t.k);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }
}
