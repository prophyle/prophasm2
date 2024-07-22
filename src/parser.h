#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <vector>

#include "kmers.h"
#include "khash_utils.h"

// Conversion table for nucleotides to integers.
// Contains 0 for A and a, 1 for C and c, 2 for G and g, 3 for T and t.
// Non-nucleotide characters have value >= 4. In particular, white space characters have value 5 and other 4.
static const int nucleotideToInt[] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 5, 5,  5, 5, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        5, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


#define INIT_PARSER(type, variant)                                                                        \
/*  Read encoded k-mers from the given fasta file.                                                        \
 *  Return unique k-mers in no particular order.                                                          \
 *  If complements is set to true, the result contains only one of the complementary k-mers               \
 *  - it is not guaranteed which one.                                                                     \
 *  This runs in O(sequence length) expected time.                                                        \
 */                                                                                                       \
void ReadKMers(kh_S##variant##_t *kMers, std::string &path, int k, bool complements) {                    \
    std::ifstream filestream;                                                                             \
    std::istream *fasta;                                                                                  \
    if (path == "-") {                                                                                    \
        fasta = &std::cin;                                                                                \
    } else {                                                                                              \
        filestream.open(path);                                                                            \
        fasta = &filestream;                                                                              \
        if (!filestream.is_open()) {                                                                      \
            std::cerr << "Error: file '" << path << "' could not be open." << std::endl;                  \
            exit(1);                                                                                      \
        }                                                                                                 \
    }                                                                                                     \
    unsigned char c;                                                                                      \
    int beforeKMerEnd = k;                                                                                \
    kmer##type##_t currentKMer = 0;                                                                       \
    kmer##type##_t complement = 0;                                                                        \
    /* mask that works even for k=32. */                                                                  \
    int k2Minus1 = 2 * k - 1;                                                                             \
    int k2Minus2 = 2 * k - 2;                                                                             \
    kmer##type##_t mask = (((kmer##type##_t) 1) <<  k2Minus1);                                            \
    mask |= mask - 1;                                                                                     \
    bool readingHeader = false;                                                                           \
    while ((*fasta) >> std::noskipws >> c) {                                                              \
        int data = nucleotideToInt[c];                                                                    \
        if (data >= 4 || readingHeader) {                                                                 \
            readingHeader |= (c == '>');                                                                  \
            readingHeader &= (c != '\n');                                                                 \
            if (data == 5) continue;                                                                      \
            currentKMer = 0;                                                                              \
            beforeKMerEnd = k;                                                                            \
            continue;                                                                                     \
        }                                                                                                 \
        currentKMer <<= 2;                                                                                \
        currentKMer &= mask;                                                                              \
        currentKMer |= data;                                                                              \
        complement >>= 2;                                                                                 \
        complement |= (kmer##type##_t((3) ^ data)) << k2Minus2;                                           \
        beforeKMerEnd -= (beforeKMerEnd > 0);                                                             \
        if (beforeKMerEnd == 0) {                                                                         \
            kmer##type##_t canonicalKMer = ((!complements) || currentKMer < complement) ?                 \
                    currentKMer : complement;                                                             \
            insertCanonicalKMer(kMers, canonicalKMer);                                                    \
        }                                                                                                 \
    }                                                                                                     \
    if (filestream.is_open()) filestream.close();                                                         \
}                                                                                                         \
                                                                                                          \
                                                                                                          \
/* Data for parallel reading of k-mers. */                                                                \
struct ReadKMersData##variant {                                                                           \
    std::vector<kh_S##variant##_t*> kMers;                                                                \
    std::vector<std::string> paths;                                                                       \
    int k;                                                                                                \
    bool complements;                                                                                     \
};                                                                                                        \
                                                                                                          \
/* Parallel wrapper for ReadKMers. */                                                                     \
void ReadKMersThread##variant(void *arg, long i, int _) {                                                 \
    auto *data = (ReadKMersData##variant *) arg;                                                          \
    ReadKMers(data->kMers[i], data->paths[i], data->k, data->complements);                                \
}                                                                                                         \

INIT_PARSER(64, 64S)
INIT_PARSER(64, 64M)
INIT_PARSER(128, 128M)
INIT_PARSER(256, 256M)
