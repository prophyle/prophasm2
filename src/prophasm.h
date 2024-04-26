#pragma once

#include <string>
#include <unordered_set>
#include <vector>
#include <deque>
#include <cstdint>
#include <algorithm>
#include <fstream>
#include <stack>

#include "kmers.h"
#include "khash_utils.h"


/// Find the right extension to the provided last k-mer from the kMers by trying to append each of {A, C, G, T}.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
template <typename KHT, typename kmer_t>
inline uint32_t RightExtension(kmer_t &last, kmer_t &complement, kmer_t &canonical, KHT *kMers, int k, bool complements) {
    kmer_t k2m1 = (k - 1) << 1;
    for (kmer_t ext = 0; ext < 4; ++ext) {
        kmer_t next = (BitSuffix(last, k - 1) << 2) | ext;
        kmer_t nextComplement = ((ext ^ 3) << k2m1) | (complement >> 2);
        canonical = next;
        if (complements) canonical = std::min(next, nextComplement);
        if (containsCanonicalKMer(kMers, canonical)) {
            last = next;
            complement = nextComplement;
            return ext;
        }
    }
    return -1;
}

/// Find the left extension to the provided first k-mer from the kMers by trying to prepend each of {A, C, G, T}.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
template <typename KHT, typename kmer_t>
inline uint32_t LeftExtension(kmer_t &first, kmer_t &complement, kmer_t &canonical, KHT *kMers, int k,  bool complements) {
    kmer_t k2m1 = (k - 1) << 1;
    for (kmer_t ext = 0; ext < 4; ++ext) {
        kmer_t next = (ext << k2m1) | (first >> 2);
        kmer_t nextComplement = (BitSuffix(complement, k - 1) << 2) | (ext ^ 3);
        canonical = next;
        if (complements) canonical = std::min(next, nextComplement);
        if (containsCanonicalKMer(kMers, canonical)) {
            first = next;
            complement = nextComplement;
            return ext;
        }
    }
    return -1;
}

/// Find the next simplitig.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
template <typename KHT, typename kmer_t>
void NextSimplitig(KHT *kMers, kmer_t begin, std::ostream& of,  int k, bool complements, int simplitigID) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    kmer_t complement = ReverseComplement(begin, k);
    kmer_t canonical;
    std::string initialKMer = NumberToKMer(begin, k);
    std::stack<char> simplitig {};
    for (int i = k - 1; i >= 0; --i) {
        simplitig.push(initialKMer[i]);
    }
    eraseKMer(kMers, last, k, complements);
    bool extendToRight = true;
    bool extendToLeft = true;
    while (extendToLeft) {
        uint32_t ext =  LeftExtension(first, complement,canonical, kMers, k, complements);
        if (ext == uint32_t(-1)) {
            // No left extension found.
            extendToLeft = false;
        } else {
            // Extend the simplitig to the left.
            eraseCanonicalKMer(kMers, canonical);
            simplitig.push(letters[ext]);
        }
    }
    complement = ReverseComplement(begin, k);
    of << ">" << simplitigID << std::endl;
    while (!simplitig.empty()) {
        of << simplitig.top();
        simplitig.pop();
    }
    while (extendToRight) {
        uint32_t ext = RightExtension(last, complement, canonical, kMers, k, complements);
        if (ext == uint32_t(-1)) {
            // No right extension found.
            extendToRight = false;
        } else {
            // Extend the simplitig to the right.
            eraseCanonicalKMer(kMers, canonical);
            of << letters[ext];
        }
    }
    of << std::endl;
}


#define INIT_PROPHASM(type, variant)                                                                                   \
                                                                                                                       \
/*  Heuristically compute simplitigs.                                                                                  \
 *                                                                                                                     \
 *  If complements are provided, treat k-mer and its complement as identical.                                          \
 *  If this is the case, k-mers are expected not to contain both k-mer and its complement.                             \
 *  Warning: this will destroy kMers.                                                                                  \
 */                                                                                                                    \
int ComputeSimplitigs(kh_S##variant##_t *kMers, std::ostream& of, int k, bool complements) {                           \
    size_t lastIndex = 0;                                                                                              \
    kmer##type##_t begin = 0;                                                                                          \
    int simplitigID = 0;                                                                                               \
    while(true) {                                                                                                      \
        bool found = nextKMer(kMers, lastIndex, begin);                                                                \
        /* No more k-mers. */                                                                                          \
        if (!found) return simplitigID;                                                                                \
        NextSimplitig(kMers, begin, of,  k, complements, simplitigID++);                                               \
    }                                                                                                                  \
}                                                                                                                      \
                                                                                                                       \
/* Data for parallel computation of simplitigs. */                                                                     \
struct ComputeSimplitigsData##variant {                                                                                \
    std::vector<kh_S##variant##_t *> kMers;                                                                            \
    std::vector<std::ostream*> ofs;                                                                                    \
    int k;                                                                                                             \
    bool complements;                                                                                                  \
    std::vector<int> simplitigsCounts;                                                                                 \
};                                                                                                                     \
                                                                                                                       \
/* Parallel wrapper for ComputeSimplitigs. */                                                                          \
void ComputeSimplitigsThread##variant(void *arg, long i, int _) {                                                      \
    auto *data = (ComputeSimplitigsData##variant *) arg;                                                               \
    data->simplitigsCounts[i] = ComputeSimplitigs(data->kMers[i], *data->ofs[i], data->k, data->complements);          \
}                                                                                                                      \


INIT_PROPHASM(64, 64S)
INIT_PROPHASM(64, 64M)
INIT_PROPHASM(128, 128M)
INIT_PROPHASM(256, 256M)
