#pragma once

#include <string>
#include <unordered_set>
#include <vector>
#include <deque>
#include <cstdint>
#include <algorithm>
#include <fstream>

#include "kmers.h"
#include "khash_utils.h"


/// Find the right extension to the provided last k-mer from the kMers by trying to append each of {A, C, G, T}.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
template <typename KHT, typename kmer_t>
std::pair<kmer_t, kmer_t> RightExtension(kmer_t last, KHT *kMers, int k, bool complements) {
    for (kmer_t ext = 0; ext < 4; ++ext) {
        kmer_t next = (BitSuffix(last, k - 1) << 2) | ext;
        if (containsKMer(kMers, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the left extension to the provided first k-mer from the kMers by trying to prepend each of {A, C, G, T}.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
template <typename KHT, typename kmer_t>
std::pair<kmer_t, kmer_t> LeftExtension(kmer_t first, KHT *kMers, int k,  bool complements) {
    for (kmer_t ext = 0; ext < 4; ++ext) {
        kmer_t next = (ext << ((k - 1) << 1)) | BitPrefix(first, k, k - 1);
        if (containsKMer(kMers, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the next simplitig.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
template <typename KHT, typename kmer_t>
void NextSimplitig(KHT *kMers, kmer_t begin, std::ostream& of,  int k, bool complements, int simplitigID) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    std::string initialKMer = NumberToKMer(begin, k);
    std::list<char> simplitig {initialKMer.begin(), initialKMer.end()};
    eraseKMer(kMers, last, k, complements);
    bool extendToRight = true;
    bool extendToLeft = true;
    while (extendToLeft) {
        auto [ext, next] =  LeftExtension(first, kMers, k, complements);
        if (ext == kmer_t(-1)) {
            // No left extension found.
            extendToLeft = false;
        } else {
            // Extend the simplitig to the left.
            eraseKMer(kMers, next, k, complements);
            simplitig.emplace_front(letters[ext]);
            first = next;
        }
    }
    of << ">" << simplitigID << std::endl;
    of << std::string(simplitig.begin(), simplitig.end());
    simplitig.resize(0);
    while (extendToRight) {
        auto [ext, next] = RightExtension(last, kMers, k, complements);
        if (ext == kmer_t(-1)) {
            // No right extension found.
            extendToRight = false;
        } else {
            // Extend the simplitig to the right.
            eraseKMer(kMers, next, k, complements);
            of << letters[ext];
            last = next;
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
INIT_PROPHASM(128, 128S)
INIT_PROPHASM(128, 128M)
