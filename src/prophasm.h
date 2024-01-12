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
std::pair<kmer_t, kmer_t> RightExtension(kmer_t last, kh_S64_t *kMers, int k, bool complements) {
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
std::pair<kmer_t, kmer_t> LeftExtension(kmer_t first, kh_S64_t *kMers, int k,  bool complements) {
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
void NextSimplitig(kh_S64_t *kMers, kmer_t begin, std::ostream& of,  int k, bool complements, int simplitigID) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    std::string initialKMer = NumberToKMer(begin, k);
    std::list<char> simplitig {initialKMer.begin(), initialKMer.end()};
    eraseKMer(kMers, last, k, complements);
    // TODO: print the right simplitig part directly after constructing the left part.
    bool extendToRight = true;
    bool extendToLeft = true;
    while (extendToRight || extendToLeft) {
        if (extendToRight) {
            auto [ext, next] = RightExtension(last, kMers, k, complements);
            if (ext == kmer_t(-1)) {
                // No right extension found.
                extendToRight = false;
            } else {
                // Extend the simplitig to the right.
                eraseKMer(kMers, next, k, complements);
                simplitig.emplace_back(letters[ext]);
                last = next;
            }
        } else {
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
    }
    of << ">" << simplitigID << std::endl;
    of << std::string(simplitig.begin(), simplitig.end()) << std::endl;
}

/// Heuristically compute simplitigs.
///
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
/// Warning: this will destroy kMers.
void ComputeSimplitigs(kh_S64_t *kMers, std::ostream& of, int k, bool complements) {
    size_t lastIndex = 0;
    int simplitigID = 0;
    while(true) {
        kmer_t begin = nextKMer(kMers, lastIndex);
        // No more k-mers.
        if (begin == kmer_t(-1)) return;
        NextSimplitig(kMers, begin, of,  k, complements, simplitigID++);
    }
}
