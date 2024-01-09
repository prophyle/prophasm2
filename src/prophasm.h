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


/// Find the right extension to the provided last k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<kmer_t, kmer_t> RightExtension(kmer_t last, kh_S64_t *kMers, int k, int d, bool complements) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    // TODO: remove kmercamel relicts.
    for (kmer_t ext = 0; ext < (kmer_t(1) << (d << 1)); ++ext) {
        kmer_t next = BitSuffix(last, k - d) << (d << 1) | ext;
        if (containsKMer(kMers, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the left extension to the provided first k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<kmer_t, kmer_t> LeftExtension(kmer_t first, kh_S64_t *kMers, int k, int d, bool complements) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    // TODO: remove kmercamel relicts.
    for (kmer_t ext = 0; ext < (kmer_t(1) << (d << 1)); ++ext) {
        kmer_t next = ext << ((k - d) << 1) | BitPrefix(first, k, k - d);
        if (containsKMer(kMers, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the next simplitig.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
void NextSimplitig(kh_S64_t *kMers, kmer_t begin, std::ostream& of,  int k, bool complements, int simplitig_id) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    std::list<char> simplitig {NucleotideAtIndex(first, k, 0)};
    eraseKMer(kMers, last, k, complements);
    // TODO: remove kmercamel relicts.
    int d_l = 1, d_r = 1;
    while (d_l <= 1 || d_r <= 1) {
        if (d_r <= d_l) {
            auto extension = RightExtension(last, kMers, k, 1, complements);
            kmer_t ext = extension.first;
            if (ext == kmer_t(-1)) {
                // No right extension found.
                ++d_r;
            } else {
                // Extend the generalized simplitig to the right.
                eraseKMer(kMers, extension.second, k, complements);
                for (int i = 0; i < d_r; ++i) simplitig.emplace_back(NucleotideAtIndex(extension.second, k, i));
                last = extension.second;
                d_r = 1;
            }
        } else {
            auto extension = LeftExtension(first, kMers, k, 1, complements);
            kmer_t ext = extension.first;
            if (ext == kmer_t(-1)) {
                // No left extension found.
                ++d_l;
            } else {
                // Extend the simplitig to the left.
                eraseKMer(kMers, extension.second, k, complements);
                for (int i = d_l - 1; i >= 0; --i) simplitig.emplace_front(NucleotideAtIndex(extension.second, k, i));
                first = extension.second;
                d_l = 1;
            }
        }
    }
    for (int i = 1; i < k; ++i) simplitig.emplace_back(NucleotideAtIndex(last, k, i));
    of << ">" << simplitig_id << std::endl;
    of << std::string(simplitig.begin(), simplitig.end()) << std::endl;
}

/// Heuristically compute simplitigs.
///
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
/// Warning: this will destroy kMers.
void ComputeSimplitigs(kh_S64_t *kMers, std::ostream& of, int k, bool complements) {
    size_t lastIndex = 0;
    int simplitig_id = 0;
    while(true) {
        kmer_t begin = nextKMer(kMers, lastIndex);
        // No more k-mers.
        if (begin == kmer_t(-1)) return;
        NextSimplitig(kMers, begin, of,  k, complements, simplitig_id++);
    }
}
