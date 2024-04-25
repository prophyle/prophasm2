#pragma once

#include <string>
#include <iostream>
#include <cstdint>

#include "uint256_t/uint256_t.h"

typedef uint64_t kmer64_t;
typedef __uint128_t kmer128_t;
typedef uint256_t kmer256_t;

/// Convert the given basic nucleotide to int so it can be used for indexing in AC.
/// If non-existing nucleotide is given, return -1.
int NucleotideToInt (char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default: return -1;
    }
}

/// Compute the prefix of size d of the given k-mer.
template <typename kmer_t>
kmer_t BitPrefix(kmer_t kMer, int k, int d) {
    return kMer >> ((k - d) << kmer_t(1));
}

/// Compute the suffix of size d of the given k-mer.
template <typename kmer_t>
kmer_t BitSuffix(kmer_t kMer, int d) {
    return kMer & ((kmer_t(1) << (d << kmer_t(1))) - kmer_t(1));
}

/// Checkered mask. cmask<uint16_t, 1> is every other bit on
/// (0x55). cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
/// Copyright: Jellyfish GPL-3.0
template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
struct cmask {
    static const U v =
            (cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
};
template<typename U, int len>
struct cmask<U, len, 0> {
    static const U v = 0;
};

/// Compute the reverse complement of a word.
/// Copyright: Jellyfish GPL-3.0
inline kmer64_t word_reverse_complement(kmer64_t w) {
    typedef kmer64_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ( w >> 32                   ) | ( w                    << 32);
    return ((U)-1) - w;
}
/// Compute the reverse complement of a word.
/// Copyright: Jellyfish GPL-3.0
inline kmer128_t word_reverse_complement(kmer128_t w) {
    typedef kmer128_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ((w >> 32) & cmask<U, 32>::v) | ((w & cmask<U, 32>::v) << 32);
    w = ( w >> 64                   ) | ( w                    << 64);
    return ((U)-1) - w;
}
/// Compute the reverse complement of a word.
/// Copyright: Jellyfish GPL-3.0
inline kmer256_t word_reverse_complement(kmer256_t w) {
    kmer128_t low = word_reverse_complement(w.lower());
    kmer128_t high = word_reverse_complement(w.upper());
    return kmer256_t(low, high);
}

constexpr int KMER_SIZE_64 = 64;
constexpr int KMER_SIZE_128 = 128;
constexpr int KMER_SIZE_256 = 256;
#define INIT_KMERS(type)                                                                                                              \
                                                                                                                                      \
/* Get the mask to mask k-mers. */                                                                                                    \
inline kmer##type##_t MaskForK##type(int k) {                                                                                         \
    kmer##type##_t mask = (((kmer##type##_t) 1) << ((k << 1) - 1));                                                                   \
    return mask | (mask - 1);                                                                                                         \
}                                                                                                                                     \
                                                                                                                                      \
/* Compute the reverse complement of the given k-mer. */                                                                              \
inline kmer##type##_t ReverseComplement(kmer##type##_t kMer, int k) {                                                                 \
    return (((kmer##type##_t)word_reverse_complement(kMer)) >> (KMER_SIZE_##type - (k << kmer##type##_t(1)))) & MaskForK##type(k);    \
}                                                                                                                                     \

INIT_KMERS(64)
INIT_KMERS(128)
INIT_KMERS(256)

/// Return the lexicographically smaller of the k-mer and its reverse complement.
template <typename kmer_t>
inline kmer_t CanonicalKMer(kmer_t kMer, int k) {
    kmer_t rev = ReverseComplement(kMer, k);
    return kMer < rev ? kMer : rev;
}

const char letters[4] {'A', 'C', 'G', 'T'};

/// Return the index-th nucleotide from the encoded k-mer.
template <typename kmer_t>
inline char NucleotideAtIndex(kmer_t encoded, int k, int index) {
    return letters[(uint32_t)(encoded >> ((k - index - kmer_t(1)) << kmer_t(1))) & kmer_t(3)];
}

/// Convert the encoded KMer representation to string.
template <typename kmer_t>
std::string NumberToKMer(kmer_t encoded, int length) {
    std::string ret(length, 'N');
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret[length - i -1] = letters[(uint32_t)(encoded & 3)];
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}
