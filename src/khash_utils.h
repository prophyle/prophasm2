#pragma once

#include <vector>
#include <list>

#include "kmers.h"
#include "khash.h"


#define kh_int128_hash_func(key) kh_int64_hash_func((khint64_t)((key)>>65^(key)^(key)<<21))
#define kh_int128_hash_equal(a, b) ((a) == (b))

#define KHASH_MAP_INIT_INT128(name, khval_t)								\
	KHASH_INIT(name, __uint128_t, khval_t, 1, kh_int128_hash_func, kh_int128_hash_equal)

#define KHASH_SET_INIT_INT128(name)										\
	KHASH_INIT(name, __uint128_t, char, 0, kh_int128_hash_func, kh_int128_hash_equal)

#ifdef LARGE_KMERS
    // Use 128-bit integers for k-mers to allow for larger k.
    KHASH_SET_INIT_INT128(S64)
    KHASH_MAP_INIT_INT128(P64, size_t)
    KHASH_MAP_INIT_INT128(O64, size_t)
#else
    // Use 64-bits integers for k-mers for faster operations and less memory usage.
    KHASH_SET_INIT_INT64(S64)
    KHASH_MAP_INIT_INT64(P64, size_t)
    KHASH_MAP_INIT_INT64(O64, size_t)
#endif


/// Determine whether the k-mer or its reverse complement is present.
bool containsKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    bool ret = kh_get_S64(kMers, kMer) != kh_end(kMers);
    if (complements) ret |= kh_get_S64(kMers, ReverseComplement(kMer, k )) != kh_end(kMers);
    return ret;
}

/// Remove the k-mer and its reverse complement.
void eraseKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    auto key = kh_get_S64(kMers, kMer);
    if (key != kh_end(kMers)) {
        kh_del_S64(kMers, key);
    }
    if (complements) {
        kmer_t reverseComplement = ReverseComplement(kMer, k);
        key = kh_get_S64(kMers, reverseComplement);
        if (key != kh_end(kMers)) kh_del_S64(kMers, key);
    }
}

/// Return the next k-mer in the k-mer set and update the index.
kmer_t nextKMer(kh_S64_t *kMers, size_t &lastIndex, kmer_t &kMer) {
    for (size_t i = kh_begin(kMers) + lastIndex; i != kh_end(kMers); ++i, ++lastIndex) {
        if (!kh_exist(kMers, i)) continue;
        kMer = kh_key(kMers, i);
        return true;
    }
    // No more k-mers.
    lastIndex = -1;
    return false;
}

/// Construct a vector of the k-mer set in an arbitrary order.
std::vector<kmer_t> kMersToVec(kh_S64_t *kMers) {
    std::vector<kmer_t> res(kh_size(kMers));
    size_t index = 0;
    for (auto i = kh_begin(kMers); i != kh_end(kMers); ++i) {
        if (!kh_exist(kMers, i)) continue;
        res[index++] = kh_key(kMers, i);
    }
    return res;
}

