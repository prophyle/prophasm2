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
    KHASH_MAP_INIT_INT64(S64, char)
    KHASH_MAP_INIT_INT64(P64, size_t)
    KHASH_MAP_INIT_INT64(O64, size_t)
#endif


/// Determine whether the canonical k-mer is present.
bool containsKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    if (complements) kMer = CanonicalKMer(kMer, k);
    return kh_get_S64(kMers, kMer) != kh_end(kMers);
}

/// Remove the canonical k-mer from the set.
void eraseKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    if (complements) kMer = CanonicalKMer(kMer, k);
    auto key = kh_get_S64(kMers, kMer);
    if (key != kh_end(kMers)) {
        kh_del_S64(kMers, key);
    }
}

/// Insert the canonical k-mer into the set.
void insertKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    if (complements) kMer = CanonicalKMer(kMer, k);
    int ret;
    kh_put_S64(kMers, kMer, &ret);
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

/// Compute the intersection of several k-mer sets.
kh_S64_t *getIntersection(std::vector<kh_S64_t*> &kMerSets, int k, bool complements) {
    kh_S64_t* result = kh_init_S64();
    if (kMerSets.size() < 2) return result;
    kh_S64_t* smallestSet = kMerSets[0];
    for (size_t i = 1; i < kMerSets.size(); ++i) {
        if (kh_size(kMerSets[i]) < kh_size(smallestSet)) smallestSet = kMerSets[i];
    }
    for (auto i = kh_begin(smallestSet); i != kh_end(smallestSet); ++i) {
        if (!kh_exist(smallestSet, i)) continue;
        auto kMer = kh_key(smallestSet, i);
        bool everywhere = true;
        for (size_t i = 0; i < kMerSets.size(); ++i) if (kMerSets[i] != smallestSet) {
            if (!containsKMer(kMerSets[i], kMer, k, complements)) {
                everywhere = false;
                break;
            }
        }
        if (everywhere) {
            insertKMer(result, kMer, k, complements);
        }
    }
    return result;
}

/// Subtract the intersection from each k-mer set.
void differenceInPlace(kh_S64_t* kMerSet, kh_S64_t* intersection, int k, bool complements) {
    for (auto i = kh_begin(intersection); i != kh_end(intersection); ++i) {
        if (!kh_exist(intersection, i)) continue;
        auto kMer = kh_key(intersection, i);
        eraseKMer(kMerSet, kMer, k, complements);
    }
}

/// Data for parallel computation of set differences.
struct DifferenceInPlaceData {
    std::vector<kh_S64_t*> kMerSets;
    kh_S64_t* intersection;
    int k;
    bool complements;
};

/// Parallel wrapper for differenceInPlace.
void DifferenceInPlaceThread(void *arg, long i, int _) {
    auto data = (DifferenceInPlaceData*)arg;
    differenceInPlace(data->kMerSets[i], data->intersection, data->k, data->complements);
}