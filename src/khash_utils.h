#pragma once

#include <vector>
#include <list>

#include "kmers.h"
#include "khash.h"

typedef unsigned char byte;

#define kh_int128_hash_func(key) kh_int64_hash_func((khint64_t)((key)>>65^(key)^(key)<<21))
#define kh_int128_hash_equal(a, b) ((a) == (b))
#define kh_int256_hash_func(key) kh_int128_hash_func((__uint128_t)((key)>>129^(key)^(key)<<35))
#define kh_int256_hash_equal(a, b) ((a) == (b))

#define KHASH_MAP_INIT_INT128(name, khval_t)								\
	KHASH_INIT(name, __uint128_t, khval_t, 1, kh_int128_hash_func, kh_int128_hash_equal)

#define KHASH_SET_INIT_INT128(name)										\
	KHASH_INIT(name, __uint128_t, byte, 0, kh_int128_hash_func, kh_int256_hash_equal)

#define KHASH_MAP_INIT_INT256(name, khval_t)								\
	KHASH_INIT(name, uint256_t, khval_t, 1, kh_int256_hash_func, kh_int128_hash_equal)

#define KHASH_SET_INIT_INT256(name)										\
	KHASH_INIT(name, uint256_t, byte, 0, kh_int256_hash_func, kh_int256_hash_equal)

KHASH_MAP_INIT_INT256(S256M, byte)
KHASH_MAP_INIT_INT128(S128M, byte)
KHASH_MAP_INIT_INT64(S64M, byte)
KHASH_SET_INIT_INT256(S256S)
KHASH_SET_INIT_INT128(S128S)
KHASH_SET_INIT_INT64(S64S)

byte MINIMUM_ABUNDANCE = 1;

// Forward definition for the macro.
template <typename KHT>
void differenceInPlace(KHT* kMerSet, KHT* intersection, int k, bool complements);

#define INIT_KHASH_UTILS(type, variant)                                                                                      \
                                                                                                                    \
/* Data for parallel computation of set differences. */                                                             \
struct DifferenceInPlaceData##variant {                                                                                \
    std::vector<kh_S##variant##_t*> kMerSets;                                                                          \
    kh_S##variant##_t* intersection;                                                                                   \
    int k;                                                                                                          \
    bool complements;                                                                                               \
};                                                                                                                  \
                                                                                                                    \
/* Determine whether the canonical k-mer is present.*/                                                              \
bool containsKMer(kh_S##variant##_t *kMers, kmer##type##_t kMer, int k, bool complements) {                                    \
    if (complements) kMer = CanonicalKMer(kMer, k);                                                                 \
    bool contains_key = kh_get_S##variant(kMers, kMer) != kh_end(kMers);                                               \
    if (MINIMUM_ABUNDANCE == 1) return contains_key;                                                                \
    if (!contains_key) return false;                                                                                \
    return kh_val(kMers, kh_get_S##variant(kMers, kMer)) >= MINIMUM_ABUNDANCE;                                         \
}                                                                                                                   \
                                                                                                                    \
/* Remove the canonical k-mer from the set.*/                                                                       \
void eraseKMer(kh_S##variant##_t *kMers, kmer##type##_t kMer, int k, bool complements) {                                       \
    if (complements) kMer = CanonicalKMer(kMer, k);                                                                 \
    auto key = kh_get_S##variant(kMers, kMer);                                                                         \
    if (key != kh_end(kMers)) {                                                                                     \
        kh_del_S##variant(kMers, key);                                                                                 \
    }                                                                                                               \
}                                                                                                                   \
                                                                                                                    \
/* Insert the canonical k-mer into the set. */                                                                      \
void insertKMer(kh_S##variant##_t *kMers, kmer##type##_t kMer, int k, bool complements, bool force = false) {                  \
    if (complements) kMer = CanonicalKMer(kMer, k);                                                                 \
    int ret;                                                                                                        \
    if (MINIMUM_ABUNDANCE == (byte)1) {                                                                             \
        kh_put_S##variant(kMers, kMer, &ret);                                                                          \
    } else {                                                                                                        \
        byte value = 0;                                                                                             \
        khint_t key = kh_get_S##variant(kMers, kMer);                                                                  \
        if (key != kh_end(kMers)) {                                                                                 \
            value = kh_val(kMers, key);                                                                             \
        } else {                                                                                                    \
            key = kh_put_S##variant(kMers, kMer, &ret);                                                                \
        }                                                                                                           \
        if (force || value == (byte)255) kh_value(kMers, key) = (byte)255;                                          \
        else kh_value(kMers, key) = value + 1;                                                                      \
    }                                                                                                               \
}                                                                                                                   \
                                                                                                                    \
/* Parallel wrapper for differenceInPlace. */                                                                       \
void DifferenceInPlaceThread##variant(void *arg, long i, int _) {                                                      \
    auto data = (DifferenceInPlaceData##variant*)arg;                                                                  \
    differenceInPlace(data->kMerSets[i], data->intersection, data->k, data->complements);                           \
}                                                                                                                   \

INIT_KHASH_UTILS(64, 64S)
INIT_KHASH_UTILS(64, 64M)
INIT_KHASH_UTILS(128, 128S)
INIT_KHASH_UTILS(128, 128M)
INIT_KHASH_UTILS(256, 256S)
INIT_KHASH_UTILS(256, 256M)

/// Return the next k-mer in the k-mer set and update the index.
template <typename KHT, typename kmer_t>
kmer_t nextKMer(KHT *kMers, size_t &lastIndex, kmer_t &kMer) {
    for (size_t i = kh_begin(kMers) + lastIndex; i != kh_end(kMers); ++i, ++lastIndex) {
        if (!kh_exist(kMers, i)) continue;
        kMer = kh_key(kMers, i);
        if (MINIMUM_ABUNDANCE > 1) if (!containsKMer(kMers, kMer, 0, false)) continue;
        return true;
    }
    // No more k-mers.
    lastIndex = -1;
    return false;
}

/// Construct a vector of the k-mer set in an arbitrary order. Only for testing.
std::vector<kmer64_t> kMersToVec(kh_S64M_t *kMers) {
    std::vector<kmer64_t> result(kh_size(kMers));
    size_t index = 0;
    for (auto i = kh_begin(kMers); i != kh_end(kMers); ++i) {
        if (!kh_exist(kMers, i)) continue;
        if (MINIMUM_ABUNDANCE > 1) if (!containsKMer(kMers, kh_key(kMers, i), 0, false)) continue;
        result[index++] = kh_key(kMers, i);
    }
    result.resize(index);
    return result;
}

/// Compute the intersection of several k-mer sets.
template <typename KHT>
KHT *getIntersection(KHT* result, std::vector<KHT*> &kMerSets, int k, bool complements) {
    if (kMerSets.size() < 2) return result;
    KHT* smallestSet = kMerSets[0];
    for (size_t i = 1; i < kMerSets.size(); ++i) {
        if (kh_size(kMerSets[i]) < kh_size(smallestSet)) smallestSet = kMerSets[i];
    }
    for (auto i = kh_begin(smallestSet); i != kh_end(smallestSet); ++i) {
        if (!kh_exist(smallestSet, i)) continue;
        auto kMer = kh_key(smallestSet, i);
        if (!containsKMer(smallestSet, kMer, k, complements)) continue;
        bool everywhere = true;
        for (size_t i = 0; i < kMerSets.size(); ++i) if (kMerSets[i] != smallestSet) {
            if (!containsKMer(kMerSets[i], kMer, k, complements)) {
                everywhere = false;
                break;
            }
        }
        if (everywhere) {
            insertKMer(result, kMer, k, complements, true);
        }
    }
    return result;
}

/// Subtract the intersection from each k-mer set.
template <typename KHT>
void differenceInPlace(KHT* kMerSet, KHT* intersection, int k, bool complements) {
    for (auto i = kh_begin(intersection); i != kh_end(intersection); ++i) {
        if (!kh_exist(intersection, i)) continue;
        auto kMer = kh_key(intersection, i);
        eraseKMer(kMerSet, kMer, k, complements);
    }
}


