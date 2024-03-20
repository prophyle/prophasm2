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

typedef unsigned char byte;

KHASH_MAP_INIT_INT128(S128, byte)
KHASH_MAP_INIT_INT64(S64, byte)

byte MINIMUM_ABUNDANCE = 1;

// Forward definition for the macro.
template <typename KHT>
void differenceInPlace(KHT* kMerSet, KHT* intersection, int k, bool complements);

#define INIT_KHASH_UTILS(type)                                                                                      \
                                                                                                                    \
/* Data for parallel computation of set differences. */                                                             \
struct DifferenceInPlaceData##type {                                                                                \
    std::vector<kh_S##type##_t*> kMerSets;                                                                          \
    kh_S##type##_t* intersection;                                                                                   \
    int k;                                                                                                          \
    bool complements;                                                                                               \
};                                                                                                                  \
                                                                                                                    \
/* Determine whether the canonical k-mer is present.*/                                                              \
bool containsKMer(kh_S##type##_t *kMers, kmer_t kMer, int k, bool complements) {                                    \
    if (complements) kMer = CanonicalKMer(kMer, k);                                                                 \
    bool contains_key = kh_get_S##type(kMers, kMer) != kh_end(kMers);                                               \
    if (MINIMUM_ABUNDANCE == 1) return contains_key;                                                                \
    if (!contains_key) return false;                                                                                \
    return kh_val(kMers, kh_get_S##type(kMers, kMer)) >= MINIMUM_ABUNDANCE;                                         \
}                                                                                                                   \
                                                                                                                    \
/* Remove the canonical k-mer from the set.*/                                                                       \
void eraseKMer(kh_S##type##_t *kMers, kmer_t kMer, int k, bool complements) {                                       \
    if (complements) kMer = CanonicalKMer(kMer, k);                                                                 \
    auto key = kh_get_S##type(kMers, kMer);                                                                         \
    if (key != kh_end(kMers)) {                                                                                     \
        kh_del_S##type(kMers, key);                                                                                 \
    }                                                                                                               \
}                                                                                                                   \
                                                                                                                    \
/* Insert the canonical k-mer into the set. */                                                                      \
void insertKMer(kh_S##type##_t *kMers, kmer_t kMer, int k, bool complements, bool force = false) {                  \
    if (complements) kMer = CanonicalKMer(kMer, k);                                                                 \
    int ret;                                                                                                        \
    if (MINIMUM_ABUNDANCE == (byte)1) {                                                                             \
        kh_put_S##type(kMers, kMer, &ret);                                                                          \
    } else {                                                                                                        \
        byte value = 0;                                                                                             \
        khint_t key = kh_get_S##type(kMers, kMer);                                                                  \
        if (key != kh_end(kMers)) {                                                                                 \
            value = kh_val(kMers, key);                                                                             \
        } else {                                                                                                    \
            key = kh_put_S##type(kMers, kMer, &ret);                                                                \
        }                                                                                                           \
        if (force || value == (byte)255) kh_value(kMers, key) = (byte)255;                                          \
        else kh_value(kMers, key) = value + 1;                                                                      \
    }                                                                                                               \
}                                                                                                                   \
                                                                                                                    \
/* Parallel wrapper for differenceInPlace. */                                                                       \
void DifferenceInPlaceThread##type(void *arg, long i, int _) {                                                      \
    auto data = (DifferenceInPlaceData##type*)arg;                                                                  \
    differenceInPlace(data->kMerSets[i], data->intersection, data->k, data->complements);                           \
}                                                                                                                   \

INIT_KHASH_UTILS(64)
INIT_KHASH_UTILS(128)

/// Return the next k-mer in the k-mer set and update the index.
template <typename KHT>
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

/// Construct a vector of the k-mer set in an arbitrary order.
template <typename KHT>
std::vector<kmer_t> kMersToVec(KHT *kMers) {
    std::vector<kmer_t> res(kh_size(kMers));
    size_t index = 0;
    for (auto i = kh_begin(kMers); i != kh_end(kMers); ++i) {
        if (!kh_exist(kMers, i)) continue;
        if (MINIMUM_ABUNDANCE > 1) if (!containsKMer(kMers, kh_key(kMers, i), 0, false)) continue;
        res[index++] = kh_key(kMers, i);
    }
    res.resize(index);
    return res;
}

/// Compute the intersection of several k-mer sets.
template <typename KHT>
KHT *getIntersection(KHT* result, std::vector<KHT*> &kMerSets, int k, bool complements) {
    if (kMerSets.size() < 2) return result;
    kh_S64_t* smallestSet = kMerSets[0];
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


