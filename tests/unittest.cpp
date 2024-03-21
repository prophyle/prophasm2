#include "../src/kmers.h"
typedef kmer64_t  kmer_t;

#include "kmers_unittest.h"
#include "prophasm_unittest.h"
#include "khash_utils_unittest.h"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
