#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include "kmers.h"
#include "khash_utils.h"


/// Read encoded k-mers from the given fasta file.
/// Return unique k-mers in no particular order.
/// If complements is set to true, the result contains only one of the complementary k-mers - it is not guaranteed which one.
/// This runs in O(sequence length) expected time.
void ReadKMers(kh_S64_t *kMers, std::string &path, int k, bool complements) {
    std::ifstream fasta(path);
    if (fasta.is_open()) {
        char c;
        int beforeKMerEnd = k;
        kmer_t currentKMer = 0;
        // mask that works even for k=32.
        kmer_t mask = (((kmer_t) 1) <<  (2 * k - 1));
        mask |= mask - 1;
        bool readingHeader = false;
        while (fasta >> std::noskipws >> c) {
            if (c == '>') {
                readingHeader = true;
                currentKMer = 0;
                beforeKMerEnd = k;
            }
            else if (c == '\n') readingHeader = false;
            if (readingHeader) continue;
            auto data = NucleotideToInt(c);
            // Disregard white space.
            if (c == '\n' || c == '\r' || c == ' ') continue;
            if (data == -1) {
                currentKMer = 0;
                beforeKMerEnd = k;
                continue;
            }
            currentKMer <<= 2;
            currentKMer &= mask;
            currentKMer |= data;
            if(beforeKMerEnd > 0) --beforeKMerEnd;
            if (beforeKMerEnd == 0 && (!complements || kh_get_S64(kMers, ReverseComplement(currentKMer, k)) == kh_end(kMers))) {
                int ret;
                // If the k-mer was masked as present.
                kh_put_S64(kMers, currentKMer, &ret);
            }
        }
        fasta.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
}
