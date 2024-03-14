#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <vector>

#include "kmers.h"
#include "khash_utils.h"


/// Read encoded k-mers from the given fasta file.
/// Return unique k-mers in no particular order.
/// If complements is set to true, the result contains only one of the complementary k-mers - it is not guaranteed which one.
/// This runs in O(sequence length) expected time.
void ReadKMers(kh_S64_t *kMers, std::string &path, int k, bool complements) {
    std::ifstream filestream;
    std::istream *fasta;
    if (path == "-") {
        fasta = &std::cin;
    } else {
        filestream.open(path);
        fasta = &filestream;
        if (!filestream.is_open()) {
            std::cerr << "Error: file '" << path << "' could not be open." << std::endl;
            exit(1);
        }
    }
    char c;
    int beforeKMerEnd = k;
    kmer_t currentKMer = 0;
    // mask that works even for k=32.
    kmer_t mask = (((kmer_t) 1) <<  (2 * k - 1));
    mask |= mask - 1;
    bool readingHeader = false;
    while ((*fasta) >> std::noskipws >> c) {
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
        if (beforeKMerEnd == 0) {
            insertKMer(kMers, currentKMer, k, complements);
        }
    }
    if (filestream.is_open()) filestream.close();
}

/// Data for parallel reading of k-mers.
struct ReadKMersData {
    std::vector<kh_S64_t*> kMers;
    std::vector<std::string> paths;
    int k;
    bool complements;
};

/// Parallel wrapper for ReadKMers.
void ReadKMersThread(void *arg, long i, int _) {
    auto *data = (ReadKMersData *) arg;
    ReadKMers(data->kMers[i], data->paths[i], data->k, data->complements);
}
