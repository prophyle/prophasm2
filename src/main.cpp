#include <iostream>
#include <string>

#include "unistd.h"
#include "version.h"
#include "prophasm.h"
#include "parser.h"
#include "kthread.h"


#ifdef LARGE_KMERS
    constexpr int MAX_K = 64;
    const std::string VARIANT = "(128bit k-mer variant)";
#else
    constexpr int MAX_K = 32;
#endif


int Help() {
    std::cerr <<
              "\n" <<
              "Program:  prophasm2 (a greedy assembler for k-mer set compression)\n" <<
              "Version:  " << VERSION << "\n" <<
              "Contact:  Ondrej Sladky <ondra.sladky@gmail.com>\n" <<
              "\n" <<
              "Usage:    prophasm2 [options]\n" <<
              "\n" <<
              "Examples: prophasm2 -k 15 -i f1.fa -i f2.fa -x fx.fa\n" <<
              "             - compute intersection of f1 and f2\n" <<
              "          prophasm2 -k 15 -i f1.fa -i f2.fa -x fx.fa -o g1.fa -o g2.fa\n" <<
              "             - compute intersection of f1 and f2, and subtract it from them\n" <<
              "          prophasm2 -k 15 -i f1.fa -o g1.fa\n" <<
              "             - re-assemble f1 to g1\n" <<
              "\n" <<
              "Command-line parameters:\n" <<
              " -k INT   K-mer size.\n" <<
              " -i FILE  Input FASTA file (can be used multiple times).\n" <<
              " -o FILE  Output FASTA file (if used, must be used as many times as -i).\n" <<
              " -x FILE  Compute intersection, subtract it, save it.\n" <<
              " -s FILE  Output file with k-mer statistics.\n" <<
              " -S       Silent mode.\n" <<
              " -u       Do not consider k-mer and its reverse complement as equivalent.\n" <<
              "\n" <<
              "Note that '-' can be used for standard input/output. \n" <<
              std::endl;
    return 1;
}

void Version() {
    std::cerr << VERSION << std::endl;
}

void TestFile(FILE *fo, std::string fn) {
    if (fo == nullptr) {
        std::cerr << "Error: file '" << fn << "' could not be open (error " << errno << ", " << strerror(errno) << ")." << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv) {
    int32_t k = -1;

    std::string intersectionPath;
    std::vector<std::string> inPaths;
    std::vector<std::string> outPaths;
    std::string statsPath;
    FILE *fstats = nullptr;

    bool computeIntersection = false;
    bool computeOutput = false;
    bool verbose = true;
    bool complements = true;
    int threads = 1;

    if (argc<2) {
        Help();
        return 1;
    }
    int c;
    while ((c = getopt(argc, (char *const *)argv, "hSi:o:x:s:k:uvt:")) >= 0) {
        switch (c) {
            case 'h': {
                return Help();
            }
            case 'i': {
                inPaths.push_back(std::string(optarg));
                break;
            }
            case 'o': {
                outPaths.push_back(std::string(optarg));
                computeOutput = true;
                break;
            }
            case 'x': {
                intersectionPath = std::string(optarg);
                computeIntersection = true;
                break;
            }
            case 's': {
                statsPath=std::string(optarg);
                if (statsPath == "-") {
                    fstats = stdin;
                } else {
                    fstats = fopen(statsPath.c_str(), "w+");
                    TestFile(fstats, statsPath);
                }
                break;
            }
            case 'S': {
                verbose = false;
                break;
            }
            case 'k': {
                k = atoi(optarg);
                break;
            }
            case 't': {
                threads = atoi(optarg);
                break;
            }
            case 'u': {
                complements = false;
                break;
            }
            case 'v': {
                Version();
                return 0;
            }
            case '?': {
                std::cerr << "Unknown error" << std::endl;
                return 1;
            }
        }
    }
    if (k == -1) {
        std::cerr << "K-mer size (-k) is required." << std::endl;
        return Help();
    } else if (k <= 0 || MAX_K < k) {
        std::cerr << "K-mer size must satisfy 1 <= k <= " << MAX_K << "." << std::endl;
        return Help();
    }
    size_t setCount = inPaths.size();
    if (computeOutput && (outPaths.size() != setCount)) {
        std::cerr << "If -o is used, it must be used as many times as -i (" << setCount << "!=" << outPaths.size() << ")." << std::endl;
        return Help();
    }
    if (computeIntersection && (setCount < 2)) {
        std::cerr << "If -x is used, at least two sets must be provided." << std::endl;
        return Help();
    }
    if (setCount < 1) {
        std::cerr << "At least one input set must be provided." << std::endl;
        return Help();
    }
    if (threads < 1) {
        std::cerr << "Number of threads must be at least 1." << std::endl;
        return Help();
    }

    if (fstats) {
        fprintf(fstats,"# cmd: %s",argv[0]);

        for (int32_t i = 1; i < argc; i++){
            fprintf(fstats," %s",argv[i]);
        }
        fprintf(fstats,"\n");
    }

    if (verbose) {
        std::cerr << "=====================" << std::endl;
        std::cerr << "1) Loading references" << std::endl;
        std::cerr << "=====================" << std::endl;
    }
    std::vector<kh_S64_t*> fullSets(setCount);
    std::vector<size_t> inSizes = std::vector<size_t>(setCount);
    std::vector<size_t> outSizes;




    for (size_t i = 0; i < setCount; i++) {
        fullSets[i] = kh_init_S64();
    }

    ReadKMersData data = {fullSets, inPaths, k, complements};
    kt_for(threads, ReadKMersThread, (void*)&data, setCount);

    for (size_t i = 0; i < setCount; i++) {
        if (verbose) {
            std::cerr << "Loaded " << inPaths[i] << std::endl;
        }
        inSizes[i] = kh_size(fullSets[i]);
        if (fstats != nullptr) {
            fprintf(fstats,"%s\t%lu\n", inPaths[i].c_str(), inSizes[i]);
        }
    }

    if (verbose) {
        std::cerr << "===============" << std::endl;
        std::cerr << "2) Intersecting" << std::endl;
        std::cerr << "===============" << std::endl;
    }
    kh_S64_t* intersection = nullptr;
    size_t intersectionSize = 0;
    if (computeIntersection) {
        if (verbose) {
            std::cerr << "2.1) Computing intersection" << std::endl;
        }
        intersection = getIntersection(fullSets, k, complements);
        intersectionSize  = kh_size(intersection);
        if (verbose) {
            std::cerr << "   intersection size: " <<  intersectionSize << std::endl;
        }
        if (computeOutput) {
            if (verbose) {
                std::cerr << "2.2) Removing this intersection from all k-mer sets" << std::endl;
            }
            DifferenceInPlaceData data = {fullSets, intersection, k, complements};
            kt_for(threads, DifferenceInPlaceThread, (void*)&data, setCount);
        }
    }
    if (computeOutput) {
        for (size_t i = 0; i < setCount; i++) {
            outSizes.push_back(kh_size(fullSets[i]));
            if (inSizes[i] != outSizes[i] + intersectionSize) {
                std::cerr << "Internal error: k-mer set sizes do not correspond " << inSizes[i] << " != " << outSizes[i] << " + " << intersectionSize << std::endl;
                return 1;
            }
            if (verbose){
                std::cerr << inSizes[i] << " " << outSizes[i] << " ...inter:" << intersectionSize << std::endl;
            }
        }
    }

    if (verbose){
        std::cerr << "=============" << std::endl;
        std::cerr << "3) Assembling" << std::endl;
        std::cerr << "=============" << std::endl;
    }
    if (computeOutput) {
        std::vector<std::ostream*> ofs (setCount);
        std::vector<std::ofstream> filestreams (setCount);

        for (size_t i = 0; i < setCount; i++) {
            if (outPaths[i] != "-") {
                filestreams[i] = std::ofstream(outPaths[i]);
                ofs[i] = &filestreams[i];
            } else {
                ofs[i] = &std::cout;
            }
            if (fstats) {
                fprintf(fstats,"%s\t%lu\n", outPaths[i].c_str(), outSizes[i]);
            }
        }
        ComputeSimplitigsData data = {fullSets, ofs, k, complements, std::vector<int>(setCount)};
        kt_for(threads, ComputeSimplitigsThread, (void*)&data, setCount);

        for (size_t i = 0; i < setCount; i++) {
            if (verbose) {
                std::cerr << "   assembly finished (" << data.simplitigsCounts[i] << " contigs)" << std::endl;
            }
            if (filestreams[i].is_open()) {
                filestreams[i].close();
            }
        }
    }
    if (computeIntersection) {
        std::ostream *of;
        std::ofstream filestream;
        if (intersectionPath != "-") {
            filestream = std::ofstream(intersectionPath);
            of = &filestream;
        }
        else {
            of = &std::cout;
        }
        if (fstats) {
            fprintf(fstats,"%s\t%lu\n", intersectionPath.c_str(), intersectionSize);
        }
        int simplitigCount = ComputeSimplitigs(intersection, *of, k, complements);
        if (verbose) {
            std::cerr << "   assembly finished (" << simplitigCount << " contigs)" << std::endl;
        }
        if (filestream.is_open()) {
            filestream.close();
        }
    }
    if (fstats){
        fclose(fstats);
    }
    return 0;
}
