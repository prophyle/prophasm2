#include <iostream>
#include <string>

#include "unistd.h"
#include "version.h"
#include "prophasm.h"
#include "parser.h"


#ifdef LARGE_KMERS
    constexpr int MAX_K = 63;
    const std::string VARIANT = "(128bit k-mer variant)";
#else
    constexpr int MAX_K = 31;
#endif


int Help() {
    std::cerr << "ProphAsm2 version " << VERSION << std::endl;
    std::cerr << "Here will be help" << std::endl;
    return 1;
}

void Version() {
    std::cerr << VERSION << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k = 0;
    std::ofstream output;
    std::ostream *of = &std::cout;
    bool complements = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "p:k:d:a:o:hcvm"))  != -1) {
            switch(opt) {
                case  'p':
                    path = optarg;
                    break;
                case 'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case  'c':
                    complements = true;
                    break;
                case 'v':
                    Version();
                    return 0;
                case 'h':
                default:
                    Help();
                    return 0;
            }
        }
    } catch (std::invalid_argument&) {
        return Help();
    }
    if (path.empty()) {
        std::cerr << "Required parameter p not set." << std::endl;
        return Help();
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return Help();
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return Help();
    } else if (k > MAX_K) {
        std::cerr << "k > " << MAX_K << " not supported." << std::endl;
        return Help();
    }

    kh_S64_t *kMers = kh_init_S64();
    ReadKMers(kMers, path, k, complements);
    if (!kh_size(kMers)) {
        std::cerr << "Path '" << path << "' contains no k-mers." << std::endl;
        return Help();
    }
    ComputeSimplitigs(kMers, *of, k, complements);
    return 0;
}
