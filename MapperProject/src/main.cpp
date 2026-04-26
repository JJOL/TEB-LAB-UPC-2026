#include <iostream>
#include <unordered_map>
#include "indexer.hpp"
#include "mapper.hpp"
#include "reader.hpp"
#define DEBUG

void parseArguments(int argc, char **argv, std::unordered_map<std::string, std::string> &argsMap);
void printUsage();
void printArguments(int argc, char **argv, const std::unordered_map<std::string, std::string> &argsMap);

int main(int argc, char **argv) {
    std::cout << "SeqMapperAligner v0.0.0" << std::endl;
    std::unordered_map<std::string, std::string> argsMap;
    parseArguments(argc, argv, argsMap);
#ifdef DEBUG
    printArguments(argc, argv, argsMap);
#endif

    // if help option is present, print usage and exit
    if (argsMap.find("help") != argsMap.end()) {
        printUsage();
        exit(0);
    }

    if (argsMap["subcommand"] == "indexer") {
        std::cout << "Indexing reference genome..." << std::endl;
        // Call indexing function here
        createIndex(argsMap);
    } else if (argsMap["subcommand"] == "mapper") {
        std::cout << "Aligning reads to reference genome..." << std::endl;
        mapReads(argsMap);
    } else if (argsMap["subcommand"] == "chromosomer") {
        std::cout << "Looking for chromosomes in reference genome file: " << argsMap["reference"] << std::endl;
        auto chromosomes = findChromosomesInFASTA(argsMap["reference"]);

        std::cout << "Found " << chromosomes.size() << " chromosomes in the reference genome." << std::endl;
        for (const auto &chromosome : chromosomes) {
            std::cout << "Chromosome: " << chromosome.name << ", Length: " << chromosome.length << ", Offset: " << chromosome.offset << std::endl;
        }
    } else {
        std::cerr << "Unknown subcommand: " << argsMap["subcommand"] << std::endl;
        printUsage();
        exit(1);
    }


    return 0;
}

void parseArguments(int argc, char **argv, std::unordered_map<std::string, std::string> &argsMap) {
    if (argc < 2) {
        printUsage();
        exit(0);
    }

    std::string subcommand = argv[1];
    argsMap["subcommand"] = subcommand;

    std::unordered_map<std::string, std::string> optionNameMap = {
        {"--help", "help"},
        {"-h", "help"},
        {"-R", "reference"},
        {"-I", "index_path"},
        {"-i", "reads"},
        {"-O", "sam_output_path"},
        {"-S", "sequences_whitelist"},
        {"--kmer-size", "kmer_size"},
        {"--ref-chunk-size", "chunk_size"},
        // Add more options as needed
    };

    for (int i = 2; i < argc-1; i += 2) {
        std::string key = argv[i];
        std::string value = argv[i + 1];
        if (optionNameMap.find(key) == optionNameMap.end()) {
            std::cerr << "Unknown option: " << key << std::endl;
            printUsage();
            exit(1);
        }
        argsMap[optionNameMap[key]] = value;
    }
}

void printUsage() {
    std::cout << "Usage: MapperProject [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --help, -h       Show this help message and exit" << std::endl;
    std::cout << "  -R               Reference genome file" << std::endl;
    std::cout << "  -I               Index file path" << std::endl;
    std::cout << "  -O               SAM output file path" << std::endl;
    std::cout << "  -i               Reads file (FASTQ)" << std::endl;
    std::cout << "  -S               Reference whitelist sequence names" << std::endl;
    std::cout << "  --kmer-size      K-mer size" << std::endl;
    std::cout << "  --ref-chunk-size Reference chunk size" << std::endl;
}

void printArguments(int argc, char **argv, const std::unordered_map<std::string, std::string> &argsMap) {
    std::cout << "Parsed Arguments:" << std::endl;
    for (const auto &[key, value] : argsMap) {
        std::cout << "  " << key << ": " << value << std::endl;
    }
}