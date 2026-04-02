#include <iostream>
#include <string>
#include <vector>
// hash map container from string to vector of integers. It doesnt need to be ordered, so we can use unordered map
#include <unordered_map>
#include <algorithm>
#include <limits>

std::unordered_map<std::string, std::vector<int>> getKmersIndex(std::string text, int k = 3) {
    std::unordered_map<std::string, std::vector<int>> hash_map;
    for (size_t i = 0; i <= text.size() - k; i++) {
        std::string kmer = text.substr(i, k);
        hash_map[kmer].push_back(i);
    }

    return hash_map;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    // 100 characters random
    std::string T = "TACGATCGCGACGAGCTAGGCTAGCTGACTGCATGCATGCATGCTAGCTCGATGCTAGCTAGCTAGCTAGCTCTACTAGCTAGCTAGCTGCATGCATGCGCGCGGACGAGCTACGATCGTAGCTAGGTAGCATCGATGATTAGATCGTAGATGATGACTAGCTCGATCGCGCTACGATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGGCGCGATGCATAGCTAGCTCGTACGATCGATCGATTATATATTATATCGAGCGCGACGGCAGTCAGCATC";
    std::string P = "GCTAGCTCGATGCTAGCTAGCTAGCTAGCTCTACTAG";

    /* Task 2.1: K-mers map */
    std::cout << "==========================================================" << std::endl;
    std::cout << " Task 2.1: K-mers map" << std::endl;
    std::cout << "==========================================================" << std::endl;
    int k = 3;
    
    std::unordered_map<std::string, std::vector<int>> hash_map = getKmersIndex(T, k);
    int totalKmersInstances = T.size() - k + 1;
    std::cout << "Extracted " << totalKmersInstances << " kmers in total!" << std::endl;

    int size = hash_map.size();
    std::cout << "We found (" << size << "/" << (4*4*4) << ") distinct kmers!" << std::endl;
    std::string topKmer;
    int topKmerCount = -1;
    int sum = 0;
    for (auto kv : hash_map) {
        int count = kv.second.size();
        std::cout << " - " << kv.first << " has " << count << " ocurrences!" << std::endl;
        if (count > topKmerCount) {
            topKmerCount = count;
            topKmer = kv.first;
        }
        sum += count;
    }

    std::cout << "Most popular kmer: " << topKmer << " with " << topKmerCount << " occurrences." << std::endl;
    std::cout << "Sum of K-mer frequencies: " << sum << std::endl;


    /* Task 2.1: K-mers map */
    std::cout << "==========================================================" << std::endl;
    std::cout << " Task 2.2: K-mers map" << std::endl;
    std::cout << "==========================================================" << std::endl;
    auto patternKmers = getKmersIndex(P, k);
    int leastFreqKmerCount = INT32_MAX;
    std::string leastFreqKmer;
    for (auto kv : patternKmers) {
        int count = kv.second.size();
        if (count < leastFreqKmerCount) {
            leastFreqKmerCount = count;
            leastFreqKmer = kv.first;
        }
    }
    

    return 0;
}