# MapperProject

### Objective:
Create a reads-sequence mapper that is able to assign optimal positions of a set of DNA reads onto a reference genome and with controllable parameters.

### Usage

Simple building and running with small test files.
```bash
cd MapperProject
make test_indexer
make test_mapper
```

Building
```bash
cd MapperProject
make main.exe
```

Running Indexer
```bash
./build/main.exe indexer -R ../Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I ../Data/GRCh38Idx --kmer-size 10 --ref-chunk-size 10000000
```

Running Mapper
```bash
./build/main.exe mapper -R ../Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I ../Data/GRCh38Idx -i ../Data/reads_1M.fastq -O ../Data/reads_1M.sam --kmer-size 10 --ref-chunk-size 10000000
```


#### Additional Parameters
* -S: For both indexer and mapper, limits to a subset of sequences in the reference genome. It should be a string with 1 or more sequence names separated by semicolon (;). Example: "NC_000008.11 Primary;NC_000003.12 Primary".
* --kmer-size: The size of K-mers used to build the indices.
* --ref-chunksize: The division used to split the reference chromosomes and improve peak memory usage.
* chromsomer: A subcommand for printing all the chromosomes found in a -R reference FASTA file, their length (in lines) and file line offset.

### Requirements
* Make (Tested with GNU Make 3.81)
* gcc/clang (Tested with Apple clang version 17.0.0 (clang-1700.6.4.2))
* python3.12 for some helper scripts

### Design Considerations
Simplicity assumptions:
* Input: A FASTA genome reference file and a FASTQ reads file.
* 1M reads in a FASTQ format of approximate 101 bases.
* Reads Qualities are ignored.
* 24 chromosomes of full human genome Chr38.
* If we divide chromsomes in chunks, reads are very small and usually do not fall between chunks. The ones that do will have slight higer penalty and slight offset to true position.
* All reads have alphabet {A,C,G,T}. If a read has 'N', it is not matched.
* All reference genome sequences have alphabet {A,C,G,T}. If another character is present it is treated as 'N'.
* Reference genome regions with 'N' are not matched.
* A simple DP edit distance scoring is used with matches costing 0, mismatches 1 and deletion/insertion 1.
* Final Output: A SAM format file.
* CIGAR operations string are in direction Reference -> Read.
* We report the TOP 2 best alignments for each read across the whole FASTA genome (if they exist).

Additional assumptions:
* Oriented for comodity hardware: personal or gaming laptos or edge-hardware acceleration like Nvidia Jetsons.
* Hopeful: Be able to map Whole Human Genome (3 giga bases) accross 23 chromosomes. Could be downscaled to simple shorter genomes like from bacteria for specific applications if 3 giga bases is too big for commodity hardware.
* Commodity hardware:
    * Gaming Laptop:
        - RAM: 16GB DDR4
        - Intel Core i7 (6 perf + 4 eff cores, 24MB-L3 Caches, 256AVX, 3GHz)
        - GPU: Nvidia RTX 3070 Ampere, 5129 CUDA cores, FLOPS, 8GB VRAM, 1.5GHz, GGDR6
    * Edge-Device: (NVIDIA Jetson Orin Nano Super Developer Kit)
        - RAM: 8 GB DDR4
        - ARM Cortex-A78AE v8.2 (6 cores, 4MB L3 Cache, 2GHz)
        - GPU: NVIDIA Ampere 1024 CUDA cores, 1.3GHz
* Primarly for long-reads technologies like Nanopore, modular for Illumina technology
* Reads Characteristics:
    - 50K - 500K bases
    - Error Rate: 5-15%
    - Delins: 60-70% of errors
    - \# reads: ~1M - 10M reads

### Additional Material
