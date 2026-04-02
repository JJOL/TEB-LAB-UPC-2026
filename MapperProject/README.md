# MapperProject

### Objective:
Create a reads-sequence mapper that is able to assign optimal positions of a set of DNA reads onto a reference genome and with controllable parameters.

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


### Requirements
* Make (Tested with GNU Make 3.81)
* gcc/clang (Tested with Apple clang version 17.0.0 (clang-1700.6.4.2))