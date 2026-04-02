
.PHONY: all clean

main.exe:
	mkdir -p build
	c++ src/main.cpp src/indexer.cpp src/reader.cpp -std=c++17 -O3 -march=native -I include -o build/main.exe

clean:
	rm -rf build
test_indexer: main.exe
	./build/main.exe indexer -R ../Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I ../Data --kmer-size 10 --ref-chunk-size 10000000

all: test_indexer