
.PHONY: all clean

main.exe:
	mkdir -p build
	c++ src/main.cpp src/indexer.cpp src/reader.cpp src/mapper.cpp src/aligner.cpp src/writer.cpp src/utils.cpp -std=c++17 -O3 -march=native -I include -o build/main.exe

clean:
	rm -rf build
test_indexer: main.exe
	#./build/main.exe indexer -R Data/small.multiline.fna -I Data/smallidx --kmer-size 10 --ref-chunk-size 320
	./build/main.exe indexer -R Data/chr22.fna -I Data/chr22idx --kmer-size 10 --ref-chunk-size 10000000
test_mapper: main.exe
	#./build/main.exe mapper -R Data/chr22.fna -I Data/GRCh38_22_idx -i Data/reads_1.fastq -O Data/reads_1M.sam -k 10 --kmer-size 10 --ref-chunk-size 10000000
	#./build/main.exe mapper -R Data/small.multiline.fna -I Data/smallidx -i Data/small_reads.fastq -O Data/small.sam -S "NC_000003.12 Homo sapiens chromosome 3, GRCh38.p14 Primary Assembly" -k 10 --kmer-size 10 --ref-chunk-size 320
	#./build/main.exe mapper -R Data/small.multiline.fna -I Data/smallidx -i Data/light3900x101_reads.fastq -O Data/light3900x101.sam -k 50 --kmer-size 10 --ref-chunk-size 320
	#python pyMapper/analyze.py verify_test_sam -O Data/light3900x101.sam
	./build/main.exe mapper -R Data/chr22.fna -I Data/chr22idx -i Data/reads_10.fastq -O Data/reads_10.sam -k 4 --kmer-size 10 --ref-chunk-size 10000000

test_chromosomer: main.exe
	./build/main.exe chromosomer -R Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa

all: test_indexer