#!/bin/bash

# List of Chromosomes Repo
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
# Chromosome in argument $1
# validate the input chromosome is of [1-22, X, Y]
if [[ ! $1 =~ ^([1-9]|1[0-9]|2[0-2]|X|Y)$ ]]; then
  echo "Invalid chromosome: $1. Please provide a valid chromosome number (1-22) or X/Y."
  exit 1
fi
CHROMOSOME_FILE="chr$1.fna.gz"


# Place file in Data/ directory

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/$CHROMOSOME_FILE \
    -O Data/$CHROMOSOME_FILE