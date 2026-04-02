#!/bin/bash

# Chromosome in argument $1
# validate the input chromosome is of [1-22, X, Y]
if [[ ! $1 =~ ^([1-9]|1[0-9]|2[0-2]|X|Y)$ ]]; then
  echo "Invalid chromosome: $1. Please provide a valid chromosome number (1-22) or X/Y."
  exit 1
fi
CHROMOSOME_FILE_ZIPPED="chr$1.fna.gz"
CHROMOSOME_FILE="chr$1.fna"

gzip -d -c Data/$CHROMOSOME_FILE_ZIPPED > Data/$CHROMOSOME_FILE