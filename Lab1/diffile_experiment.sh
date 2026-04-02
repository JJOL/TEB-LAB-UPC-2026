#!/bin/bash

# run 5 times sequentially ./t21 with different files (alternating chr8 and chrX) and output all into diffile.timestamp.txt

timestamp=$(date +%Y%m%d_%H%M%S)
output_file="diffile_${timestamp}.txt"
files=("../Data/chr8.fna" "../Data/chrX.fna")

echo "Running t21 5 times with alternating files (chr8, chrX), outputting to $output_file"
for i in {1..10}; do
    # Determine which file to use (0-indexed, so 0=chr8, 1=chrX)
    file_index=$(( (i - 1) % 2 ))
    input_file="${files[$file_index]}"
    
    echo "Run $i: $input_file" >> "$output_file"
    ./t21 "$input_file" >> "$output_file"
    echo -e "\n\n" >> "$output_file" # separate runs by blank lines
done
echo "Done. Output written to $output_file"
