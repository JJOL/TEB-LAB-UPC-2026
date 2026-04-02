#!/bin/bash

# run 5 times sequentially ./t21 with chr3.fna and output all into samefile.timestamp.txt

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.fna>"
    exit 1
fi

input_file="$1"
timestamp=$(date +%Y%m%d_%H%M%S)
output_file="samefile_${timestamp}.txt"
echo "Running t21 on $input_file 5 times, outputting to $output_file"
for i in {1..5}; do
    echo "Run $i:" >> "$output_file"
    ./t21 "$input_file" >> "$output_file"
    echo -e "\n\n" >> "$output_file" # separate runs by blank lines
done
echo "Done. Output written to $output_file"X