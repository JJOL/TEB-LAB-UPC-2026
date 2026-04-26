"""
analyze.py
Takes a reference file and command name and performs analysis on the reference file.
Commmands: get_alphabet
Usage: python analyze.py <command> -R <reference_file>
"""
import sys
import argparse

def print_usage():
    print("Usage: python analyze.py <command> -R <reference_file>")
    print("Commands:")
    print("  get_alphabet: Get the alphabet used in the reference file.")
    print("  verify_test_sam: Verify a SAM file against the reference.")
    print()

def main():
    parser = argparse.ArgumentParser(description="Analyze a reference file.")
    parser.add_argument("command", help="Command to execute")
    parser.add_argument("-R", "--reference", help="Reference FASTA file", required=False)
    parser.add_argument("-i", "--reads", help="Reads FASTQ file", required=False)
    parser.add_argument("-O", "--sam_output", help="SAM file to verify", required=False)
    args = parser.parse_args()

    if args.command == "get_alphabet":
        # if refrence file, use it as fasta with 1 read. If reads file, is a fastq with arbitrary number of reads.
        if not args.reference and not args.reads:
            print("Error: You must provide either a reference file or a reads file.")
            print_usage()
            sys.exit(1)
        elif args.reference and args.reads:
            print("Error: You cannot provide both a reference file and a reads file.")
            print_usage()
            sys.exit(1)
        if args.reference:
            print("Getting alphabet for reference file:", args.reference)
            alphabet = get_alphabet_fasta(args.reference)
            print("Alphabet:", alphabet)
        else:
            print("Getting alphabet for reads file:", args.reads)
            alphabet = get_alphabet_fastq(args.reads)
            # some reads have N. Count them
            total_reads, n_reads = count_n_reads(args.reads)
            print(f"Total reads: {total_reads}, Reads with N: {n_reads}, Percentage with N: {n_reads/total_reads*100:.2f}%")
            print("Alphabet:", alphabet)
    elif args.command == "verify_test_sam":
        # -O SAM file to verify. Read each sequence read with format:
        # @read0000011|NC_000008.11:1123	NC_000008.11 Homo sapiens chromosome 8, GRCh38.p14 Primary Assembly	1123	61M1X39M	GGTTCTGTGAGACTGGTAGAAAGCACAGACCCCTTAGACTTCTCCCCAAGGAGAATACGTGGGACTAGTGGAGGAAAGGAGAGTAATGAAATATGCATTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	ALT:NC_000003.12 Homo sapiens chromosome 3, GRCh38.p14 Primary Assembly,149,1M1X4M1X2M3X2M1I2M1I1M1I2X1M1X2M1X1M4X1M1I1M1I1M1I5M1X2M1X1M3X3M1X1M1D1M1X4M1D3X3M1X1M1X2M2X2M4X1M1X2M1X1M1X3M1X1M1X1M1I1X1M2X1M
        # the read name has first read index (11), the reference sequence (8) and offset (1123). Then we have the matched ref seq name and then the matched offset.
        # count the mount of reads that match for the first option. The amount of reads that do not match for the first option, but do in the second.
        # finally, the number of reads that do not match nor the first neither second best matches.
        if not args.sam_output:
            print("Error: You must provide a SAM file to verify.")
            print_usage()
            sys.exit(1)
        else:
            print("Verifying SAM file:", args.sam_output)
            matches_1st, matches_2nd, no_matches = verify_sam(args.sam_output)
            print(f"Matches for first option: {matches_1st}")
            print(f"Matches for second option: {matches_2nd}")
            print(f"No matches: {no_matches}")
    else:
        print_usage()

def get_alphabet_fasta(reference_file):
    alphabet = set()
    # first line is the header, so skip it
    with open(reference_file, 'r') as f:
        next(f)  # Skip the header line
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                alphabet.update(line)
    return sorted(alphabet)

def get_alphabet_fastq(reads_file):
    alphabet = set()
    with open(reads_file, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line in FASTQ format
                line = line.strip()
                if line and not line.startswith("#"):
                    alphabet.update(line)
    return sorted(alphabet)

def count_n_reads(reads_file):
    total_reads = 0
    n_reads = 0
    with open(reads_file, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line in FASTQ format
                total_reads += 1
                if 'N' in line:
                    n_reads += 1
    return total_reads, n_reads

def verify_sam(sam_file):
    matches_1st = 0
    matches_2nd = 0
    no_matches = 0
    with open(sam_file, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            read_name = fields[0]
            ref_name = fields[1]
            ref_offset = int(fields[2])
            # get the reference name from the read name
            read_ref_name = read_name.split("|")[1].split(":")[0]
            read_ref_offset = int(read_name.split("|")[1].split(":")[1])
            ref_seq_name = ref_name.split(" ")[0]
            if read_ref_name == ref_seq_name and read_ref_offset == ref_offset:
                matches_1st += 1
            elif read_ref_name != ref_seq_name and len(fields) > 6 and "ALT:" in fields[6]:
                alt_ref_name = fields[6].split("ALT:")[1].split(" ")[0]
                alt_ref_offset = int(fields[6].split(",")[-2])
                if alt_ref_name == read_ref_name and alt_ref_offset == read_ref_offset:
                    matches_2nd += 1
                else:
                    no_matches += 1
            else:
                no_matches += 1

    return matches_1st, matches_2nd, no_matches

if __name__ == "__main__":
    main()