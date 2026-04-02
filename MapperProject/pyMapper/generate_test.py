"""Generate synthetic FASTA references and FASTQ reads for mapper testing.

Examples:
	python generate_test.py reference -o Data/test_ref.fna -n 3 -l 5000
	python generate_test.py reads -r Data/test_ref.fna -o Data/test_reads.fastq -n 1000 -l 120 -e 0.01
"""

import argparse
import gzip
import random
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

DNA_ALPHABET = "ACGT"
QUALITY_CHAR = "I"  # Phred+33 ~= 40, good synthetic default


def build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description="Synthetic data generator for ProjectMapper (FASTA + FASTQ)."
	)
	parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")

	subparsers = parser.add_subparsers(dest="command", required=True)

	reference = subparsers.add_parser("reference", help="Generate a reference FASTA file")
	reference.add_argument("-o", "--output", required=True, type=Path, help="Output FASTA path")
	reference.add_argument(
		"-n", "--num-sequences", type=int, default=1, help="Number of FASTA entries"
	)
	reference.add_argument(
		"-l", "--length", type=int, default=10000, help="Length of each sequence"
	)
	reference.add_argument(
		"--line-width",
		type=int,
		default=80,
		help="FASTA sequence line width (default: 80)",
	)

	reads = subparsers.add_parser("reads", help="Generate FASTQ reads")
	reads.add_argument(
		"-r",
		"--reference",
		required=True,
		type=Path,
		help="Reference FASTA(.gz) used as source",
	)
	reads.add_argument("-o", "--output", required=True, type=Path, help="Output FASTQ path")
	reads.add_argument("-n", "--num-reads", type=int, default=1000, help="Number of reads")
	reads.add_argument("-l", "--read-length", type=int, default=150, help="Nominal read length")
	reads.add_argument(
		"-e",
		"--error-rate",
		type=float,
		default=0.01,
		help="Per-base error probability (substitutions/insertions/deletions)",
	)
	reads.add_argument(
		"--indel-fraction",
		type=float,
		default=0.2,
		help=(
			"Fraction of errors that become indels (rest are substitutions). "
			"Default: 0.2"
		),
	)

	return parser


def parse_args() -> argparse.Namespace:
	parser = build_parser()
	args = parser.parse_args()

	if args.command == "reference":
		if args.num_sequences <= 0:
			parser.error("--num-sequences must be > 0")
		if args.length <= 0:
			parser.error("--length must be > 0")
		if args.line_width <= 0:
			parser.error("--line-width must be > 0")

	if args.command == "reads":
		if args.num_reads <= 0:
			parser.error("--num-reads must be > 0")
		if args.read_length <= 0:
			parser.error("--read-length must be > 0")
		if not (0.0 <= args.error_rate <= 1.0):
			parser.error("--error-rate must be in [0, 1]")
		if not (0.0 <= args.indel_fraction <= 1.0):
			parser.error("--indel-fraction must be in [0, 1]")

	return args


def _open_text_auto(path: Path, mode: str):
	if path.suffix == ".gz":
		return gzip.open(path, mode, encoding="utf-8")
	return path.open(mode, encoding="utf-8")


def _random_dna(length: int) -> str:
	return "".join(random.choices(DNA_ALPHABET, k=length))


def write_fasta(records: Iterable[Tuple[str, str]], output_path: Path, line_width: int) -> None:
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with _open_text_auto(output_path, "wt") as handle:
		for name, seq in records:
			handle.write(f">{name}\n")
			for i in range(0, len(seq), line_width):
				handle.write(seq[i : i + line_width] + "\n")


def read_fasta(path: Path) -> List[Tuple[str, str]]:
	sequences: List[Tuple[str, str]] = []
	current_name = None
	current_chunks: List[str] = []

	with _open_text_auto(path, "rt") as handle:
		for raw_line in handle:
			line = raw_line.strip()
			if not line:
				continue
			if line.startswith(">"):
				if current_name is not None:
					sequences.append((current_name, "".join(current_chunks).upper()))
				current_name = line[1:].split()[0]
				current_chunks = []
			else:
				current_chunks.append(line)

	if current_name is not None:
		sequences.append((current_name, "".join(current_chunks).upper()))

	if not sequences:
		raise ValueError(f"No FASTA sequences found in {path}")
	return sequences


def _mutate_read(seq: str, error_rate: float, indel_fraction: float) -> str:
	out: List[str] = []
	i = 0
	while i < len(seq):
		base = seq[i]
		if random.random() >= error_rate:
			out.append(base)
			i += 1
			continue

		if random.random() < indel_fraction:
			if random.random() < 0.5:
				# Deletion: skip reference base.
				i += 1
			else:
				# Insertion: add random base before current one.
				out.append(random.choice(DNA_ALPHABET))
		else:
			# Substitution: replace with a different base.
			choices = [b for b in DNA_ALPHABET if b != base]
			out.append(random.choice(choices))
			i += 1

	# Keep reads close to nominal length for easier downstream tests.
	if len(out) == 0:
		return random.choice(DNA_ALPHABET)
	return "".join(out)


def write_fastq(records: Iterable[Tuple[str, str]], output_path: Path) -> None:
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with _open_text_auto(output_path, "wt") as handle:
		for name, seq in records:
			qual = QUALITY_CHAR * len(seq)
			handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


def generate_reference(args: argparse.Namespace) -> None:
	records = [
		(f"chr{i + 1}", _random_dna(args.length))
		for i in range(args.num_sequences)
	]
	write_fasta(records, args.output, line_width=args.line_width)
	print(
		f"Reference generated: {args.output} "
		f"(num_sequences={args.num_sequences}, length={args.length})"
	)


def generate_reads(args: argparse.Namespace) -> None:
	references = read_fasta(args.reference)
	eligible = [(name, seq) for name, seq in references if len(seq) >= args.read_length]
	if not eligible:
		raise ValueError(
			f"No reference sequence is long enough for read_length={args.read_length}"
		)

	records: List[Tuple[str, str]] = []
	for i in range(args.num_reads):
		ref_name, ref_seq = random.choice(eligible)
		start = random.randint(0, len(ref_seq) - args.read_length)
		clean = ref_seq[start : start + args.read_length]
		noisy = _mutate_read(clean, args.error_rate, args.indel_fraction)
		read_name = f"read{i + 1:07d}|{ref_name}:{start}"
		records.append((read_name, noisy))

	write_fastq(records, args.output)
	print(
		f"Reads generated: {args.output} "
		f"(num_reads={args.num_reads}, read_length={args.read_length}, error_rate={args.error_rate})"
	)


def main() -> None:
	args = parse_args()
	random.seed(args.seed)

	if args.command == "reference":
		generate_reference(args)
		return

	if args.command == "reads":
		generate_reads(args)
		return


if __name__ == "__main__":
	main()