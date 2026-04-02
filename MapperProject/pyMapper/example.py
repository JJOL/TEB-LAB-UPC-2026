import argparse
from array import array
from bisect import bisect_right
from collections import Counter, deque
import gzip
import math
import pickle
from pathlib import Path
from typing import Deque, Dict, Iterable, List, Optional, Sequence, Tuple


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mapper",
        description=(
            "PoC CLI for the TEB mapper project. "
            "Supports mapper and indexer interfaces from the assignment."
        ),
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    mapper = subparsers.add_parser("mapper", help="Map FASTQ reads to a reference genome")
    mapper.add_argument(
        "-R",
        "--reference",
        type=Path,
        help="Reference genome in FASTA format (e.g. chr16.fna)",
    )
    mapper.add_argument(
        "-I",
        "--index",
        type=Path,
        help="Prebuilt index file for the reference genome",
    )
    mapper.add_argument(
        "-i",
        "--input-reads",
        required=True,
        type=Path,
        help="Input reads in FASTQ(.gz) format",
    )
    mapper.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output SAM file path",
    )
    mapper.add_argument(
        "-k",
        "--max-errors",
        required=True,
        type=int,
        help="Maximum allowed edit distance errors (typically 0..3)",
    )
    mapper.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="Optional cap of reads to process (0 means all reads)",
    )
    mapper.add_argument(
        "--max-candidates",
        type=int,
        default=5000,
        help="Cap candidate positions per read after seeding (default: 5000)",
    )

    indexer = subparsers.add_parser("indexer", help="Build index from a reference genome")
    indexer.add_argument(
        "-R",
        "--reference",
        required=True,
        type=Path,
        help="Reference genome in FASTA format",
    )
    indexer.add_argument(
        "-I",
        "--index",
        required=True,
        type=Path,
        help="Output index file path",
    )
    indexer.add_argument(
        "--kmer-size",
        type=int,
        default=15,
        help="K-mer size used to build the hash index (default: 15)",
    )
    indexer.add_argument(
        "--max-postings-per-kmer",
        type=int,
        default=0,
        help=(
            "Optional cap for number of stored positions per k-mer. "
            "0 means unlimited (default)."
        ),
    )
    indexer.add_argument(
        "--kmers-per-chunk",
        type=int,
        default=4000,
        help="Number of distinct k-mers per chunk file (default: 4000)",
    )

    return parser


def parse_args() -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "mapper":
        if args.reference is None or args.index is None:
            parser.error("mapper currently requires both -R/--reference and -I/--index")
        if args.max_errors < 0:
            parser.error("-k/--max-errors must be >= 0")
        if args.max_reads < 0:
            parser.error("--max-reads must be >= 0")
        if args.max_candidates <= 0:
            parser.error("--max-candidates must be > 0")

    if args.command == "indexer":
        if args.kmer_size <= 0:
            parser.error("--kmer-size must be > 0")
        if args.max_postings_per_kmer < 0:
            parser.error("--max-postings-per-kmer must be >= 0")
        if args.kmers_per_chunk <= 0:
            parser.error("--kmers-per-chunk must be > 0")

    return args


def _open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta_sequences(path: Path) -> List[Tuple[str, str]]:
    """Read FASTA(.gz) and return [(sequence_name, sequence)]."""
    sequences: List[Tuple[str, str]] = []
    current_name = None
    current_chunks: List[str] = []

    with _open_text_auto(path) as handle:
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


def read_fastq_reads(path: Path, max_reads: int = 0) -> Iterable[Tuple[str, str, str]]:
    """Yield FASTQ records as (read_name, sequence, quality)."""
    count = 0
    with _open_text_auto(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not seq or not plus or not qual:
                raise ValueError(f"Malformed FASTQ record near read {count + 1} in {path}")

            header = header.strip()
            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header at read {count + 1} in {path}")

            name = header[1:].split()[0]
            yield name, seq.strip().upper(), qual.strip()

            count += 1
            if max_reads > 0 and count >= max_reads:
                break


def build_kmer_index(
    references: Sequence[Tuple[str, str]],
    kmer_size: int,
    max_postings_per_kmer: int = 0,
) -> Dict[str, object]:
    """Build compact index: k-mer -> array('I') of global genome positions."""
    index: Dict[str, array] = {}
    total_bases = 0
    ref_names: List[str] = []
    ref_offsets: List[int] = []

    for ref_name, sequence in references:
        print(f"Indexing {ref_name} (length={len(sequence)})...")
        ref_names.append(ref_name)
        ref_offsets.append(total_bases)
        total_bases += len(sequence)
        if len(sequence) < kmer_size:
            continue

        ref_start = ref_offsets[-1]
        for i in range(0, len(sequence) - kmer_size + 1):
            # chromosome has around 87 to 150 million bases. Print each 5 million bases to track progress.
            if i % 5_000_000 == 0:
                print(f"  Processed {i} bases...")
            kmer = sequence[i : i + kmer_size]
            postings = index.setdefault(kmer, array("I"))

            if max_postings_per_kmer == 0 or len(postings) < max_postings_per_kmer:
                postings.append(ref_start + i)

    print(f"Indexing complete: {len(index)} distinct k-mers, total_bases={total_bases}")

    return {
        "format": "kmer_hash_u32_v2",
        "kmer_size": kmer_size,
        "total_kmers": len(index),
        "total_bases": total_bases,
        "ref_names": ref_names,
        "ref_offsets": ref_offsets,
        "index": index,
    }


def _chunk_sort_key(path: Path) -> int:
    parts = path.name.split(".")
    for i in range(len(parts) - 1):
        if parts[i].isdigit():
            return int(parts[i])
    return -1


def load_index(index_path: Path) -> Dict[str, object]:
    """Load single index file or multi-chunk index generated by save_index()."""
    if index_path.exists():
        opener = gzip.open if index_path.name.endswith(".gz") else open
        with opener(index_path, "rb") as handle:
            data = pickle.load(handle)
        if data.get("format") == "kmer_hash_u32_v2_chunk":
            # If a single chunk file was passed directly, still load all sibling chunks.
            pass
        else:
            return data

    name = index_path.name
    if name.endswith(".idx.gz"):
        prefix = name[: -len(".idx.gz")]
        chunk_pattern = f"{prefix}.*.idx.gz"
    elif name.endswith(".idx"):
        prefix = name[: -len(".idx")]
        chunk_pattern = f"{prefix}.*.idx"
    else:
        prefix = name
        chunk_pattern = f"{prefix}.*.idx*"

    chunk_files = sorted(index_path.parent.glob(chunk_pattern), key=_chunk_sort_key)
    if not chunk_files:
        raise FileNotFoundError(f"Could not find index file/chunks for base path: {index_path}")

    merged_index: Dict[str, array] = {}
    meta = None
    for chunk_file in chunk_files:
        opener = gzip.open if chunk_file.name.endswith(".gz") else open
        with opener(chunk_file, "rb") as handle:
            chunk_data = pickle.load(handle)

        if chunk_data.get("format") != "kmer_hash_u32_v2_chunk":
            raise ValueError(f"Unsupported chunk format in {chunk_file}")

        if meta is None:
            meta = {
                "format": "kmer_hash_u32_v2",
                "kmer_size": chunk_data["kmer_size"],
                "total_kmers": chunk_data["total_kmers"],
                "total_bases": chunk_data["total_bases"],
                "ref_names": chunk_data["ref_names"],
                "ref_offsets": chunk_data["ref_offsets"],
            }

        merged_index.update(chunk_data["index"])

    if meta is None:
        raise ValueError("No metadata found while loading index chunks")

    meta["index"] = merged_index
    return meta


def seed_candidates(
    read: str,
    index_data: Dict[str, object],
    max_errors: int,
    max_candidates: int,
) -> Deque[int]:
    """Stage 1: q-gram seeding returning candidate global start positions."""
    k = int(index_data["kmer_size"])
    if len(read) < k:
        return deque()

    votes: Counter = Counter()
    total_qgrams = len(read) - k + 1
    idx: Dict[str, array] = index_data["index"]

    for qpos in range(total_qgrams):
        qgram = read[qpos : qpos + k]
        postings = idx.get(qgram)
        if postings is None:
            continue
        for hit_pos in postings:
            start_pos = int(hit_pos) - qpos
            votes[start_pos] += 1

    if not votes:
        return deque()

    min_hits = max(1, total_qgrams // (max_errors + 1))
    selected = [pos for pos, hits in votes.items() if hits >= min_hits]
    if not selected:
        ranked = [pos for pos, _ in votes.most_common(max_candidates)]
        return deque(ranked)

    selected.sort(key=lambda p: votes[p], reverse=True)
    return deque(selected[:max_candidates])


def bounded_edit_distance(a: str, b: str, max_dist: int) -> int:
    """Levenshtein distance with early stop; returns max_dist+1 if exceeded."""
    n, m = len(a), len(b)
    if abs(n - m) > max_dist:
        return max_dist + 1

    prev = list(range(m + 1))
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        row_min = cur[0]
        ai = a[i - 1]

        for j in range(1, m + 1):
            cost = 0 if ai == b[j - 1] else 1
            cur[j] = min(
                prev[j] + 1,
                cur[j - 1] + 1,
                prev[j - 1] + cost,
            )
            if cur[j] < row_min:
                row_min = cur[j]

        if row_min > max_dist:
            return max_dist + 1
        prev = cur

    return prev[m]


def extract_reference_window(
    references: Sequence[Tuple[str, str]],
    ref_offsets: Sequence[int],
    global_start: int,
    read_len: int,
) -> Optional[Tuple[str, int, str]]:
    """Map global start to (ref_name, local_start, window_sequence)."""
    if global_start < 0:
        return None

    idx = bisect_right(ref_offsets, global_start) - 1
    if idx < 0 or idx >= len(references):
        return None

    ref_name, ref_seq = references[idx]
    local_start = global_start - ref_offsets[idx]
    local_end = local_start + read_len
    if local_start < 0 or local_end > len(ref_seq):
        return None

    return ref_name, local_start, ref_seq[local_start:local_end]


def map_read(
    read_seq: str,
    index_data: Dict[str, object],
    references: Sequence[Tuple[str, str]],
    max_errors: int,
    max_candidates: int,
) -> List[Tuple[int, str, int, str]]:
    """Return up to two best alignments as (dist, ref_name, local_start, cigar)."""
    candidates = seed_candidates(read_seq, index_data, max_errors, max_candidates)
    ref_offsets = index_data["ref_offsets"]

    hits: List[Tuple[int, str, int, str]] = []
    seen = set()
    while candidates:
        start = candidates.popleft()
        if start in seen:
            continue
        seen.add(start)

        ref_window = extract_reference_window(references, ref_offsets, start, len(read_seq))
        if ref_window is None:
            continue

        ref_name, local_start, ref_subseq = ref_window
        dist = bounded_edit_distance(read_seq, ref_subseq, max_errors)
        if dist <= max_errors:
            # PoC CIGAR: fixed-length match block for now.
            hits.append((dist, ref_name, local_start, f"{len(read_seq)}M"))

    hits.sort(key=lambda x: x[0])
    return hits[:2]


def write_mapper_output(
    output_path: Path,
    mapped_records: Iterable[Tuple[str, str, str, List[Tuple[int, str, int, str]]]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for read_name, read_seq, read_qual, hits in mapped_records:
            if not hits:
                handle.write(f"{read_name}\t*\t0\t*\t{read_seq}\t{read_qual}\n")
                continue

            best = hits[0]
            line = (
                f"{read_name}\t{best[1]}\t{best[2]}\t{best[3]}\t"
                f"{read_seq}\t{read_qual}"
            )
            if len(hits) > 1:
                alt = hits[1]
                line += f"\tALT:{alt[1]},{alt[2]},{alt[3]}"
            handle.write(line + "\n")


def _chunk_path(base_output_path: Path, chunk_id: int) -> Path:
    name = base_output_path.name
    if name.endswith(".idx.gz"):
        prefix = name[: -len(".idx.gz")]
        return base_output_path.with_name(f"{prefix}.{chunk_id}.idx.gz")
    if name.endswith(".idx"):
        prefix = name[: -len(".idx")]
        return base_output_path.with_name(f"{prefix}.{chunk_id}.idx")
    return base_output_path.with_name(f"{name}.{chunk_id}.idx")


def save_index(index_data: Dict[str, object], output_path: Path, kmers_per_chunk: int) -> None:
    print(f"Saving chunked index to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if output_path.name.endswith(".gz") else open

    all_items = list(index_data["index"].items())
    total_kmers = len(all_items)
    total_chunks = max(1, math.ceil(total_kmers / kmers_per_chunk))

    for chunk_id in range(total_chunks):
        start = chunk_id * kmers_per_chunk
        end = min(start + kmers_per_chunk, total_kmers)
        chunk_items = all_items[start:end]
        chunk_index = dict(chunk_items)

        chunk_data = {
            "format": "kmer_hash_u32_v2_chunk",
            "kmer_size": index_data["kmer_size"],
            "total_kmers": index_data["total_kmers"],
            "total_bases": index_data["total_bases"],
            "ref_names": index_data["ref_names"],
            "ref_offsets": index_data["ref_offsets"],
            "chunk_id": chunk_id,
            "chunk_count": total_chunks,
            "chunk_start_kmer": start,
            "chunk_end_kmer": end,
            "index": chunk_index,
        }

        chunk_output = _chunk_path(output_path, chunk_id)
        with opener(chunk_output, "wb") as handle:
            pickle.dump(chunk_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print(
            f"  Saved chunk {chunk_id + 1}/{total_chunks}: {chunk_output} "
            f"(kmers={len(chunk_index)}, size={chunk_output.stat().st_size} bytes)"
        )


def run_indexer(args: argparse.Namespace) -> None:
    references = read_fasta_sequences(args.reference)
    index_data = build_kmer_index(
        references=references,
        kmer_size=args.kmer_size,
        max_postings_per_kmer=args.max_postings_per_kmer,
    )
    save_index(index_data, args.index, args.kmers_per_chunk)

    print(
        f"Index built: {args.index} | k={index_data['kmer_size']} "
        f"| distinct_kmers={index_data['total_kmers']} | kmers_per_chunk={args.kmers_per_chunk}"
    )


def run_mapper(args: argparse.Namespace) -> None:
    print(f"Loading reference: {args.reference}")
    references = read_fasta_sequences(args.reference)
    print(f"Loading index: {args.index}")
    index_data = load_index(args.index)

    print(
        f"Mapper start | reads={args.input_reads} | k={args.max_errors} "
        f"| seed_k={index_data['kmer_size']}"
    )

    mapped: List[Tuple[str, str, str, List[Tuple[int, str, int, str]]]] = []
    for read_id, read_seq, read_qual in read_fastq_reads(args.input_reads, args.max_reads):
        hits = map_read(
            read_seq,
            index_data=index_data,
            references=references,
            max_errors=args.max_errors,
            max_candidates=args.max_candidates,
        )
        mapped.append((read_id, read_seq, read_qual, hits))

    write_mapper_output(args.output, mapped)
    total = len(mapped)
    aligned = sum(1 for _, _, _, hits in mapped if hits)
    print(f"Mapper finished | output={args.output} | reads={total} | aligned={aligned}")


def main() -> None:
    # run_indexer(argparse.Namespace(
    #     reference=Path("Data/chr16.fna"),
    #     index=Path("Data/chr16.idx.gz"),
    #     kmer_size=7,
    #     max_postings_per_kmer=0,
    #     kmers_per_chunk=4000,
    # ))
    args = parse_args()
    if args.command == "indexer":
        run_indexer(args)
        return

    if args.command == "mapper":
        run_mapper(args)
        return



if __name__ == "__main__":
    # Parsing PoC: next steps can consume this namespace for real pipeline work.
    main()