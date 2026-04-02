"""Benchmark script: measure how indexer + mapper time and memory scale with genome size.

Each configuration generates a synthetic reference of increasing size, builds the
k-mer index and maps a fixed number of reads.  Results are collected and plotted as:
  - genome size vs indexer wall time + peak memory
  - genome size vs mapper  wall time + peak memory

Usage:
    python test_growth.py                  # run all configs and plot
    python test_growth.py --plot-only      # re-plot from saved results.json
    python test_growth.py --configs 3      # run only first N configs (quick test)
"""

import argparse
import importlib.util
import json
import os
import sys
import time
import tracemalloc
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Locate example.py and generate_test.py next to this file
# ---------------------------------------------------------------------------
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import generate_test as gen

spec = importlib.util.spec_from_file_location("example", HERE / "example.py")
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)


# ---------------------------------------------------------------------------
# Benchmark configurations
# ---------------------------------------------------------------------------
CONFIGS: List[Dict[str, Any]] = [
    {"genome_len": 10_000,    "n_reads": 20,  "read_len": 150, "kmer_size": 11},
    {"genome_len": 50_000,    "n_reads": 50,  "read_len": 150, "kmer_size": 11},
    {"genome_len": 100_000,   "n_reads": 50,  "read_len": 150, "kmer_size": 11},
    {"genome_len": 200_000,   "n_reads": 50,  "read_len": 150, "kmer_size": 11},
    {"genome_len": 500_000,   "n_reads": 100, "read_len": 150, "kmer_size": 13},
    {"genome_len": 1_000_000, "n_reads": 100, "read_len": 150, "kmer_size": 13},
    {"genome_len": 2_000_000, "n_reads": 200, "read_len": 150, "kmer_size": 15},
    {"genome_len": 5_000_000, "n_reads": 200, "read_len": 150, "kmer_size": 15},
]

WORK_DIR = Path("/tmp/mapper_benchmark")
ERROR_RATE = 0.02   # 2% -> ~3 errors in 150bp reads
MAX_ERRORS = 6      # generous k for benchmark correctness
MAX_CANDIDATES = 2000
SEED = 42
RESULTS_JSON = WORK_DIR / "results.json"
PLOT_PATH = WORK_DIR / "benchmark.png"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mb(bytes_: int) -> float:
    return bytes_ / (1024 ** 2)


def run_and_measure(fn, *args, **kwargs):
    """Run fn(*args, **kwargs) and return (result, wall_seconds, peak_mb)."""
    tracemalloc.start()
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    elapsed = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, elapsed, _mb(peak_bytes)


def _fake_namespace(**kwargs):
    import argparse
    ns = argparse.Namespace()
    for k, v in kwargs.items():
        setattr(ns, k, v)
    return ns


# ---------------------------------------------------------------------------
# Per-configuration benchmark
# ---------------------------------------------------------------------------

def run_config(cfg: Dict[str, Any], idx: int) -> Dict[str, Any]:
    genome_len = cfg["genome_len"]
    n_reads    = cfg["n_reads"]
    read_len   = cfg["read_len"]
    kmer_size  = cfg["kmer_size"]

    label = f"{genome_len // 1000}K"
    print(f"\n[{idx+1}/{len(CONFIGS)}] genome={label}  reads={n_reads}  k={kmer_size}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    ref_path   = WORK_DIR / f"ref_{genome_len}.fna"
    idx_path   = WORK_DIR / f"ref_{genome_len}.idx"
    reads_path = WORK_DIR / f"reads_{genome_len}.fastq"
    out_path   = WORK_DIR / f"out_{genome_len}.sam"

    # ---- generate data ----
    import random
    random.seed(SEED)
    gen_ref_ns = _fake_namespace(
        output=ref_path, num_sequences=1, length=genome_len,
        line_width=80,
    )
    gen.generate_reference(gen_ref_ns)

    random.seed(SEED + 1)
    gen_reads_ns = _fake_namespace(
        reference=ref_path, output=reads_path,
        num_reads=n_reads, read_length=read_len,
        error_rate=ERROR_RATE, indel_fraction=0.2,
    )
    gen.generate_reads(gen_reads_ns)

    # ---- indexer ----
    # Remove stale chunk files
    for stale in WORK_DIR.glob(f"ref_{genome_len}.*.idx"):
        stale.unlink()

    idx_ns = _fake_namespace(
        reference=ref_path, index=idx_path,
        kmer_size=kmer_size, max_postings_per_kmer=0,
        kmers_per_chunk=50_000,
    )
    print(f"  Indexing...", end=" ", flush=True)
    _, idx_time, idx_mem = run_and_measure(ex.run_indexer, idx_ns)
    print(f"  {idx_time:.2f}s  {idx_mem:.1f} MB")

    # ---- mapper ----
    references = ex.read_fasta_sequences(ref_path)
    index_data = ex.load_index(idx_path)

    map_ns = _fake_namespace(
        reference=ref_path, index=idx_path, input_reads=reads_path,
        output=out_path, max_errors=MAX_ERRORS,
        max_reads=0, max_candidates=MAX_CANDIDATES,
    )
    print(f"  Mapping...", end=" ", flush=True)
    _, map_time, map_mem = run_and_measure(ex.run_mapper, map_ns)
    print(f"  {map_time:.2f}s  {map_mem:.1f} MB")

    # ---- count aligned ----
    aligned = 0
    if out_path.exists():
        with out_path.open() as f:
            for line in f:
                if "\t*\t" not in line.split("\t")[1:2]:
                    aligned += 1

    return {
        "genome_len":   genome_len,
        "n_reads":      n_reads,
        "kmer_size":    kmer_size,
        "idx_time_s":   round(idx_time, 4),
        "idx_mem_mb":   round(idx_mem, 2),
        "map_time_s":   round(map_time, 4),
        "map_mem_mb":   round(map_mem, 2),
        "aligned":      aligned,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(results: List[Dict[str, Any]], out_path: Path) -> None:
    sizes  = [r["genome_len"] / 1_000 for r in results]  # in Kbp
    labels = [f"{int(s)}K" for s in sizes]

    idx_times = [r["idx_time_s"] for r in results]
    map_times = [r["map_time_s"] for r in results]
    idx_mems  = [r["idx_mem_mb"] for r in results]
    map_mems  = [r["map_mem_mb"] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Indexer + Mapper scaling vs genome size", fontsize=13)

    # ---- time plot ----
    ax = axes[0]
    ax.plot(sizes, idx_times, "o-", color="steelblue",  label="indexer")
    ax.plot(sizes, map_times, "s-", color="darkorange", label="mapper")
    ax.set_xlabel("Genome size (Kbp)")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Time scaling")
    ax.set_xticks(sizes)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)

    # ---- memory plot ----
    ax = axes[1]
    ax.plot(sizes, idx_mems, "o-", color="steelblue",  label="indexer")
    ax.plot(sizes, map_mems, "s-", color="darkorange", label="mapper")
    ax.set_xlabel("Genome size (Kbp)")
    ax.set_ylabel("Peak memory (MB)")
    ax.set_title("Memory scaling")
    ax.set_xticks(sizes)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    print(f"\nPlot saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Mapper scaling benchmark")
    p.add_argument("--plot-only", action="store_true",
                   help="Skip benchmarking, re-plot from saved results.json")
    p.add_argument("--configs", type=int, default=0,
                   help="Run only first N configurations (0 = all)")
    return p.parse_args()


def main():
    args = parse_args()

    if args.plot_only:
        if not RESULTS_JSON.exists():
            sys.exit(f"No results file found at {RESULTS_JSON}")
        with RESULTS_JSON.open() as f:
            results = json.load(f)
        plot_results(results, PLOT_PATH)
        return

    configs = CONFIGS[: args.configs] if args.configs > 0 else CONFIGS
    results = []
    for i, cfg in enumerate(configs):
        record = run_config(cfg, i)
        results.append(record)
        # Save incrementally so a crash doesn't lose earlier data
        WORK_DIR.mkdir(parents=True, exist_ok=True)
        with RESULTS_JSON.open("w") as f:
            json.dump(results, f, indent=2)

    print("\n=== Summary ===")
    print(f"{'Genome':>10}  {'idx_t(s)':>9}  {'idx_mem(MB)':>11}  "
          f"{'map_t(s)':>9}  {'map_mem(MB)':>11}  {'aligned':>7}")
    for r in results:
        print(f"{r['genome_len']:>10,}  {r['idx_time_s']:>9.2f}  "
              f"{r['idx_mem_mb']:>11.1f}  {r['map_time_s']:>9.2f}  "
              f"{r['map_mem_mb']:>11.1f}  {r['aligned']:>7}")

    plot_results(results, PLOT_PATH)


if __name__ == "__main__":
    main()
