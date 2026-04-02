import importlib.util
import random
import statistics
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent / "Lab.S5"
sys.path.insert(0, str(BASE))


def load_module(name: str, path: Path):
  spec = importlib.util.spec_from_file_location(name, str(path))
  if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load module: {path}")
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


mod_full = load_module("ex3_1", BASE / "ex3_1.edit_alignment.py")
mod_lin = load_module("ex5_1", BASE / "ex5_1.lineal_mem.py")
mod_hir = load_module("ex5_2", BASE / "ex5_2.hirschberg.py")

RNG = random.Random(12345)
DNA = "ACGT"
SIZES = [768, 1024, 1280, 1536, 2048]
REPEATS = 5
OUTPUT_PLOT = Path(__file__).resolve().parent / "benchmark_alignment_scaling2.png"


def rand_seq(n: int) -> str:
  return "".join(RNG.choice(DNA) for _ in range(n))


def plot_growth(rows) -> None:
  grouped = {}
  for n, name, med, _, _, c_n2 in rows:
    if name not in grouped:
      grouped[name] = {"n": [], "time": [], "norm": []}
    grouped[name]["n"].append(n)
    grouped[name]["time"].append(med)
    grouped[name]["norm"].append(c_n2)

  fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
  for name, series in grouped.items():
    axes[0].plot(series["n"], series["time"], marker="o", linewidth=1.8, label=name)
    axes[1].plot(series["n"], series["norm"], marker="o", linewidth=1.8, label=name)

  axes[0].set_title("Median Runtime Growth")
  axes[0].set_xlabel("Sequence length n (m=n)")
  axes[0].set_ylabel("Median time (s)")
  axes[0].grid(alpha=0.3)

  axes[1].set_title("Quadratic Normalization Check")
  axes[1].set_xlabel("Sequence length n (m=n)")
  axes[1].set_ylabel("Median time / n^2")
  axes[1].grid(alpha=0.3)

  handles, labels = axes[0].get_legend_handles_labels()
  fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
  fig.tight_layout(rect=[0, 0, 1, 0.92])
  fig.savefig(OUTPUT_PLOT, dpi=180)
  plt.close(fig)


def run_benchmark() -> None:
  algorithms = [
    ("full_dp_score_only", lambda a, b: mod_full.compute_dp_matrix(a, b)[-1][-1]),
    ("linear_space_score", lambda a, b: mod_lin.global_alignment_linear_space(a, b)),
    ("hirschberg_align", lambda a, b: mod_hir.global_alignment_hirschberg(a, b)[0]),
  ]

  rows = []
  print("size,algo,median_s,mean_s,std_s,time_per_n2")

  for n in SIZES:
    x = rand_seq(n)
    y = rand_seq(n)

    for name, fn in algorithms:
      samples = []
      for _ in range(REPEATS):
        t0 = time.perf_counter()
        fn(x, y)
        t1 = time.perf_counter()
        samples.append(t1 - t0)

      med = statistics.median(samples)
      mean = statistics.mean(samples)
      std = statistics.pstdev(samples)
      c_n2 = med / (n * n)
      rows.append((n, name, med, mean, std, c_n2))

      print(f"{n},{name},{med:.6f},{mean:.6f},{std:.6f},{c_n2:.10f}")

  print("\nsummary: average time_per_n2 by algorithm")
  for name in ["full_dp_score_only", "linear_space_score", "hirschberg_align"]:
    values = [r[5] for r in rows if r[1] == name]
    print(f"{name},{statistics.mean(values):.10f}")

  plot_growth(rows)
  print(f"\nSaved growth plot to: {OUTPUT_PLOT}")


if __name__ == "__main__":
  run_benchmark()
