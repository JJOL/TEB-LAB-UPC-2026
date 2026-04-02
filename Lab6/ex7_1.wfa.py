
"""
Edit-distance Wavefront Alignment Algorithm (WFA) with adaptive visualization.

Key improvements over the original pedagogical version:
  - compact array-based wavefront representation (lower Python overhead)
  - faster extend loop using local variables and fewer dictionary lookups
  - adaptive plotting defaults for large sequence pairs
  - optional heavy outputs (wavefront dumps and DP table printing)
"""

import argparse
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from wfa_plot import build_wfa_dp_matrix, plot_wfa_dp_matrix


@dataclass
class WavefrontLevel:
  score: int
  lo: int
  hi: int
  offsets: List[int]
  starts: Optional[List[int]] = None


def wavefront_v(k: int, offset: int) -> int:
  return offset - k


def wavefront_h(offset: int) -> int:
  return offset


def diagonal(h: int, v: int) -> int:
  return h - v


def offset_of(h: int, v: int) -> int:
  del v
  return h


def get_offset(level: WavefrontLevel, k: int) -> int:
  idx = k - level.lo
  if 0 <= idx < len(level.offsets):
    return level.offsets[idx]
  return -1


def normalize_sequence(sequence: str) -> str:
  # Support multiline FASTA/plain sequence content.
  return "".join(sequence.split())


def read_sequence_from_file(path: str) -> str:
  with open(path, "r", encoding="utf-8") as handle:
    lines = handle.readlines()

  if len(lines) > 0 and lines[0].startswith(">"):
    body = "".join(line.strip() for line in lines if not line.startswith(">"))
    return normalize_sequence(body)

  return normalize_sequence("".join(lines))


def extend_wavefront(
  level: WavefrontLevel,
  pattern: str,
  text: str,
  keep_starts: bool
) -> WavefrontLevel:
  lo = level.lo
  old_offsets = level.offsets
  n = len(pattern)
  m = len(text)

  extended = [-1] * len(old_offsets)
  starts = [-1] * len(old_offsets) if keep_starts else None

  for idx, start in enumerate(old_offsets):
    if start < 0:
      continue

    k = lo + idx
    current = start
    if keep_starts:
      starts[idx] = start

    while True:
      h = current
      v = h - k
      if h >= m or v >= n:
        break
      if pattern[v] != text[h]:
        break
      current += 1

    extended[idx] = current

  return WavefrontLevel(
    score=level.score,
    lo=lo,
    hi=level.hi,
    offsets=extended,
    starts=starts
  )


def compute_next_wavefront(level: WavefrontLevel) -> WavefrontLevel:
  prev = level.offsets
  old_lo = level.lo
  old_hi = level.hi

  new_lo = old_lo - 1
  new_hi = old_hi + 1
  next_offsets = [-1] * (len(prev) + 2)

  for k in range(new_lo, new_hi + 1):
    left = get_offset(level, k - 1)
    mid = get_offset(level, k)
    right = get_offset(level, k + 1)

    insertion = left + 1 if left >= 0 else -1
    substitution = mid + 1 if mid >= 0 else -1
    deletion = right if right >= 0 else -1

    next_offsets[k - new_lo] = max(insertion, substitution, deletion)

  return WavefrontLevel(
    score=level.score + 1,
    lo=new_lo,
    hi=new_hi,
    offsets=next_offsets,
    starts=None
  )


def level_to_dict(level: WavefrontLevel) -> Dict[int, int]:
  out = {}
  for idx, value in enumerate(level.offsets):
    if value >= 0:
      out[level.lo + idx] = value
  return out


def level_starts_to_dict(level: WavefrontLevel) -> Dict[int, int]:
  if level.starts is None:
    return {}

  out = {}
  for idx, value in enumerate(level.starts):
    if value >= 0:
      out[level.lo + idx] = value
  return out


def print_wavefronts(wavefronts: List[WavefrontLevel], pattern: str, text: str) -> None:
  del pattern, text
  print("Wavefronts:")
  print()

  for level in wavefronts:
    score = level.score
    offsets = level_to_dict(level)
    starts = level_starts_to_dict(level)
    keys = sorted(offsets.keys())

    print(f"s = {score}")
    for k in keys:
      start = starts.get(k, offsets[k])
      end = offsets[k]
      v = wavefront_v(k, end)
      h = wavefront_h(end)
      print(
        f"  k = {k:3d}   start = {start:3d}   end = {end:3d}   cell = ({v:3d}, {h:3d})"
      )
    print()


def backtrace_edit_wavefronts(
  pattern: str,
  text: str,
  wavefronts: List[WavefrontLevel]
) -> Optional[str]:
  target_k = diagonal(len(text), len(pattern))
  score = len(wavefronts) - 1
  offsets = wavefronts[score]

  if get_offset(offsets, target_k) < 0:
    return None

  offset = get_offset(offsets, target_k)
  cigar_rev = []
  k = target_k

  while score > 0:
    prev = wavefronts[score - 1]

    if get_offset(prev, k + 1) >= 0 and offset == get_offset(prev, k + 1):
      cigar_rev.append("D")
      k += 1
      score -= 1
    elif get_offset(prev, k - 1) >= 0 and offset == get_offset(prev, k - 1) + 1:
      cigar_rev.append("I")
      k -= 1
      offset -= 1
      score -= 1
    elif get_offset(prev, k) >= 0 and offset == get_offset(prev, k) + 1:
      cigar_rev.append("X")
      offset -= 1
      score -= 1
    else:
      cigar_rev.append("M")
      offset -= 1

  while offset > 0:
    cigar_rev.append("M")
    offset -= 1

  cigar_rev.reverse()
  return "".join(cigar_rev)


def cigar_to_alignment(pattern: str, text: str, cigar: str) -> Tuple[str, str, str]:
  aligned_pattern = []
  aligned_text = []
  middle = []

  i = 0
  j = 0

  for op in cigar:
    if op == "M":
      aligned_pattern.append(pattern[i])
      aligned_text.append(text[j])
      middle.append("|")
      i += 1
      j += 1
    elif op == "X":
      aligned_pattern.append(pattern[i])
      aligned_text.append(text[j])
      middle.append(" ")
      i += 1
      j += 1
    elif op == "D":
      aligned_pattern.append(pattern[i])
      aligned_text.append("-")
      middle.append(" ")
      i += 1
    elif op == "I":
      aligned_pattern.append("-")
      aligned_text.append(text[j])
      middle.append(" ")
      j += 1

  return "".join(aligned_pattern), "".join(middle), "".join(aligned_text)


def cigar_to_path(cigar: str) -> List[Tuple[int, int]]:
  path = []
  i = 0
  j = 0
  path.append((i, j))

  for op in cigar:
    if op == "M" or op == "X":
      i += 1
      j += 1
    elif op == "D":
      i += 1
    elif op == "I":
      j += 1
    path.append((i, j))

  return path


def align_edit_wavefronts(
  pattern: str,
  text: str,
  keep_starts_for_plot: bool = True,
  max_score: Optional[int] = None
) -> Optional[dict]:
  pattern = normalize_sequence(pattern)
  text = normalize_sequence(text)

  target_k = diagonal(len(text), len(pattern))
  target_offset = offset_of(len(text), len(pattern))
  max_distance = len(pattern) + len(text)
  if max_score is not None:
    max_distance = min(max_distance, max_score)

  current = WavefrontLevel(score=0, lo=0, hi=0, offsets=[0], starts=None)
  wavefronts: List[WavefrontLevel] = []

  while current.score <= max_distance:
    extended = extend_wavefront(
      current,
      pattern,
      text,
      keep_starts=keep_starts_for_plot
    )
    wavefronts.append(extended)

    if get_offset(extended, target_k) == target_offset:
      break

    current = compute_next_wavefront(extended)

  if len(wavefronts) == 0 or get_offset(wavefronts[-1], target_k) != target_offset:
    return None

  cigar = backtrace_edit_wavefronts(pattern, text, wavefronts)
  if cigar is None:
    return None

  aligned_pattern, middle, aligned_text = cigar_to_alignment(pattern, text, cigar)
  path = cigar_to_path(cigar)
  wavefront_elements = sum(level.hi - level.lo + 1 for level in wavefronts)

  result = {
    "distance": len(wavefronts) - 1,
    "cigar": cigar,
    "aligned_pattern": aligned_pattern,
    "middle": middle,
    "aligned_text": aligned_text,
    "path": path,
    "wavefronts": wavefronts,
    "wavefront_elements": wavefront_elements,
    "pattern": pattern,
    "text": text
  }
  return result


def compact_wavefronts_to_plot_wavefronts(
  wavefronts: List[WavefrontLevel],
  include_starts: bool
) -> List[dict]:
  out = []
  for level in wavefronts:
    offsets = level_to_dict(level)
    if include_starts and level.starts is not None:
      starts = level_starts_to_dict(level)
    else:
      # Fallback: plot only the frontier if segment starts were not stored.
      starts = dict(offsets)

    out.append(
      {
        "score": level.score,
        "offsets": offsets,
        "starts": starts
      }
    )
  return out


def maybe_build_equivalent_dp(
  result: dict,
  max_plot_cells: int,
  force_frontier_only: bool = False
) -> None:
  pattern = result["pattern"]
  text = result["text"]
  n = len(pattern) + 1
  m = len(text) + 1
  total_cells = n * m

  if force_frontier_only or total_cells > max_plot_cells:
    result["dp"] = None
    result["computed_dp_cells"] = 0
    result["plot_wavefronts"] = compact_wavefronts_to_plot_wavefronts(
      result["wavefronts"],
      include_starts=False
    )
    result["frontier_only_plot"] = True
    return

  plot_wavefronts = compact_wavefronts_to_plot_wavefronts(
    result["wavefronts"],
    include_starts=not force_frontier_only
  )
  dp, visited = build_wfa_dp_matrix(pattern, text, plot_wavefronts)
  result["dp"] = dp
  result["computed_dp_cells"] = len(visited)
  result["plot_wavefronts"] = plot_wavefronts
  result["frontier_only_plot"] = force_frontier_only


def print_alignment(aligned_pattern: str, middle: str, aligned_text: str) -> None:
  print("Alignment:")
  print(aligned_pattern)
  print(middle)
  print(aligned_text)
  print()


def print_dp_table(dp: List[List[Optional[int]]]) -> None:
  print("Equivalent DP matrix:")
  i = 0
  while i < len(dp):
    pieces = []
    j = 0
    while j < len(dp[0]):
      value = dp[i][j]
      if value is None:
        pieces.append("  .")
      else:
        pieces.append(f"{value:3d}")
      j += 1
    print("".join(pieces))
    i += 1
  print()


def run_example(
  pattern: str,
  text: str,
  show_plot: bool = True,
  print_all_wavefronts: bool = False,
  print_dp: bool = False,
  max_plot_cells: int = 120_000
) -> Optional[dict]:
  normalized_pattern = normalize_sequence(pattern)
  normalized_text = normalize_sequence(text)
  total_cells = (len(normalized_pattern) + 1) * (len(normalized_text) + 1)

  keep_starts_for_plot = show_plot and total_cells <= max_plot_cells

  result = align_edit_wavefronts(
    normalized_pattern,
    normalized_text,
    keep_starts_for_plot=keep_starts_for_plot
  )
  if result is None:
    print("No alignment found")
    return None

  maybe_build_equivalent_dp(
    result,
    max_plot_cells=max_plot_cells,
    force_frontier_only=(show_plot and not keep_starts_for_plot)
  )

  print("Pattern length:", len(result["pattern"]))
  print("Text length:   ", len(result["text"]))
  print()
  print("Edit distance:", result["distance"])
  print("CIGAR:", result["cigar"])
  print("Wavefront elements:", result["wavefront_elements"])
  if result["dp"] is None:
    print("Computed DP cells: skipped (matrix too large, plotting frontier only)")
  else:
    print("Computed DP cells:", result["computed_dp_cells"])
  print()

  if print_all_wavefronts:
    print_wavefronts(result["wavefronts"], pattern, text)

  print_alignment(
    result["aligned_pattern"],
    result["middle"],
    result["aligned_text"]
  )

  if print_dp and result["dp"] is not None:
    print_dp_table(result["dp"])

  if show_plot:
    annotate = total_cells <= 2500
    title = f"WFA Equivalent DP Matrix (distance={result['distance']})"

    if result["dp"] is None:
      dp, _ = build_wfa_dp_matrix(
        result["pattern"],
        result["text"],
        result["plot_wavefronts"]
      )
      plot_wfa_dp_matrix(
        dp,
        result["pattern"],
        result["text"],
        path=result["path"],
        title=title + " [frontier]",
        annotate=False,
        missing_label="",
        show_colorbar=True
      )
    else:
      plot_wfa_dp_matrix(
        result["dp"],
        result["pattern"],
        result["text"],
        path=result["path"],
        title=title,
        annotate=annotate,
        missing_label="",
        show_colorbar=True
      )

  return result


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(
    description="Edit-distance WFA with adaptive large-sequence visualization."
  )
  parser.add_argument("--pattern", type=str, help="Pattern sequence")
  parser.add_argument("--text", type=str, help="Text sequence")
  parser.add_argument("--pattern-file", type=str, help="Pattern file (FASTA/plain)")
  parser.add_argument("--text-file", type=str, help="Text file (FASTA/plain)")
  parser.add_argument("--no-plot", action="store_true", help="Disable plotting")
  parser.add_argument(
    "--print-wavefronts",
    action="store_true",
    help="Print all wavefront levels"
  )
  parser.add_argument(
    "--print-dp",
    action="store_true",
    help="Print equivalent DP matrix (small inputs only)"
  )
  parser.add_argument(
    "--max-plot-cells",
    type=int,
    default=120000,
    help="Maximum DP cells to fully materialize before switching to frontier-only"
  )
  return parser.parse_args()


def resolve_input_sequences(args: argparse.Namespace) -> Tuple[str, str]:
  if args.pattern_file:
    pattern = read_sequence_from_file(args.pattern_file)
  else:
    pattern = args.pattern if args.pattern is not None else "ACTATTTACGTACT"

  if args.text_file:
    text = read_sequence_from_file(args.text_file)
  else:
    text = args.text if args.text is not None else "ACTATTACCTACT"

  return pattern, text


if __name__ == "__main__":
  args = parse_args()
  pattern, text = resolve_input_sequences(args)
  run_example(
    pattern,
    text,
    show_plot=not args.no_plot,
    print_all_wavefronts=args.print_wavefronts,
    print_dp=args.print_dp,
    max_plot_cells=args.max_plot_cells
  )


# What is the research gap
# What is the main proposal
# How does it compare to previous work (state of the art)