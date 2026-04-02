
"""
Banded global alignment under edit distance.

This version reuses the external plotting module dp_plot.py.

It:
  - computes a banded dynamic programming table
  - counts how many cells were actually computed
  - computes the edit distance if the optimal path fits in the band
  - backtraces one optimal alignment
  - prints the alignment using a gap representation
  - plots the partially computed DP table using plot_dp_matrix()
"""

from dp_plot import plot_dp_matrix

def in_band(i, j, band):
  return abs(i - j) <= band


def banded_edit_distance(x, y, band):
  n = len(x)
  m = len(y)

  dp = []
  trace = []
  for i in range(n + 1):
    dp_row = []
    trace_row = []
    for j in range(m + 1):
      dp_row.append(None)
      trace_row.append(None)
    dp.append(dp_row)
    trace.append(trace_row)

  computed_cells = 0

  for i in range(n + 1):
    if in_band(i, 0, band):
      dp[i][0] = i
      computed_cells += 1
      if i > 0:
        trace[i][0] = "U"

  for j in range(m + 1):
    if in_band(0, j, band):
      if dp[0][j] is None:
        computed_cells += 1
      dp[0][j] = j
      if j > 0:
        trace[0][j] = "L"

  for i in range(1, n + 1):
    j_start = max(1, i - band)
    j_end = min(m, i + band)

    for j in range(j_start, j_end + 1):
      best = None
      move = None

      if dp[i - 1][j] is not None:
        value = dp[i - 1][j] + 1
        best = value
        move = "U"

      if dp[i][j - 1] is not None:
        value = dp[i][j - 1] + 1
        if best is None or value < best:
          best = value
          move = "L"

      if dp[i - 1][j - 1] is not None:
        mismatch = 0
        if x[i - 1] != y[j - 1]:
          mismatch = 1
        value = dp[i - 1][j - 1] + mismatch
        if best is None or value < best:
          best = value
          move = "D"

      if best is not None:
        dp[i][j] = best
        trace[i][j] = move
        computed_cells += 1

  return dp, trace, dp[n][m], computed_cells


def backtrace_alignment(x, y, dp, trace):
  i = len(x)
  j = len(y)

  if dp[i][j] is None:
    return None, None, None, None

  aligned_x = []
  aligned_y = []
  path = []

  while True:
    path.append((i, j))
    if i == 0 and j == 0:
      break

    move = trace[i][j]

    if move == "D":
      aligned_x.append(x[i - 1])
      aligned_y.append(y[j - 1])
      i -= 1
      j -= 1
    elif move == "U":
      aligned_x.append(x[i - 1])
      aligned_y.append("-")
      i -= 1
    elif move == "L":
      aligned_x.append("-")
      aligned_y.append(y[j - 1])
      j -= 1
    else:
      return None, None, None, None

  aligned_x.reverse()
  aligned_y.reverse()
  path.reverse()

  middle = []
  for k in range(len(aligned_x)):
    if aligned_x[k] == aligned_y[k]:
      middle.append("|")
    else:
      middle.append(" ")

  return "".join(aligned_x), "".join(middle), "".join(aligned_y), path


def print_alignment(aligned_x, middle, aligned_y):
  print(aligned_x)
  print(middle)
  print(aligned_y)


def print_dp_table(dp):
  print("DP table:")
  for row in dp:
    pieces = []
    for value in row:
      if value is None:
        pieces.append("  .")
      else:
        pieces.append(f"{value:3d}")
    print("".join(pieces))


def run_example(x, y, band):
  dp, trace, distance, computed_cells = banded_edit_distance(x, y, band)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print("Band:", band)
  print("Matrix size:", (len(x) + 1), "x", (len(y) + 1))
  print("Computed cells:", computed_cells)

  total_cells = (len(x) + 1) * (len(y) + 1)
  print("Total cells:", total_cells)

  if distance is None:
    print("Distance: not available inside this band")
    print_dp_table(dp)
    plot_dp_matrix(
      dp,
      x,
      y,
      title=f"Banded DP Table (band={band}, no full solution)",
      annotate=True,
      missing_value=None,
      missing_label=""
    )
    return

  print("Edit distance:", distance)

  aligned_x, middle, aligned_y, path = backtrace_alignment(x, y, dp, trace)

  print()
  print("Alignment:")
  print_alignment(aligned_x, middle, aligned_y)
  print()
  print_dp_table(dp)

  plot_dp_matrix(
    dp,
    x,
    y,
    path=path,
    title=f"Banded Edit Distance DP Table (band={band}, distance={distance})",
    annotate=True,
    missing_value=None,
    missing_label=""
  )


if __name__ == "__main__":
  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"
  band = 2
  run_example(x, y, band)

  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"
  band = 6
  run_example(x, y, band)

