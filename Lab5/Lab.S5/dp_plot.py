"""
DP matrix visualization for edit distance.

Given a DP matrix, draw:
  - a heatmap of matrix values
  - the numeric value inside each cell
  - optionally, an optimal backtrace path

Time:
  - plotting the matrix values: O(n * m)
  - plotting the path: O(n + m)

Space:
  - O(1) extra, excluding matplotlib internal structures
"""

import math
from typing import List, Tuple, Optional

import matplotlib.pyplot as plt

def plot_dp_matrix(
  dp: List[List[int]],
  x: str,
  y: str,
  path: Optional[List[Tuple[int, int]]] = None,
  title: str = "Edit Distance DP Matrix",
  annotate: bool = True,
  max_tick_labels: int = 60
) -> None:
  """
  Plot a DP matrix as a heatmap and optionally overlay a backtrace path.

  Parameters:
    dp: DP matrix of size (len(x) + 1) x (len(y) + 1)
    x: first sequence
    y: second sequence
    path: optional list of cells (i, j) from (0, 0) to (n, m)
    title: plot title
    annotate: whether to print numeric values inside cells

  Time: O(n * m)
  Space: O(1) extra
  """
  n = len(dp)
  m = len(dp[0])
  total_cells = n * m

  # Keep the image readable and avoid creating extremely large figures.
  max_fig_w = 16
  max_fig_h = 12
  cell_size = 0.28

  fig_w = min(max_fig_w, max(6, cell_size * m + 2))
  fig_h = min(max_fig_h, max(5, cell_size * n + 2))
  fig, ax = plt.subplots(figsize=(fig_w, fig_h))

  image = ax.imshow(dp, cmap="RdYlGn", aspect="equal", origin="upper")
  cbar = plt.colorbar(
    image,
    ax=ax,
    shrink=0.7,     # reduce altura
    aspect=25,      # hace la barra más fina
    pad=0.02        # separación del gráfico
  )
  cbar.set_label("DP score", fontsize=10)
  cbar.ax.tick_params(labelsize=8)

  # Axis labels
  x_labels = ["-"] + list(y)
  y_labels = ["-"] + list(x)

  x_step = max(1, math.ceil(m / max_tick_labels))
  y_step = max(1, math.ceil(n / max_tick_labels))
  xticks = list(range(0, m, x_step))
  yticks = list(range(0, n, y_step))

  if xticks[-1] != m - 1:
    xticks.append(m - 1)
  if yticks[-1] != n - 1:
    yticks.append(n - 1)

  ax.set_xticks(xticks)
  ax.set_yticks(yticks)
  ax.set_xticklabels([x_labels[j] for j in xticks])
  ax.set_yticklabels([y_labels[i] for i in yticks])

  ax.set_xlabel("Sequence X")
  ax.set_ylabel("Sequence Y")
  ax.set_title(title)

  # Grid lines are useful for small matrices but expensive for large ones.
  show_grid = total_cells <= 8000
  if show_grid:
    ax.set_xticks([j - 0.5 for j in range(1, m)], minor=True)
    ax.set_yticks([i - 0.5 for i in range(1, n)], minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

  # Cell annotations
  effective_annotate = annotate and total_cells <= 2500
  if effective_annotate:
    max_value = max(max(row) for row in dp)
    threshold = max_value / 2 if max_value > 0 else 0

    for i in range(n):
      for j in range(m):
        text_color = "white" if dp[i][j] > threshold else "black"
        ax.text(
          j,
          i,
          str(dp[i][j]),
          ha="center",
          va="center",
          color=text_color,
          fontsize=9
        )

  # Overlay optimal path
  if path is not None and len(path) > 0:
    xs = [j for i, j in path]
    ys = [i for i, j in path]

    if len(path) <= 400:
      ax.plot(xs, ys, marker="o", linewidth=2)
    else:
      ax.plot(xs, ys, linewidth=1.6)

  plt.tight_layout()
  plt.show()