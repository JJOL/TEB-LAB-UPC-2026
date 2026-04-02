"""
DP matrix visualization for edit distance and related alignment variants.

Given a DP matrix, draw:
  - a heatmap of matrix values
  - the numeric value inside each computed cell
  - optionally, an optimal backtrace path
  - optionally, a mask for cells that have not been computed yet

This version is designed to work well with:
  - full DP tables
  - partially computed DP tables
  - banded alignment tables
  - custom plot titles

Convention for non-computed cells:
  - use None, or
  - use a configurable sentinel value such as -1

Time:
  - plotting the matrix values: O(n * m)
  - plotting the path: O(len(path))

Space:
  - O(n * m) for the plotting copy/mask
"""

import math
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def _is_missing(value, missing_value):
  if value is None:
    return True
  if missing_value is None:
    return False
  return value == missing_value


def _build_plot_matrix(dp, missing_value):
  n = len(dp)
  m = len(dp[0])

  values = []
  missing_mask = []
  computed_values = []

  for i in range(n):
    value_row = []
    mask_row = []
    for j in range(m):
      cell = dp[i][j]
      missing = _is_missing(cell, missing_value)
      mask_row.append(missing)
      if missing:
        value_row.append(float("nan"))
      else:
        numeric = float(cell)
        value_row.append(numeric)
        computed_values.append(numeric)
      values.append(value_row) if False else None
    values.append(value_row)
    missing_mask.append(mask_row)

  return values, missing_mask, computed_values


def plot_dp_matrix(
  dp,
  x,
  y,
  path=None,
  title="DP Matrix",
  annotate=True,
  missing_value=None,
  missing_label="",
  show_colorbar=True,
  cmap="YlOrRd",
  computed_alpha=1.0,
  missing_color="#d9d9d9",
  path_style="-o"
):
  """
  Plot a DP matrix as a heatmap and optionally overlay a backtrace path.

  Parameters:
    dp:
      Matrix of size (len(x) + 1) x (len(y) + 1).
      Non-computed cells may be represented by None or by missing_value.
    x:
      First sequence (rows).
    y:
      Second sequence (columns).
    path:
      Optional list of cells (i, j) to overlay.
    title:
      Plot title.
    annotate:
      Whether to print values inside computed cells.
    missing_value:
      Sentinel used for non-computed cells, e.g. -1.
      If None, only cells equal to None are treated as missing.
    missing_label:
      Text shown inside missing cells if annotate=True.
      Typical choices: "", ".", "?", "X".
    show_colorbar:
      Whether to display the colorbar.
    cmap:
      Colormap used for computed cells.
    computed_alpha:
      Opacity of computed cells.
    missing_color:
      Color used for non-computed cells.
    path_style:
      Matplotlib line style for the optional path.

  Time: O(n * m)
  Space: O(n * m)
  """
  if dp is None or len(dp) == 0 or len(dp[0]) == 0:
    raise ValueError("dp must be a non-empty rectangular matrix")

  n = len(dp)
  m = len(dp[0])

  for row in dp:
    if len(row) != m:
      raise ValueError("dp must be rectangular")

  if n != len(x) + 1 or m != len(y) + 1:
    raise ValueError("dp dimensions must be (len(x)+1) x (len(y)+1)")

  values, missing_mask, computed_values = _build_plot_matrix(dp, missing_value)

  cell_size = 0.6
  fig, ax = plt.subplots(
    figsize=(max(6, cell_size * m + 2), max(5, cell_size * n + 2))
  )

  image = None
  if len(computed_values) > 0:
    vmin = min(computed_values)
    vmax = max(computed_values)
    if math.isclose(vmin, vmax):
      vmax = vmin + 1.0

    image = ax.imshow(
      values,
      cmap=cmap,
      aspect="equal",
      origin="upper",
      vmin=vmin,
      vmax=vmax,
      alpha=computed_alpha
    )

    if show_colorbar:
      cbar = plt.colorbar(
        image,
        ax=ax,
        shrink=0.7,
        aspect=25,
        pad=0.02
      )
      cbar.set_label("DP score", fontsize=10)
      cbar.ax.tick_params(labelsize=8)

  missing_overlay = []
  for i in range(n):
    row = []
    for j in range(m):
      row.append(1 if missing_mask[i][j] else 0)
    missing_overlay.append(row)

  missing_cmap = ListedColormap([(1, 1, 1, 0), missing_color])
  ax.imshow(
    missing_overlay,
    cmap=missing_cmap,
    aspect="equal",
    origin="upper",
    vmin=0,
    vmax=1
  )

  x_labels = ["-"] + list(y)
  y_labels = ["-"] + list(x)

  ax.set_xticks(range(m))
  ax.set_yticks(range(n))
  ax.set_xticklabels(x_labels)
  ax.set_yticklabels(y_labels)

  ax.set_xlabel("Sequence Y")
  ax.set_ylabel("Sequence X")
  ax.set_title(title)

  ax.set_xticks([j - 0.5 for j in range(1, m)], minor=True)
  ax.set_yticks([i - 0.5 for i in range(1, n)], minor=True)
  ax.grid(which="minor", color="white", linestyle="-", linewidth=1)
  ax.tick_params(which="minor", bottom=False, left=False)

  if annotate:
    threshold = 0
    if len(computed_values) > 0:
      threshold = (min(computed_values) + max(computed_values)) / 2.0

    for i in range(n):
      for j in range(m):
        if missing_mask[i][j]:
          if missing_label != "":
            ax.text(
              j,
              i,
              missing_label,
              ha="center",
              va="center",
              color="black",
              fontsize=10
            )
        else:
          value = dp[i][j]
          text_color = "white" if float(value) > threshold else "black"
          ax.text(
            j,
            i,
            str(value),
            ha="center",
            va="center",
            color=text_color,
            fontsize=10
          )

  if path is not None and len(path) > 0:
    xs = [j for i, j in path]
    ys = [i for i, j in path]
    ax.plot(xs, ys, path_style, linewidth=2)

  plt.tight_layout()
  plt.show()


if __name__ == "__main__":
  x = "GATTACA"
  y = "GAATA"

  dp_banded = [
    [0, 1, 2, None, None, None],
    [1, 0, 1, 2, None, None],
    [2, 1, 0, 1, 2, None],
    [None, 2, 1, 1, 1, 2],
    [None, None, 2, 2, 1, 2],
    [None, None, None, 2, 2, 1],
    [None, None, None, None, 3, 2],
    [None, None, None, None, None, 3],
  ]

  plot_dp_matrix(
    dp_banded,
    x,
    y,
    title="Banded Edit Distance DP Matrix",
    annotate=True,
    missing_value=None,
    missing_label="",
    cmap="YlGnBu"
  )