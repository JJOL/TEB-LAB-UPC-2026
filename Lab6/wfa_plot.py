
"""
Plotting utilities for the edit-distance Wavefront Alignment Algorithm.

This module translates wavefront information into an equivalent partially
computed dynamic programming table and plots it as a heatmap.

Convention:
  - rows correspond to prefixes of the pattern
  - columns correspond to prefixes of the text
  - a cell (v, h) stores the minimum score s at which WFA reached it
  - cells not reached by WFA are stored as None
"""

import math
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def wavefront_cell(k, offset):
  v = offset - k
  h = offset
  return v, h


def build_wfa_dp_matrix(pattern, text, wavefronts):
  n = len(pattern)
  m = len(text)

  dp = []
  for i in range(n + 1):
    row = []
    for j in range(m + 1):
      row.append(None)
    dp.append(row)

  visited = set()

  for wavefront in wavefronts:
    score = wavefront["score"]
    starts = wavefront["starts"]
    offsets = wavefront["offsets"]

    for k in offsets:
      start = starts[k]
      end = offsets[k]

      h = start
      while h <= end:
        v = h - k
        if 0 <= v <= n and 0 <= h <= m:
          visited.add((v, h))
          if dp[v][h] is None or score < dp[v][h]:
            dp[v][h] = score
        h += 1

  return dp, visited


def build_plot_matrix(dp):
  values = []
  missing_mask = []
  computed_values = []

  i = 0
  while i < len(dp):
    value_row = []
    mask_row = []

    j = 0
    while j < len(dp[0]):
      cell = dp[i][j]
      missing = (cell is None)
      mask_row.append(missing)

      if missing:
        value_row.append(float("nan"))
      else:
        numeric = float(cell)
        value_row.append(numeric)
        computed_values.append(numeric)

      j += 1

    values.append(value_row)
    missing_mask.append(mask_row)
    i += 1

  return values, missing_mask, computed_values


def plot_wfa_dp_matrix(
  dp,
  pattern,
  text,
  path=None,
  title="WFA Equivalent DP Matrix",
  annotate=True,
  missing_label="",
  show_colorbar=True,
  cmap="YlOrRd",
  missing_color="#d9d9d9",
  path_style="-o",
  max_tick_labels=60
):
  n = len(dp)
  m = len(dp[0])
  total_cells = n * m

  values, missing_mask, computed_values = build_plot_matrix(dp)

  max_fig_w = 16
  max_fig_h = 12
  cell_size = 0.28
  fig, ax = plt.subplots(
    figsize=(
      min(max_fig_w, max(6, cell_size * m + 2)),
      min(max_fig_h, max(5, cell_size * n + 2))
    )
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
      vmax=vmax
    )

    if show_colorbar:
      cbar = plt.colorbar(
        image,
        ax=ax,
        shrink=0.7,
        aspect=25,
        pad=0.02
      )
      cbar.set_label("WFA score", fontsize=10)
      cbar.ax.tick_params(labelsize=8)

  missing_overlay = []
  i = 0
  while i < n:
    row = []
    j = 0
    while j < m:
      if missing_mask[i][j]:
        row.append(1)
      else:
        row.append(0)
      j += 1
    missing_overlay.append(row)
    i += 1

  missing_cmap = ListedColormap([(1, 1, 1, 0), missing_color])
  ax.imshow(
    missing_overlay,
    cmap=missing_cmap,
    aspect="equal",
    origin="upper",
    vmin=0,
    vmax=1
  )

  x_labels = ["-"] + list(text)
  y_labels = ["-"] + list(pattern)

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

  ax.set_xlabel("Text")
  ax.set_ylabel("Pattern")
  ax.set_title(title)

  show_grid = total_cells <= 8000
  if show_grid:
    ax.set_xticks([j - 0.5 for j in range(1, m)], minor=True)
    ax.set_yticks([i - 0.5 for i in range(1, n)], minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

  effective_annotate = annotate and total_cells <= 2500
  if effective_annotate:
    threshold = 0
    if len(computed_values) > 0:
      threshold = (min(computed_values) + max(computed_values)) / 2.0

    i = 0
    while i < n:
      j = 0
      while j < m:
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
        j += 1
      i += 1

  if path is not None and len(path) > 0:
    xs = []
    ys = []
    for i, j in path:
      xs.append(j)
      ys.append(i)
    if len(path) <= 400:
      ax.plot(xs, ys, path_style, linewidth=2)
    else:
      ax.plot(xs, ys, linewidth=1.6)

  plt.tight_layout()
  plt.show()


def plot_wfa_from_wavefronts(
  pattern,
  text,
  wavefronts,
  path=None,
  title="WFA Equivalent DP Matrix",
  annotate=True,
  missing_label=""
):
  dp, visited = build_wfa_dp_matrix(pattern, text, wavefronts)
  plot_wfa_dp_matrix(
    dp,
    pattern,
    text,
    path=path,
    title=title,
    annotate=annotate,
    missing_label=missing_label
  )
  return dp, visited
