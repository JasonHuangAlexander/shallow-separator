# ShallowSeparator

A Python implementation of the Plotkin-Rao-Smith (PRS) algorithm for finding balanced vertex separators and clique minors in real-world graphs.

This code accompanies the paper:

> Jason Alexander and Hung Le. *Large Clique Minors Or Balanced Separators in (Road) Networks: An Experimental Study.* University of Massachusetts Amherst.

---

## What it does

Given a graph *G* and an integer parameter *h*, **ShallowSeparator** outputs one of two things:

- A **balanced vertex separator** *S* — a set of vertices whose removal leaves every connected component with at most 2*n*/3 vertices, or
- A **K_h clique minor model** — *h* vertex-disjoint connected subgraphs of *G* such that every pair has an edge between them, certifying that *G* is not K_h-minor-free.

This is the *win-win framework*: the algorithm either solves the problem (separator) or explains why the standard theory predicts it cannot (large clique minor).

To find the **largest** h for which a K_h minor exists, the implementation uses an exponential search (doubling h until a separator is found) followed by a binary search over the resulting interval.

---

## Installation

```bash
pip install -r requirements.txt
```

Python 3.10 or later is required.

---

## Usage

1. Place your edge-list CSV files in a `datasets/` subdirectory.
2. Edit the `datasets` list and `STRATEGY` / `CONST_VALUES` constants at the bottom of `shallow_separator.py`.
3. Run:

```bash
python shallow_separator.py
```

Results are appended to `outputs.csv`.

---

## Input format

Each dataset must be a CSV file with:
- A **header row** (skipped automatically).
- One edge per subsequent row, with the two endpoint node IDs in **columns 1 and 2** (0-indexed), using **1-based node IDs**.

This matches the format of the Li et al. \[LCH+05\] road network datasets. If your files use different column positions or 0-based IDs, update `find_n` and `initialize_G_and_H` in `shallow_separator.py` accordingly.

---

## Output columns

| Column | Description |
|--------|-------------|
| `dataset` | Path to the input file |
| `n` | Number of vertices |
| `date` | Date the experiment was run |
| `constant_multiplier` | The constant *c* used to set ℓ = *c* · √*n* / (*h* √ln *n*) |
| `largest_h_minor_model` | Largest *h** for which a K_h minor was found |
| `smallest_h_separator` | *h** + 1; the value of *h* used to obtain the reported separator |
| `separator_size` | Number of vertices in the separator *S* |
| `max_bfs_depth` | Maximum BFS tree depth observed across all iterations |
| `iterations` | Total number of while-loop iterations |
| `first_condition_calls` | Iterations that took the shallow-tree (Case 1) branch |
| `line11_calls` | Iterations that took the deep-tree (Case 2) branch |
| `line11_median_returns` | (Strategy C only) Times the median layer itself was the chosen separator |
| `time_elapsed_seconds` | Wall-clock time for the entire run at this constant multiplier |

---

## Layer-selection strategies

Line 11 of Algorithm 1 selects a BFS layer to add to the separator. Three strategies are implemented (set `STRATEGY` in the `__main__` block):

| Strategy | Description |
|----------|-------------|
| `'A'` | **Earliest valid** — first layer from the root satisfying Equation (3) |
| `'B'` | **Smallest valid** — fewest-vertex layer satisfying Equation (3) |
| `'C'` | **Median-preferred** (default) — search outward from the median layer; return the first valid layer found, breaking ties by size |

Strategy C is theoretically motivated (median layers yield better balance) and empirically matched or outperformed A and B across all tested networks.

---

## Datasets

The experiments in the paper use:

- **Road networks (Li et al.)** — 5 classical datasets (Oldenburg, San Francisco, San Joaquin, North America, California). Available from the authors of: *F. Li, D. Cheng, M. Hadjieleftheriou, G. Kollios, S. Teng. On Trip Planning Queries in Spatial Databases. SSTD 2005.*
- **Road networks (Network Repository)** — 11 larger OSM-derived datasets. Available at <https://networkrepository.com>.
- **Social networks (SNAP)** — 10 datasets. Available at <https://snap.stanford.edu/data>.

Datasets are not included in this repository.

---

## Citation

If you use this code, please cite:

```
@article{alexander2025cliqueminor,
  title   = {Large Clique Minors Or Balanced Separators in (Road) Networks: An Experimental Study},
  author  = {Alexander, Jason and Le, Hung},
  year    = {2026},
  institution = {University of Massachusetts Amherst}
}
```
