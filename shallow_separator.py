"""
ShallowSeparator: Implementation of the Plotkin-Rao-Smith (PRS) algorithm.

Given a graph G and integer parameter h, ShallowSeparator outputs either:
  - A balanced vertex separator S (every component of G \\ S has <= 2n/3 vertices), or
  - A Kh clique minor model of G.

This implements the win-win framework described in:
  Alexander, J. and Le, H. "Large Clique Minors Or Balanced Separators in
  (Road) Networks: An Experimental Study."

Usage:
  Edit the dataset list and output path in the __main__ block, then run:
      python shallow_separator.py
  Requires Python 3.10+ and the sortedcontainers package.

Input format:
  CSV files with a header row. Edge endpoints are expected in columns 1 and 2
  (0-indexed), with node IDs starting at 1. This matches the Li et al.
  [LCH+05] road network datasets. See initialize_G_and_H to adapt for other
  formats.
"""

import csv
import math
import os
import time
from collections import deque
from datetime import date

from sortedcontainers import SortedSet


class Node:
    """A vertex in the graph, tracking both the original edges (G) and the
    current working subgraph (H) as the algorithm progresses."""

    def __init__(self, node_id):
        self.ID = node_id
        self.neighbors_in_G = []   # fixed for the lifetime of the run
        self.neighbors_in_H = []   # shrinks as the algorithm removes vertices
        self.C_reference_node = None   # representative of this node's subgraph in K
        self.parent_in_Tv = None       # BFS-tree parent in the current iteration
        self.belonging_to_X = False    # True if this node was added to separator S

    def __eq__(self, other):
        return isinstance(other, Node) and self.ID == other.ID

    def __lt__(self, other):
        return isinstance(other, Node) and self.ID < other.ID

    def __hash__(self):
        return hash(self.ID)


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def find_n(edges_file_path):
    """Return the number of distinct vertices referenced in the edge-list CSV."""
    unique_nodes = SortedSet()
    with open(edges_file_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            unique_nodes.add(int(row[1]))
            unique_nodes.add(int(row[2]))
    return len(unique_nodes)


def create_node_objects(n):
    """Return a list of n Node objects with IDs 0 .. n-1."""
    return [Node(i) for i in range(n)]


def initialize_G_and_H(edges_file_path, node_object_list):
    """
    Populate neighbors_in_G and neighbors_in_H for every node from the CSV.
    Node IDs in the file are 1-indexed; they map to node_object_list[id - 1].
    """
    with open(edges_file_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            u = node_object_list[int(row[1]) - 1]
            v = node_object_list[int(row[2]) - 1]
            u.neighbors_in_G.append(v)
            v.neighbors_in_G.append(u)
            u.neighbors_in_H.append(v)
            v.neighbors_in_H.append(u)


def initialize_nodes_in_H_list(node_object_list):
    """Return a list of all nodes, representing the initial subgraph H = G."""
    return list(node_object_list)


# ---------------------------------------------------------------------------
# BFS and core subroutines
# ---------------------------------------------------------------------------

def BFSTree(v):
    """
    Build a BFS tree of the current subgraph H rooted at v.

    Returns:
        distances_list: list where distances_list[d] = [cumulative_node_count, [nodes at depth d]]
        Tv_depth: index of the deepest layer
    """
    v.parent_in_Tv = None
    distances_list = []
    visited = SortedSet([v])
    queue = deque([(v, 0)])
    total = 0

    while queue:
        current, depth = queue.popleft()
        if len(distances_list) <= depth:
            distances_list.append([0, []])
        distances_list[depth][1].append(current)
        total += 1
        distances_list[depth][0] = total
        for neighbor in current.neighbors_in_H:
            if neighbor not in visited:
                visited.add(neighbor)
                neighbor.parent_in_Tv = current
                queue.append((neighbor, depth + 1))

    return distances_list, len(distances_list) - 1


def minimalSubtree(Tv, K, v, nodes_in_H_list):
    """
    Compute Cv: the minimal subtree of Tv that connects v to at least one
    neighbor of every subgraph in K. Adds Cv to K and returns updated K.
    """
    Cv = SortedSet()
    if not K:
        Cv.add(nodes_in_H_list[0])
    else:
        done = SortedSet()
        for layer in Tv:
            for node in layer[1]:
                for neighbor in node.neighbors_in_G:
                    ref = neighbor.C_reference_node
                    if ref in K and ref not in done:
                        current = node
                        if current == v:
                            Cv.add(current)
                        while current.parent_in_Tv is not None:
                            Cv.add(current)
                            current = current.parent_in_Tv
                        done.add(ref)

    first_node = next(iter(Cv))
    for node in Cv:
        node.C_reference_node = first_node
    K.append(first_node)
    return Cv, K


def largestConnectedComponent(removed, nodes_in_H_list):
    """
    Return the largest connected component of H after removing vertices in `removed`.
    Trims neighbors_in_H so every node only lists neighbors within the component.
    """
    def bfs(start, visited):
        component = SortedSet()
        queue = deque([start])
        visited.add(start)
        while queue:
            current = queue.popleft()
            component.add(current)
            for nb in current.neighbors_in_H:
                if nb not in removed and nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        return component

    largest = SortedSet()
    visited = SortedSet()
    for node in nodes_in_H_list:
        if node not in removed and node not in visited:
            component = bfs(node, visited)
            if len(component) > len(largest):
                largest = component

    for node in largest:
        node.neighbors_in_H = [nb for nb in node.neighbors_in_H if nb in largest]
    return list(largest)


def trim(nodes_in_H_list, K):
    """
    Restore Invariant 3: remove from K every subgraph no longer adjacent to H,
    and rebuild neighbors_in_H for the surviving subgraph nodes.
    """
    visited_nodes = SortedSet()
    adjacent_refs = SortedSet()
    for node in nodes_in_H_list:
        for nb in node.neighbors_in_G:
            ref = nb.C_reference_node
            if ref in K:
                if nb not in visited_nodes:
                    visited_nodes.add(nb)
                    nb.neighbors_in_H = []
                adjacent_refs.add(ref)
                nb.neighbors_in_H.append(node)
    return nodes_in_H_list, list(adjacent_refs)


# ---------------------------------------------------------------------------
# Layer-selection strategies (Line 11 of Algorithm 1)
# ---------------------------------------------------------------------------

def _layer_valid(layer_nodes, nodes_above, nodes_below, l):
    """Return True if `layer_nodes` satisfies Equation (3) of the paper."""
    size = len(layer_nodes)
    return size <= nodes_above * (1 / l) and size <= nodes_below * (1 / l)


def findMedian(Tv, nodes_in_H_list):
    """
    Return (layer, index) for the median layer: the first layer whose cumulative
    node count reaches n/2, minimising the vertex-count imbalance above/below.
    """
    total = 0
    for i, layer in enumerate(Tv):
        total += len(layer[1])
        if total >= len(nodes_in_H_list) // 2:
            return layer[1], i


def find_X_A(Tv, l, nodes_in_H_list, Tv_depth):
    """Approach A: return the earliest layer from the root satisfying Eq. (3)."""
    for i in range(1, Tv_depth):
        layer = Tv[i][1]
        if _layer_valid(layer, Tv[i - 1][0], len(nodes_in_H_list) - Tv[i][0], l):
            return layer
    raise RuntimeError("No valid layer found (Approach A)")


def find_X_B(Tv, l, nodes_in_H_list, Tv_depth):
    """Approach B: return the layer with fewest vertices satisfying Eq. (3)."""
    best = None
    for i in range(1, Tv_depth):
        layer = Tv[i][1]
        if _layer_valid(layer, Tv[i - 1][0], len(nodes_in_H_list) - Tv[i][0], l):
            if best is None or len(layer) < len(best):
                best = layer
    if best is None:
        raise RuntimeError("No valid layer found (Approach B)")
    return best


def find_X_C(Tv, l, nodes_in_H_list, Tv_depth, median_layer, median_idx):
    """
    Approach C (default): search outward from the median layer, returning the
    first valid layer found (above or below). Ties are broken by choosing the
    smaller layer. See paper Section 2 for motivation.
    """
    nodes_above = Tv[median_idx - 1][0] if median_idx > 0 else 0
    nodes_below = len(nodes_in_H_list) - Tv[median_idx][0]
    if _layer_valid(median_layer, nodes_above, nodes_below, l):
        return median_layer

    i, j = median_idx, median_idx
    while i > 0 or j < Tv_depth:
        i -= 1
        j += 1
        above_valid = below_valid = False
        above_layer = below_layer = None

        if i > 0:
            above_layer = Tv[i][1]
            above_valid = _layer_valid(above_layer, Tv[i - 1][0],
                                       len(nodes_in_H_list) - Tv[i][0], l)
        if j < Tv_depth:
            below_layer = Tv[j][1]
            below_valid = _layer_valid(below_layer, Tv[j - 1][0],
                                       len(nodes_in_H_list) - Tv[j][0], l)

        if above_valid and below_valid:
            return above_layer if len(above_layer) <= len(below_layer) else below_layer
        if above_valid:
            return above_layer
        if below_valid:
            return below_layer

    raise RuntimeError("No valid layer found (Approach C)")


# ---------------------------------------------------------------------------
# ShallowSeparator (Algorithm 1)
# ---------------------------------------------------------------------------

def shallowSeparator(n, l, nodes_in_H_list, h, strategy='C'):
    """
    Run the PRS ShallowSeparator algorithm on the graph encoded in nodes_in_H_list.

    Parameters:
        n        : total number of vertices in the original graph G
        l        : depth/size tradeoff parameter (see Equation 1 in the paper)
        h        : target clique minor size; algorithm stops if Kh is found
        strategy : layer-selection strategy — 'A' (earliest), 'B' (smallest),
                   or 'C' (median-preferred, default)

    Returns a tuple (result, nodes_in_H, max_bfs_depth, iterations,
                     first_condition_calls, line11_calls, line11_median_returns)
    where `result` is either the string "MINOR MODEL" (a Kh minor was found)
    or a SortedSet S of separator vertices.
    """
    S = SortedSet()
    K = []  # reference nodes for subgraphs in the current clique minor
    iterations = max_depth = first_calls = line11_calls = line11_median = 0

    while len(nodes_in_H_list) >= (2 * n) / 3:
        iterations += 1
        v = nodes_in_H_list[0]
        Tv, Tv_depth = BFSTree(v)
        max_depth = max(max_depth, Tv_depth)

        if Tv_depth <= 2 * l * math.log(n):
            # Case 1: shallow tree — extend the clique minor
            first_calls += 1
            Cv, K = minimalSubtree(Tv, K, v, nodes_in_H_list)
            if len(K) == h:
                return ("MINOR MODEL", nodes_in_H_list, max_depth,
                        iterations, first_calls, line11_calls, line11_median)
            nodes_in_H_list = largestConnectedComponent(Cv, nodes_in_H_list)
        else:
            # Case 2: deep tree — find a separator layer X
            line11_calls += 1
            if strategy == 'A':
                X = SortedSet(find_X_A(Tv, l, nodes_in_H_list, Tv_depth))
            elif strategy == 'B':
                X = SortedSet(find_X_B(Tv, l, nodes_in_H_list, Tv_depth))
            else:
                median_layer, median_idx = findMedian(Tv, nodes_in_H_list)
                X = SortedSet(find_X_C(Tv, l, nodes_in_H_list, Tv_depth,
                                       median_layer, median_idx))
                if X == SortedSet(median_layer):
                    line11_median += 1

            for node in X:
                node.belonging_to_X = True
                S.add(node)
            nodes_in_H_list = largestConnectedComponent(X, nodes_in_H_list)

        nodes_in_H_list, K = trim(nodes_in_H_list, K)

    # Add all vertices from subgraphs in K to S (line 15 of Algorithm 1)
    for ref_node in K:
        visited = SortedSet([ref_node])
        queue = deque([ref_node])
        while queue:
            current = queue.popleft()
            S.add(current)
            for nb in current.neighbors_in_G:
                if nb.C_reference_node == ref_node and nb not in visited:
                    visited.add(nb)
                    queue.append(nb)

    return (S, nodes_in_H_list, max_depth, iterations,
            first_calls, line11_calls, line11_median)


# ---------------------------------------------------------------------------
# Exponential search + binary search for the largest h
# ---------------------------------------------------------------------------

def _run_at_h(dataset_file_path, n, const, h, strategy):
    """Build the graph fresh and run ShallowSeparator at a given h."""
    l = (const * math.sqrt(n)) / (h * math.sqrt(math.log(n)))
    node_object_list = create_node_objects(n)
    initialize_G_and_H(dataset_file_path, node_object_list)
    nodes_in_H_list = initialize_nodes_in_H_list(node_object_list)
    return shallowSeparator(n, l, nodes_in_H_list, h, strategy)


def find_largest_minor_h(dataset_file_path, n, const, strategy='C'):
    """
    Identify the largest h* for which ShallowSeparator returns a Kh minor model,
    using an exponential search followed by a binary search over [H/2, H].

    Returns:
        h_star        (int)   : largest h yielding a minor model (0 if none at h=1)
        final_result  (tuple) : return value of shallowSeparator run at h*+1,
                                which is guaranteed to produce a balanced separator
    """
    # Exponential search: double h until a separator is returned
    h = 1
    last_minor_h = 0
    while True:
        result = _run_at_h(dataset_file_path, n, const, h, strategy)
        if result[0] == "MINOR MODEL":
            last_minor_h = h
            h *= 2
        else:
            separator_h = h
            break

    # Binary search over [last_minor_h, separator_h] for the exact boundary
    lo, hi = last_minor_h, separator_h
    while hi - lo > 1:
        mid = (lo + hi) // 2
        result = _run_at_h(dataset_file_path, n, const, mid, strategy)
        if result[0] == "MINOR MODEL":
            lo = mid
        else:
            hi = mid

    h_star = lo  # largest h that returned a Kh minor model
    final_result = _run_at_h(dataset_file_path, n, const, h_star + 1, strategy)
    return h_star, final_result


# ---------------------------------------------------------------------------
# Experiment runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Dataset paths
    # Edge CSV format: header row; node IDs in columns 1 and 2 (1-indexed).
    # Compatible with the Li et al. [LCH+05] road network datasets.
    # ------------------------------------------------------------------
    datasets = [
        "datasets/Oldenburg Road Network's Edges.csv",
        "datasets/Oldenburg Road Network's Edges.csv",
        "datasets/San Francisco Road Network's Edges.csv",
        "datasets/San Joaquin Road Network's Edges.csv",
        "datasets/North America Road Network's Edges.csv",
        "datasets/California Road Network's Edges.csv",
    ]

    STRATEGY = 'C'                          # 'A', 'B', or 'C' (see paper Section 2)
    CONST_VALUES = [1, 2, 3, 5, 10, 20, 50, 100]
    OUTPUT_CSV = "outputs.csv"

    header = [
        "dataset", "n", "date", "constant_multiplier",
        "largest_h_minor_model", "smallest_h_separator",
        "separator_size", "max_bfs_depth", "iterations",
        "first_condition_calls", "line11_calls", "line11_median_returns",
        "time_elapsed_seconds",
    ]
    if not os.path.exists(OUTPUT_CSV):
        with open(OUTPUT_CSV, "w", newline="") as f:
            csv.writer(f).writerow(header)

    for dataset_path in datasets:
        n = find_n(dataset_path)
        for const in CONST_VALUES:
            start_time = time.time()
            h_star, result = find_largest_minor_h(dataset_path, n, const, STRATEGY)
            elapsed = time.time() - start_time

            S, _, max_depth, iters, first_calls, line11_calls, line11_median = result
            with open(OUTPUT_CSV, "a", newline="") as f:
                csv.writer(f).writerow([
                    dataset_path, n, date.today(), const,
                    h_star, h_star + 1,
                    len(S), max_depth, iters,
                    first_calls, line11_calls, line11_median,
                    round(elapsed, 3),
                ])
