#!/usr/bin/env python3
"""Independent exact reference for the dimension-six K6 Lorentz CSP.

Fix a required unit ``K6`` seed.  For every outside vertex ``x``, ``D_x`` is
the set of seed coordinates at which its defect vector is *allowed* to be
nonzero.  It is only an allowed support: candidate nonedges are never required
to be non-unit.

The K6 Schur identity has a Lorentz factor ``ell_x=(c_x,z_x)``.  Let ``L``
contain a required outside edge ``xy`` exactly when ``D_x`` and ``D_y`` are
disjoint.  Such an edge forces ``ell_x`` and ``ell_y`` to be Lorentz
orthogonal.  The finite necessary condition implemented here is:

* choose ``Z0={x: ell_x=0}``; its masks have size at least three, ``|Z0|<=6``,
  and admit a matching into the six seed coordinates;
* delete ``Z0`` from ``L`` and find its non-bipartite components;
* put every such component on one of the two lightlike projective lines;
* in either lightlike bin, the component vertices together with ``Z0`` have
  pairwise-orthogonal nonzero defect vectors, so their allowed masks must
  admit a matching into the six coordinates.

The search is existential over every possible ``Z0`` and both colors of every
non-bipartite component.  A seed is rejected only if this exhaustive search
fails.  Bipartite components are deliberately unconstrained.  The proof uses
``alpha(G)<=2``; every public graph-level entry point checks that precondition.

This is intentionally plain Python and structurally independent of the C
profiler.  Graphs are lists of integer adjacency masks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence


COORDINATES = 6


def vertices(mask: int) -> list[int]:
    """Return the increasing positions of the one bits in ``mask``."""

    answer: list[int] = []
    while mask:
        bit = mask & -mask
        answer.append(bit.bit_length() - 1)
        mask ^= bit
    return answer


def validate_graph(adj: Sequence[int], require_alpha_two: bool = True) -> None:
    """Validate a simple undirected graph and, by default, ``alpha<=2``."""

    n = len(adj)
    full = (1 << n) - 1
    for u, row in enumerate(adj):
        if row & ~full:
            raise ValueError(f"adjacency row {u} has a bit outside the graph")
        if row & (1 << u):
            raise ValueError(f"loop at vertex {u}")
        for v in range(u):
            if bool(row & (1 << v)) != bool(adj[v] & (1 << u)):
                raise ValueError(f"asymmetric adjacency at {u},{v}")
    if require_alpha_two and has_independent_triple(adj):
        raise ValueError("K6 Lorentz rule requires alpha(G)<=2")


def has_independent_triple(adj: Sequence[int]) -> bool:
    """Return whether the graph contains an independent set of size three."""

    n = len(adj)
    full = (1 << n) - 1
    for u in range(n):
        later = full & ~((1 << (u + 1)) - 1)
        nonneighbors = later & ~adj[u]
        remaining = nonneighbors
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            v = bit.bit_length() - 1
            if remaining & ~adj[v]:
                return True
    return False


def add_edge(adj: list[int], u: int, v: int) -> None:
    """Add a required undirected edge."""

    if u == v:
        raise ValueError("simple graphs have no loops")
    adj[u] |= 1 << v
    adj[v] |= 1 << u


def add_clique(adj: list[int], vs: Iterable[int]) -> None:
    """Add every required edge on ``vs``."""

    for u, v in combinations(tuple(vs), 2):
        add_edge(adj, u, v)


def is_clique(adj: Sequence[int], vs: Sequence[int]) -> bool:
    """Return whether all pairs in ``vs`` are required edges."""

    return all(adj[u] & (1 << v) for u, v in combinations(vs, 2))


def clique_masks(
    adj: Sequence[int], size: int, candidate_mask: int | None = None
) -> Iterator[int]:
    """Enumerate each clique of ``size`` once, as an absolute vertex mask."""

    if candidate_mask is None:
        candidate_mask = (1 << len(adj)) - 1

    def visit(candidates: int, need: int, chosen: int) -> Iterator[int]:
        if need == 0:
            yield chosen
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            yield from visit(candidates & adj[vertex], need - 1, chosen | bit)

    yield from visit(candidate_mask, size, 0)


def find_clique_mask(
    adj: Sequence[int], size: int, candidate_mask: int | None = None
) -> int:
    """Return one clique mask of the requested size, or zero."""

    return next(clique_masks(adj, size, candidate_mask), 0)


@dataclass(frozen=True)
class K6LorentzInstance:
    """The exact finite CSP induced by one K6 seed."""

    seed: tuple[int, ...]
    outside: tuple[int, ...]
    defects: tuple[int, ...]
    l_adj: tuple[int, ...]
    eligible_z0_mask: int

    def jsonable(self) -> dict:
        return {
            "seed": list(self.seed),
            "outside": list(self.outside),
            "allowed_defects": {
                str(x): [self.seed[q] for q in vertices(mask)]
                for x, mask in zip(self.outside, self.defects)
            },
            "eligible_Z0": [
                self.outside[i] for i in vertices(self.eligible_z0_mask)
            ],
            "L_edges": [
                [self.outside[i], self.outside[j]]
                for i in range(len(self.outside))
                for j in range(i)
                if self.l_adj[i] & (1 << j)
            ],
        }


@dataclass
class SeedDecision:
    """Exhaustive result for one fixed required K6 seed."""

    feasible: bool
    z0_subsets_considered: int
    z0_matchable: int
    colorings_considered: int
    chosen_z0: list[int] | None
    odd_components: list[list[int]] | None
    component_colors: list[int] | None
    bin_matchings: list[dict[str, int]] | None
    failure_counts: dict[str, int]
    witness: dict | None

    def jsonable(self) -> dict:
        return asdict(self)


@dataclass
class GraphDecision:
    """Graph-level result: one impossible K6 seed rejects the graph."""

    applicable: bool
    rejected: bool
    seeds_checked: int
    impossible_seeds: int
    total_z0_subsets_considered: int
    total_colorings_considered: int
    first_witness: dict | None
    feasible_z0_size_histogram: dict[str, int]
    note: str | None = None

    def jsonable(self) -> dict:
        return asdict(self)


def build_instance(adj: Sequence[int], seed: Sequence[int]) -> K6LorentzInstance:
    """Build relative defect masks and the disjoint-defect required-edge L."""

    if len(seed) != COORDINATES or not is_clique(adj, seed):
        raise ValueError("seed is not a required K6")
    seed = tuple(sorted(seed))
    seed_mask = sum(1 << q for q in seed)
    seed_position = {q: i for i, q in enumerate(seed)}
    outside = tuple(v for v in range(len(adj)) if not seed_mask & (1 << v))
    defects: list[int] = []
    for x in outside:
        relative = 0
        for q in seed:
            if not adj[x] & (1 << q):
                relative |= 1 << seed_position[q]
        defects.append(relative)

    l_adj = [0] * len(outside)
    for i, j in combinations(range(len(outside)), 2):
        if adj[outside[i]] & (1 << outside[j]) and not (
            defects[i] & defects[j]
        ):
            l_adj[i] |= 1 << j
            l_adj[j] |= 1 << i
    eligible = sum(
        1 << i for i, mask in enumerate(defects) if mask.bit_count() >= 3
    )
    return K6LorentzInstance(
        seed, outside, tuple(defects), tuple(l_adj), eligible
    )


def support_matching(
    selected_mask: int, defects: Sequence[int]
) -> dict[int, int] | None:
    """Find an injection from selected vertices into their allowed masks."""

    selected = vertices(selected_mask)
    if len(selected) > COORDINATES:
        return None
    selected.sort(key=lambda i: (defects[i].bit_count(), i))
    assignment: dict[int, int] = {}

    def visit(at: int, used: int) -> bool:
        if at == len(selected):
            return True
        local_vertex = selected[at]
        choices = defects[local_vertex] & ~used & ((1 << COORDINATES) - 1)
        while choices:
            coordinate = choices & -choices
            choices ^= coordinate
            assignment[local_vertex] = coordinate.bit_length() - 1
            if visit(at + 1, used | coordinate):
                return True
        assignment.pop(local_vertex, None)
        return False

    return assignment.copy() if visit(0, 0) else None


def nonbipartite_components(
    instance: K6LorentzInstance, deleted_mask: int
) -> list[int]:
    """Return the local masks of non-bipartite components of ``L-Z0``."""

    remaining = ((1 << len(instance.outside)) - 1) & ~deleted_mask
    answer: list[int] = []
    while remaining:
        root = remaining & -remaining
        root_index = root.bit_length() - 1
        component = 0
        frontier = root
        color_one = 0
        colored = root
        nonbipartite = False
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            component |= bit
            neighbors = instance.l_adj[vertex] & remaining
            want_one = not bool(color_one & bit)
            same_color = neighbors & (color_one if not want_one else ~color_one)
            if same_color & colored:
                nonbipartite = True
            new = neighbors & ~colored
            if want_one:
                color_one |= new
            colored |= new
            frontier |= new
        remaining &= ~component
        if nonbipartite:
            answer.append(component)
    return answer


def feasible_light_ray_assignments(
    instance: K6LorentzInstance,
    z0_mask: int,
    odd_components: Sequence[int] | None = None,
) -> list[tuple[tuple[int, ...], tuple[dict[int, int], dict[int, int]]]]:
    """Enumerate all two-color assignments satisfying both Hall matchings."""

    if odd_components is None:
        odd_components = nonbipartite_components(instance, z0_mask)
    answer = []
    for coloring in range(1 << len(odd_components)):
        bins = [z0_mask, z0_mask]
        for i, component in enumerate(odd_components):
            bins[(coloring >> i) & 1] |= component
        left = support_matching(bins[0], instance.defects)
        if left is None:
            continue
        right = support_matching(bins[1], instance.defects)
        if right is None:
            continue
        answer.append(
            (
                tuple((coloring >> i) & 1 for i in range(len(odd_components))),
                (left, right),
            )
        )
    return answer


def _z0_subsets(eligible: int) -> Iterator[int]:
    """Enumerate every eligible subset of size at most six, empty first."""

    bits = vertices(eligible)
    for size in range(min(COORDINATES, len(bits)) + 1):
        for chosen in combinations(bits, size):
            yield sum(1 << i for i in chosen)


def solve_seed(instance: K6LorentzInstance) -> SeedDecision:
    """Exhaustively solve the necessary Lorentz CSP for one K6 seed."""

    subsets_considered = 0
    matchable = 0
    colorings_considered = 0
    failures = {
        "Z0_support_matching": 0,
        "light_ray_bin_matching": 0,
    }
    failed_examples = []
    for z0_mask in _z0_subsets(instance.eligible_z0_mask):
        subsets_considered += 1
        zero_matching = support_matching(z0_mask, instance.defects)
        if zero_matching is None:
            failures["Z0_support_matching"] += 1
            if len(failed_examples) < 4:
                failed_examples.append(
                    {
                        "Z0": [instance.outside[i] for i in vertices(z0_mask)],
                        "reason": "Z0_support_matching",
                    }
                )
            continue
        matchable += 1
        odd = nonbipartite_components(instance, z0_mask)
        assignments = feasible_light_ray_assignments(instance, z0_mask, odd)
        colorings_considered += 1 << len(odd)
        if assignments:
            colors, bin_matchings = assignments[0]
            return SeedDecision(
                True,
                subsets_considered,
                matchable,
                colorings_considered,
                [instance.outside[i] for i in vertices(z0_mask)],
                [
                    [instance.outside[i] for i in vertices(component)]
                    for component in odd
                ],
                list(colors),
                [
                    {
                        str(instance.outside[i]): instance.seed[q]
                        for i, q in matching.items()
                    }
                    for matching in bin_matchings
                ],
                failures,
                None,
            )
        failures["light_ray_bin_matching"] += 1
        if len(failed_examples) < 4:
            failed_examples.append(
                {
                    "Z0": [instance.outside[i] for i in vertices(z0_mask)],
                    "reason": "light_ray_bin_matching",
                    "odd_components": [
                        [instance.outside[i] for i in vertices(component)]
                        for component in odd
                    ],
                    "colorings_checked": 1 << len(odd),
                }
            )

    witness = instance.jsonable()
    witness.update(
        {
            "failure_kind": "no_Z0_and_two_light_ray_assignment",
            "Z0_subsets_considered": subsets_considered,
            "Z0_support_matchable": matchable,
            "colorings_considered": colorings_considered,
            "failure_counts": failures,
            "failed_examples": failed_examples,
        }
    )
    return SeedDecision(
        False,
        subsets_considered,
        matchable,
        colorings_considered,
        None,
        None,
        None,
        None,
        failures,
        witness,
    )


def solve_k6_seed(adj: Sequence[int], seed: Sequence[int]) -> SeedDecision:
    """Validate the alpha-two precondition and solve one displayed K6."""

    validate_graph(adj, require_alpha_two=True)
    return solve_seed(build_instance(adj, seed))


def evaluate_graph(
    adj: Sequence[int], *, require_k6_only: bool = True, scan_all: bool = True
) -> GraphDecision:
    """Reject iff at least one required K6 seed has no feasible CSP choice."""

    validate_graph(adj, require_alpha_two=True)
    if require_k6_only and find_clique_mask(adj, 7):
        return GraphDecision(
            False,
            False,
            0,
            0,
            0,
            0,
            None,
            {},
            note="contains a required K7; target rule is restricted to K6-only",
        )

    checked = 0
    impossible = 0
    z0_total = 0
    coloring_total = 0
    first_witness = None
    z0_sizes: dict[str, int] = {}
    for seed_mask in clique_masks(adj, COORDINATES):
        checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_seed(instance)
        z0_total += decision.z0_subsets_considered
        coloring_total += decision.colorings_considered
        if decision.feasible:
            assert decision.chosen_z0 is not None
            key = str(len(decision.chosen_z0))
            z0_sizes[key] = z0_sizes.get(key, 0) + 1
        else:
            impossible += 1
            if first_witness is None:
                first_witness = decision.witness
            if not scan_all:
                break
    return GraphDecision(
        checked > 0,
        impossible > 0,
        checked,
        impossible,
        z0_total,
        coloring_total,
        first_witness,
        z0_sizes,
        note=None if checked else "graph has no required K6 seed",
    )


def _evaluate_record(record: dict) -> dict:
    decision = evaluate_graph(record["adjacency"], scan_all=True)
    return {
        "index": record["index"],
        "selection_priority": record["selection_priority"],
        "decision": decision.jsonable(),
    }


def evaluate_sample(sample: dict, workers: int = 1) -> dict:
    """Evaluate a committed sample, optionally across worker processes."""

    records = sample["graphs"]
    started = time.perf_counter()
    if workers == 1:
        graph_results = [_evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            graph_results = list(pool.map(_evaluate_record, records, chunksize=4))
    wall = time.perf_counter() - started

    rejected = [r for r in graph_results if r["decision"]["rejected"]]
    seeds = sum(r["decision"]["seeds_checked"] for r in graph_results)
    impossible_seeds = sum(
        r["decision"]["impossible_seeds"] for r in graph_results
    )
    z0_subsets = sum(
        r["decision"]["total_z0_subsets_considered"] for r in graph_results
    )
    colorings = sum(
        r["decision"]["total_colorings_considered"] for r in graph_results
    )
    impossible_histogram: dict[str, int] = {}
    chosen_z0_histogram: dict[str, int] = {}
    for result in graph_results:
        decision = result["decision"]
        key = str(decision["impossible_seeds"])
        impossible_histogram[key] = impossible_histogram.get(key, 0) + 1
        for size, count in decision["feasible_z0_size_histogram"].items():
            chosen_z0_histogram[size] = chosen_z0_histogram.get(size, 0) + count

    seed_witnesses = [
        {
            "index": result["index"],
            "witness": result["decision"]["first_witness"],
        }
        for result in rejected[:16]
    ]
    compact_graph_results = []
    for result in graph_results:
        decision = result["decision"].copy()
        witness = decision.pop("first_witness")
        decision["first_impossible_seed"] = (
            None if witness is None else witness["seed"]
        )
        compact_graph_results.append(
            {
                "index": result["index"],
                "selection_priority": result["selection_priority"],
                "decision": decision,
            }
        )

    return {
        "schema": 1,
        "filter": "K6_odd_component_two_light_ray_CSP",
        "mathematical_status": "exact necessary-condition rejection",
        "arithmetic": "integer bitsets and exhaustive finite search only",
        "candidate_nonedge_semantics": "allowed support only; never forced non-unit",
        "sample_schema": sample["schema"],
        "sample_graphs": len(records),
        "workers": workers,
        "wall_seconds": wall,
        "graphs_rejected": len(rejected),
        "rejection_rate": len(rejected) / len(records),
        "K6_seeds_checked": seeds,
        "impossible_K6_seeds": impossible_seeds,
        "Z0_subsets_considered": z0_subsets,
        "two_colorings_considered": colorings,
        "impossible_seed_count_per_graph_histogram": impossible_histogram,
        "chosen_Z0_size_per_feasible_seed_histogram": chosen_z0_histogram,
        "rejected_indices": [r["index"] for r in rejected],
        "seed_witnesses": seed_witnesses,
        "graph_results": compact_graph_results,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    with args.sample.open(encoding="utf-8") as stream:
        sample = json.load(stream)
    report = evaluate_sample(sample, args.workers)
    source = Path(__file__).resolve()
    report["inputs"] = {
        "sample": {"file": args.sample.name, "sha256": _sha256(args.sample)},
        "reference_source": {"file": source.name, "sha256": _sha256(source)},
        "corpus": sample["sources"]["corpus"],
        "old_kill_log": sample["sources"]["old_kill_log"],
        "previous_exact_profile": sample["sources"]["previous_exact_profile"],
    }
    report["runtime"] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "logical_cpu_count": os.cpu_count(),
        "workers": args.workers,
        "command": (
            f"python3 {source.name} {args.sample.name} --workers {args.workers}"
            + (f" --output {args.output.name}" if args.output else "")
        ),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(
            f"wrote {args.output}: {report['graphs_rejected']}/"
            f"{report['sample_graphs']} rejected; {report['wall_seconds']:.3f}s"
        )
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
