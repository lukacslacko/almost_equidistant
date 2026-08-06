#!/usr/bin/env python3
"""Exact K6 Lorentz-component rank refinement in dimension six.

Fix a required K6 seed and a possible actual zero-factor set ``Z0``.  The
disjoint-defect required-edge graph ``L-Z0`` has bipartite and non-bipartite
components.  The old filter constrains only the latter, because they are
forced onto the two Lorentz light rays.  This refinement also uses every
bipartite component ``C=A union B``.

If its nonzero Lorentz factors are generic, factors on each side are scalar
multiples of one non-lightlike vector and the two side directions are
Lorentz-orthogonal.  The Gram matrices of the defect vectors on ``A`` and
``B`` have off-diagonal graphs exactly ``G[A]`` and ``G[B]``.  Their spans,
and the span of the orthonormal vectors indexed by ``Z0``, are mutually
orthogonal in R^6.  Ordinary zero forcing therefore gives the necessary
condition

    |Z0| + |A|-Zf(G[A]) + |B|-Zf(G[B]) <= 6.

If the component direction is lightlike, every defect vector in the whole
component is instead mutually orthogonal, which is stronger and implies the
same displayed condition.  Isolated vertices contribute zero to the stated
zero-forcing lower bound, so they cause no exceptional case.

The search is existential over *every* eligible ``Z0`` of size at most six.
In particular, unlike the prior actual-support prototype, it does not omit
vertices lying in initially bipartite components: making one of those
vertices zero can change or split the very components tested here.  For each
remaining choice this module also applies the old non-bipartite light-ray
coloring and joint actual-support rules.  Candidate nonedges remain wholly
unconstrained; allowed defect coordinates may be zero.
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
from typing import Iterator, Sequence

from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    build_instance,
    clique_masks,
    find_clique_mask,
    support_matching,
    validate_graph,
    vertices,
)
from d6_k6_support_reference import actual_supports_for_bins


@dataclass(frozen=True)
class ZeroForcingResult:
    """One exact minimum zero-forcing set and its deterministic force list."""

    number: int
    initial: int
    forces: tuple[tuple[int, int], ...]


class ZeroForcingSolver:
    """Exact ordinary zero forcing on the small induced side graphs."""

    def __init__(self) -> None:
        self.cache: dict[tuple[int, ...], ZeroForcingResult] = {}
        self.cache_hits = 0
        self.initial_sets_checked = 0

    @staticmethod
    def closure(
        adj: Sequence[int], initial: int
    ) -> tuple[int, tuple[tuple[int, int], ...]]:
        black = initial
        forces: list[tuple[int, int]] = []
        while True:
            for vertex, row in enumerate(adj):
                if not (black & (1 << vertex)):
                    continue
                white = row & ~black
                if white and not (white & (white - 1)):
                    target = white.bit_length() - 1
                    black |= white
                    forces.append((vertex, target))
                    break
            else:
                return black, tuple(forces)

    def solve(self, adj: Sequence[int]) -> ZeroForcingResult:
        key = tuple(adj)
        cached = self.cache.get(key)
        if cached is not None:
            self.cache_hits += 1
            return cached
        n = len(key)
        full = (1 << n) - 1
        for size in range(n + 1):
            for chosen in combinations(range(n), size):
                initial = sum(1 << vertex for vertex in chosen)
                self.initial_sets_checked += 1
                closure, forces = self.closure(key, initial)
                if closure == full:
                    answer = ZeroForcingResult(size, initial, forces)
                    self.cache[key] = answer
                    return answer
        raise AssertionError("the full vertex set is always zero forcing")


@dataclass(frozen=True)
class LorentzComponent:
    """One connected component of L-Z0, with a valid coloring if bipartite."""

    component: int
    bipartite: bool
    side_a: int
    side_b: int


def lorentz_components(
    instance: K6LorentzInstance, deleted: int
) -> tuple[LorentzComponent, ...]:
    """Return every component of ``L-deleted`` and its bipartite sides."""

    full = (1 << len(instance.outside)) - 1
    if deleted & ~full:
        raise ValueError("deleted mask has a bit outside the K6 remainder")
    remaining = full & ~deleted
    answer: list[LorentzComponent] = []
    while remaining:
        root = remaining & -remaining
        component = 0
        colored = root
        side_b = 0
        frontier = root
        bipartite = True
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            component |= bit
            neighbors = instance.l_adj[vertex] & remaining
            vertex_in_b = bool(side_b & bit)
            same_side = side_b if vertex_in_b else (colored & ~side_b)
            if neighbors & same_side:
                bipartite = False
            new = neighbors & ~colored
            if not vertex_in_b:
                side_b |= new
            colored |= new
            frontier |= new
        remaining &= ~component
        side_b &= component
        answer.append(
            LorentzComponent(component, bipartite, component & ~side_b, side_b)
        )
    return tuple(answer)


def induced_required_graph(
    adj: Sequence[int], instance: K6LorentzInstance, selected: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Relabel ``G`` induced on one local outside-vertex mask."""

    local_vertices = tuple(vertices(selected))
    absolute = tuple(instance.outside[local] for local in local_vertices)
    rows = tuple(
        sum(
            1 << j
            for j, other in enumerate(absolute)
            if adj[vertex] & (1 << other)
        )
        for vertex in absolute
    )
    return rows, absolute


@dataclass(frozen=True)
class BipartiteRankCheck:
    """The exact zero-forcing lower bound for one bipartite L component."""

    component: tuple[int, ...]
    side_a: tuple[int, ...]
    side_b: tuple[int, ...]
    zero_forcing_a: int
    zero_forcing_b: int
    rank_lower_a: int
    rank_lower_b: int
    zero_factor_dimension: int
    required_dimension: int
    passed: bool

    def jsonable(self) -> dict:
        return asdict(self)


def check_bipartite_components(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    components: Sequence[LorentzComponent],
    solver: ZeroForcingSolver,
) -> tuple[bool, list[BipartiteRankCheck]]:
    """Apply the dimension-six inequality to every bipartite component."""

    checks: list[BipartiteRankCheck] = []
    zero_dimension = z0.bit_count()
    passed = True
    for component in components:
        if not component.bipartite:
            continue
        graph_a, absolute_a = induced_required_graph(
            adj, instance, component.side_a
        )
        graph_b, absolute_b = induced_required_graph(
            adj, instance, component.side_b
        )
        zf_a = solver.solve(graph_a).number
        zf_b = solver.solve(graph_b).number
        lower_a = len(absolute_a) - zf_a
        lower_b = len(absolute_b) - zf_b
        required = zero_dimension + lower_a + lower_b
        item = BipartiteRankCheck(
            tuple(instance.outside[i] for i in vertices(component.component)),
            absolute_a,
            absolute_b,
            zf_a,
            zf_b,
            lower_a,
            lower_b,
            zero_dimension,
            required,
            required <= COORDINATES,
        )
        checks.append(item)
        if not item.passed:
            passed = False
    return passed, checks


def _z0_subsets(eligible: int) -> Iterator[int]:
    """Enumerate every eligible Z0 of size at most six, with no component prune."""

    candidates = vertices(eligible)
    for size in range(min(COORDINATES, len(candidates)) + 1):
        for chosen in combinations(candidates, size):
            yield sum(1 << vertex for vertex in chosen)


def _pure_colorings(
    instance: K6LorentzInstance, z0: int, odd_components: Sequence[int]
) -> Iterator[tuple[tuple[int, ...], tuple[int, int]]]:
    """Yield symmetry-reduced old light-ray colorings passing mask Hall."""

    if not odd_components:
        if support_matching(z0, instance.defects) is not None:
            yield (), (z0, z0)
        return
    # Exchanging the two light rays is a global symmetry.
    for encoded in range(1 << (len(odd_components) - 1)):
        bins = [z0 | odd_components[0], z0]
        colors = [0]
        for position, component in enumerate(odd_components[1:]):
            color = (encoded >> position) & 1
            colors.append(color)
            bins[color] |= component
        if support_matching(bins[0], instance.defects) is None:
            continue
        if support_matching(bins[1], instance.defects) is None:
            continue
        yield tuple(colors), (bins[0], bins[1])


@dataclass
class RankSeedDecision:
    """Exhaustive combined decision for one fixed K6 seed."""

    feasible: bool
    z0_subsets_considered: int
    z0_support_matchable: int
    z0_bipartite_rank_passed: int
    bipartite_components_checked: int
    bipartite_rank_failures: int
    pure_colorings_considered: int
    pure_colorings_matchable: int
    support_searches: int
    support_dfs_nodes: int
    chosen_z0: list[int] | None
    chosen_bipartite_checks: list[dict] | None
    chosen_odd_components: list[list[int]] | None
    chosen_component_colors: list[int] | None
    chosen_supports: dict[str, list[int]] | None
    witness: dict | None

    def jsonable(self) -> dict:
        return asdict(self)


def solve_rank_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    zero_forcing: ZeroForcingSolver | None = None,
) -> RankSeedDecision:
    """Exhaust all Z0, bipartite rank checks, light rays, and actual supports."""

    if zero_forcing is None:
        zero_forcing = ZeroForcingSolver()
    counts = {
        "z0": 0,
        "z0_matchable": 0,
        "z0_rank_passed": 0,
        "components": 0,
        "rank_failures": 0,
        "colorings": 0,
        "colorings_matchable": 0,
        "searches": 0,
        "nodes": 0,
    }
    failed_examples: list[dict] = []

    # Deliberately use the full eligible mask.  The old optimization
    # ``eligible & initial_nonbipartite_union`` is invalid for this new rule.
    for z0 in _z0_subsets(instance.eligible_z0_mask):
        counts["z0"] += 1
        if support_matching(z0, instance.defects) is None:
            if len(failed_examples) < 16:
                failed_examples.append(
                    {
                        "Z0": [instance.outside[i] for i in vertices(z0)],
                        "reason": "Z0_support_matching",
                    }
                )
            continue
        counts["z0_matchable"] += 1

        components = lorentz_components(instance, z0)
        rank_ok, rank_checks = check_bipartite_components(
            adj, instance, z0, components, zero_forcing
        )
        counts["components"] += len(rank_checks)
        if not rank_ok:
            counts["rank_failures"] += 1
            if len(failed_examples) < 16:
                failed_examples.append(
                    {
                        "Z0": [instance.outside[i] for i in vertices(z0)],
                        "reason": "bipartite_component_rank",
                        "failed_components": [
                            check.jsonable() for check in rank_checks if not check.passed
                        ],
                    }
                )
            continue
        counts["z0_rank_passed"] += 1

        odd = [
            component.component
            for component in components
            if not component.bipartite
        ]
        raw_colorings = 1 if not odd else 1 << (len(odd) - 1)
        counts["colorings"] += raw_colorings
        for colors, bins in _pure_colorings(instance, z0, odd):
            counts["colorings_matchable"] += 1
            counts["searches"] += 1
            support = actual_supports_for_bins(instance, z0, bins)
            counts["nodes"] += support.dfs_nodes
            if not support.feasible:
                if len(failed_examples) < 16:
                    failed_examples.append(
                        {
                            "Z0": [instance.outside[i] for i in vertices(z0)],
                            "reason": "joint_actual_support",
                            "odd_components": [
                                [instance.outside[i] for i in vertices(component)]
                                for component in odd
                            ],
                            "component_colors": list(colors),
                            "support_DFS_nodes": support.dfs_nodes,
                        }
                    )
                continue
            assert support.supports is not None
            return RankSeedDecision(
                True,
                counts["z0"],
                counts["z0_matchable"],
                counts["z0_rank_passed"],
                counts["components"],
                counts["rank_failures"],
                counts["colorings"],
                counts["colorings_matchable"],
                counts["searches"],
                counts["nodes"],
                [instance.outside[i] for i in vertices(z0)],
                [check.jsonable() for check in rank_checks],
                [
                    [instance.outside[i] for i in vertices(component)]
                    for component in odd
                ],
                list(colors),
                {
                    str(instance.outside[i]): [
                        instance.seed[q] for q in vertices(actual)
                    ]
                    for i, actual in support.supports.items()
                },
                None,
            )

    witness = instance.jsonable()
    witness.update(
        {
            "failure_kind": "no_Z0_passing_bipartite_rank_and_light_support",
            "eligible_Z0_without_component_pruning": [
                instance.outside[i]
                for i in vertices(instance.eligible_z0_mask)
            ],
            "counts": counts.copy(),
            "failed_examples": failed_examples,
        }
    )
    return RankSeedDecision(
        False,
        counts["z0"],
        counts["z0_matchable"],
        counts["z0_rank_passed"],
        counts["components"],
        counts["rank_failures"],
        counts["colorings"],
        counts["colorings_matchable"],
        counts["searches"],
        counts["nodes"],
        None,
        None,
        None,
        None,
        None,
        witness,
    )


@dataclass
class RankGraphDecision:
    """A graph is rejected iff any required K6 seed exhausts every choice."""

    applicable: bool
    rejected: bool
    seeds_checked: int
    impossible_seeds: int
    z0_subsets_considered: int
    z0_support_matchable: int
    z0_bipartite_rank_passed: int
    bipartite_components_checked: int
    bipartite_rank_failures: int
    pure_colorings_considered: int
    pure_colorings_matchable: int
    support_searches: int
    support_dfs_nodes: int
    zero_forcing_initial_sets_checked: int
    zero_forcing_cache_entries: int
    zero_forcing_cache_hits: int
    feasible_z0_size_histogram: dict[str, int]
    first_witness: dict | None
    note: str | None = None

    def jsonable(self) -> dict:
        return asdict(self)


def evaluate_rank_graph(
    adj: Sequence[int], scan_all: bool = True
) -> RankGraphDecision:
    """Apply the combined refinement to every K6 seed of a K6-only graph."""

    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        return RankGraphDecision(
            False, False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, {}, None,
            note="contains a required K7; target population is K6-only",
        )

    names = (
        "z0_subsets_considered",
        "z0_support_matchable",
        "z0_bipartite_rank_passed",
        "bipartite_components_checked",
        "bipartite_rank_failures",
        "pure_colorings_considered",
        "pure_colorings_matchable",
        "support_searches",
        "support_dfs_nodes",
    )
    totals = {name: 0 for name in names}
    zero_forcing = ZeroForcingSolver()
    checked = 0
    impossible = 0
    sizes: dict[str, int] = {}
    first_witness = None
    for seed_mask in clique_masks(adj, COORDINATES):
        checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_rank_seed(adj, instance, zero_forcing)
        for name in names:
            totals[name] += getattr(decision, name)
        if decision.feasible:
            assert decision.chosen_z0 is not None
            key = str(len(decision.chosen_z0))
            sizes[key] = sizes.get(key, 0) + 1
        else:
            impossible += 1
            if first_witness is None:
                first_witness = decision.witness
            if not scan_all:
                break

    return RankGraphDecision(
        checked > 0,
        impossible > 0,
        checked,
        impossible,
        *(totals[name] for name in names),
        zero_forcing.initial_sets_checked,
        len(zero_forcing.cache),
        zero_forcing.cache_hits,
        sizes,
        first_witness,
        note=None if checked else "graph has no required K6 seed",
    )


def _evaluate_record(record: dict) -> dict:
    decision = evaluate_rank_graph(record["adjacency"], scan_all=True)
    return {"index": record["index"], "decision": decision.jsonable()}


def evaluate_input(sample: dict, workers: int) -> dict:
    """Evaluate a complete input artifact and return a compact exact report."""

    started = time.perf_counter()
    if workers == 1:
        raw = [_evaluate_record(record) for record in sample["graphs"]]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            raw = list(pool.map(_evaluate_record, sample["graphs"], chunksize=1))
    wall = time.perf_counter() - started
    rejected = [item for item in raw if item["decision"]["rejected"]]

    count_keys = (
        "seeds_checked",
        "impossible_seeds",
        "z0_subsets_considered",
        "z0_support_matchable",
        "z0_bipartite_rank_passed",
        "bipartite_components_checked",
        "bipartite_rank_failures",
        "pure_colorings_considered",
        "pure_colorings_matchable",
        "support_searches",
        "support_dfs_nodes",
        "zero_forcing_initial_sets_checked",
        "zero_forcing_cache_entries",
        "zero_forcing_cache_hits",
    )
    totals = {
        key: sum(item["decision"][key] for item in raw) for key in count_keys
    }
    witnesses = [
        {"index": item["index"], "witness": item["decision"]["first_witness"]}
        for item in rejected[:16]
    ]
    compact = []
    for item in raw:
        decision = item["decision"].copy()
        witness = decision.pop("first_witness")
        decision["first_impossible_seed"] = (
            None if witness is None else witness["seed"]
        )
        compact.append({"index": item["index"], "decision": decision})

    return {
        "schema": 1,
        "filter": "K6_bipartite_Lorentz_component_zero_forcing_rank",
        "status": "exact necessary-condition refinement",
        "input_graphs": len(raw),
        "workers": workers,
        "wall_seconds": wall,
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(raw) - len(rejected),
        "rejection_rate": len(rejected) / len(raw) if raw else 0.0,
        **totals,
        "rejected_indices": [item["index"] for item in rejected],
        "seed_witnesses": witnesses,
        "graph_results": compact,
        "soundness_scope": {
            "Z0_quantifier": (
                "every eligible subset of size at most 6; no initial-component prune"
            ),
            "bipartite_component_rule": (
                "|Z0|+|A|-Zf(G[A])+|B|-Zf(G[B]) <= 6"
            ),
            "generic_Lorentz_case": (
                "exact Gram patterns on mutually orthogonal side spans"
            ),
            "lightlike_Lorentz_case": (
                "whole component is orthonormal; stronger than tested rule"
            ),
            "nonbipartite_components": (
                "old two-light-ray bin matching plus joint actual-support CSP"
            ),
            "candidate_nonedges": (
                "unconstrained distances; alpha(G)<=2 only proves disjoint allowed masks"
            ),
            "floating_point_decisions": False,
        },
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    with args.input.open(encoding="utf-8") as stream:
        sample = json.load(stream)
    report = evaluate_input(sample, args.workers)
    source = Path(__file__).resolve()
    dependencies = (
        source.with_name("d6_k6_lorentz_reference.py"),
        source.with_name("d6_k6_support_reference.py"),
    )
    report["inputs"] = {
        "input": {"file": args.input.name, "sha256": sha256(args.input)},
        "source": {"file": source.name, "sha256": sha256(source)},
        "dependencies": [
            {"file": path.name, "sha256": sha256(path)} for path in dependencies
        ],
        **sample["sources"],
    }
    report["runtime"] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "logical_cpu_count": os.cpu_count(),
        "workers": args.workers,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(rendered, encoding="utf-8")
        temporary.replace(args.output)
        print(
            f"wrote {args.output}: {report['graphs_rejected']}/"
            f"{report['input_graphs']} rejected in {report['wall_seconds']:.3f}s"
        )
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
