#!/usr/bin/env python3
"""Exact pilot for PSD Z-matrix nullity inside the K6 fused-side layer.

This is deliberately a profiling prototype, not a theorem-level artifact.
It reuses the frozen K6 quantifiers but strengthens each generic side rank.
If

    Gram = I + kappa D (I + Adj(F)) D

is positive semidefinite, congruence by the nonsingular diagonal D gives
``diag(t)-H`` when ``kappa<0``.  This is a PSD Z-matrix, so each connected
block has nullity at most one.  When ``kappa>0`` the same conclusion follows
on every bipartite component of F after a signature switch; inertia remains
the bound on non-bipartite components.  All arithmetic below is exact.
"""

from __future__ import annotations

import argparse
import itertools
import json
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import d6_k6_fused_side_rank as frozen
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
    build_instance,
    clique_masks,
    support_matching,
    vertices,
)
from d6_k6_normal_block_support import check_block_subsets, z0_subsets
from d6_k6_normal_inertia import InertiaCache


class PsdZeroForcing:
    """Exact positive-semidefinite zero forcing by subset enumeration."""

    def __init__(self) -> None:
        self.cache: dict[tuple[int, ...], int] = {}

    @staticmethod
    def closure(rows: Sequence[int], initial: int) -> int:
        full = (1 << len(rows)) - 1
        black = initial
        while black != full:
            white = full & ~black
            white_components = []
            unseen = white
            while unseen:
                start = unseen & -unseen
                unseen ^= start
                component = start
                frontier = start
                while frontier:
                    bit = frontier & -frontier
                    frontier ^= bit
                    vertex = bit.bit_length() - 1
                    new = rows[vertex] & unseen
                    unseen &= ~new
                    component |= new
                    frontier |= new
                white_components.append(component)
            forced = 0
            for component in white_components:
                for vertex in vertices(black):
                    neighbors = rows[vertex] & component
                    if neighbors and not (neighbors & (neighbors - 1)):
                        forced |= neighbors
            forced &= ~black
            if not forced:
                break
            black |= forced
        return black

    def number(self, rows: Sequence[int]) -> int:
        key = tuple(rows)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        full = (1 << len(key)) - 1
        for size in range(len(key) + 1):
            for chosen in itertools.combinations(range(len(key)), size):
                initial = sum(1 << vertex for vertex in chosen)
                if self.closure(key, initial) == full:
                    self.cache[key] = size
                    return size
        raise AssertionError("the full set is PSD zero forcing")


@dataclass(frozen=True)
class TermRankResult:
    passed: bool
    first_failure: dict | None


def matching_size_masks(masks: Sequence[int]) -> int:
    reachable = {0}
    for mask in masks:
        following = set(reachable)
        for used in reachable:
            available = mask & ~used
            while available:
                bit = available & -available
                available ^= bit
                following.add(used | bit)
        reachable = following
    return max((used.bit_count() for used in reachable), default=0)


def term_rank_block_result(
    instance, z0: int, side_a: int, side_b: int, rank_a: int, rank_b: int
) -> TermRankResult:
    groups: list[tuple[str, int, int]] = [
        ("A", rank_a, side_a),
        ("B", rank_b, side_b),
    ]
    groups.extend(
        (f"Z0:{instance.outside[local]}", 1, 1 << local)
        for local in vertices(z0)
    )
    for selected in range(1, 1 << len(groups)):
        rank = 0
        vertex_mask = 0
        labels = []
        for position, (label, lower, group_vertices) in enumerate(groups):
            if selected & (1 << position):
                labels.append(label)
                rank += lower
                vertex_mask |= group_vertices
        masks = [instance.defects[local] for local in vertices(vertex_mask)]
        capacity = matching_size_masks(masks)
        if rank > capacity:
            return TermRankResult(False, {
                "blocks": labels,
                "rank_lower": rank,
                "term_rank_capacity": capacity,
                "vertices": [instance.outside[local] for local in vertices(vertex_mask)],
                "allowed_masks": masks,
            })
    return TermRankResult(True, None)


def graph_components(rows: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    unseen = (1 << len(rows)) - 1
    answer = []
    while unseen:
        start = unseen & -unseen
        unseen ^= start
        component = start
        frontier = start
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            new = rows[vertex] & unseen
            unseen &= ~new
            component |= new
            frontier |= new
        answer.append(tuple(vertices(component)))
    return tuple(answer)


def induced_rows(rows: Sequence[int], selected: Sequence[int]) -> tuple[int, ...]:
    position = {vertex: local for local, vertex in enumerate(selected)}
    return tuple(
        sum(1 << position[other] for other in selected if rows[vertex] & (1 << other))
        for vertex in selected
    )


def bipartite(rows: Sequence[int]) -> bool:
    colors: dict[int, int] = {}
    for start in range(len(rows)):
        if start in colors:
            continue
        colors[start] = 0
        stack = [start]
        while stack:
            vertex = stack.pop()
            for other in vertices(rows[vertex]):
                wanted = 1 - colors[vertex]
                if other in colors:
                    if colors[other] != wanted:
                        return False
                else:
                    colors[other] = wanted
                    stack.append(other)
    return True


def psd_side_rank_lower(
    graph_rows: Sequence[int], lorentz_sign: str, inertia: InertiaCache
) -> int:
    """Rank lower bound for one actual side Gram matrix."""

    total = 0
    for selected in graph_components(graph_rows):
        local = induced_rows(graph_rows, selected)
        size = len(local)
        pattern = tuple(row | (1 << index) for index, row in enumerate(local))
        positive, negative, _ = inertia.solve(pattern)
        if lorentz_sign == "negative":
            # diag(t)-H is an irreducible PSD Z-matrix on this block.
            lower = max(size - positive, size - 1)
        elif lorentz_sign == "positive":
            lower = size - negative
            if bipartite(local):
                # Signature-switch all positive edges to negative edges.
                lower = max(lower, size - 1)
        else:
            raise ValueError("unknown Lorentz sign")
        total += lower
    return total


@dataclass(frozen=True)
class ComponentDecision:
    passed: bool
    details: tuple[dict, ...]
    lightlike_passed: bool


def check_component(
    adj: Sequence[int], instance, z0: int, component, inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver, psd_zero_forcing: PsdZeroForcing,
) -> ComponentDecision:
    graph_a, absolute_a = induced_required_graph(adj, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adj, instance, component.side_b)
    zf_a = len(absolute_a) - zero_forcing.solve(graph_a).number
    zf_b = len(absolute_b) - zero_forcing.solve(graph_b).number
    psdzf_a = len(absolute_a) - psd_zero_forcing.number(graph_a)
    psdzf_b = len(absolute_b) - psd_zero_forcing.number(graph_b)

    details = []
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        psd_a = psd_side_rank_lower(graph_a, sign_a, inertia)
        psd_b = psd_side_rank_lower(graph_b, sign_b, inertia)
        rank_a = max(zf_a, psdzf_a, psd_a)
        rank_b = max(zf_b, psdzf_b, psd_b)
        support = term_rank_block_result(
            instance, z0, component.side_a, component.side_b, rank_a, rank_b
        )
        details.append({
            "case": name,
            "zero_forcing_rank_lower_A": zf_a,
            "zero_forcing_rank_lower_B": zf_b,
            "psd_zero_forcing_rank_lower_A": psdzf_a,
            "psd_zero_forcing_rank_lower_B": psdzf_b,
            "psd_rank_lower_A": psd_a,
            "psd_rank_lower_B": psd_b,
            "fused_rank_lower_A": rank_a,
            "fused_rank_lower_B": rank_b,
            "support_passed": support.passed,
            "first_failure": support.first_failure,
        })

    light_selected = z0 | component.component
    lightlike = (
        light_selected.bit_count() <= COORDINATES
        and support_matching(light_selected, instance.defects) is not None
    )
    return ComponentDecision(
        any(item["support_passed"] for item in details) or lightlike,
        tuple(details),
        lightlike,
    )


def solve_seed(adj: Sequence[int], instance) -> tuple[bool, dict]:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    psd_zero_forcing = PsdZeroForcing()
    counts = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_psd_passed": 0,
        "z0_psd_failed": 0,
        "z0_nonbipartite_failed": 0,
        "components": 0,
    }
    first_failure = None
    for z0 in z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            continue
        counts["z0_matchable"] += 1
        components = lorentz_components(instance, z0)
        decisions = []
        passed = True
        for component in components:
            if not component.bipartite:
                continue
            decision = check_component(
                adj, instance, z0, component, inertia, zero_forcing,
                psd_zero_forcing,
            )
            counts["components"] += 1
            decisions.append(decision)
            if not decision.passed:
                passed = False
        if not passed:
            counts["z0_psd_failed"] += 1
            if first_failure is None:
                first_failure = {
                    "Z0": frozen.absolute_vertices(instance, z0),
                    "components": [
                        {
                            "details": list(item.details),
                            "lightlike_passed": item.lightlike_passed,
                        }
                        for item in decisions if not item.passed
                    ],
                }
            continue
        counts["z0_psd_passed"] += 1
        nonbipartite = frozen.check_nonbipartite_system(
            instance, z0, components
        )
        if nonbipartite.passed:
            return True, counts
        counts["z0_nonbipartite_failed"] += 1
    counts["first_failure"] = first_failure
    return False, counts


def evaluate_record(record: dict) -> dict:
    adj = record["adjacency"]
    seeds = 0
    totals = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_psd_passed": 0,
        "z0_psd_failed": 0,
        "z0_nonbipartite_failed": 0,
        "components": 0,
    }
    witness = None
    rejected = False
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds += 1
        instance = build_instance(adj, vertices(seed_mask))
        feasible, counts = solve_seed(adj, instance)
        for name in totals:
            totals[name] += counts[name]
        if not feasible:
            rejected = True
            witness = {
                "seed": list(instance.seed),
                "first_failure": counts.get("first_failure"),
            }
            break
    return {
        "index": record["index"],
        "rejected": rejected,
        "seeds_checked": seeds,
        **totals,
        "witness": witness,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    _, records, indices = frozen.verify_inputs()
    records = [
        record for record in records
        if record["index"] not in set(
            json.loads(
                (frozen.ROOT / "d6_k6_fused_side_rank_report.json").read_text()
            )["rejected_indices"]
        )
    ]
    if args.limit is not None:
        records = records[: args.limit]
    started = time.perf_counter()
    if args.workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [item["index"] for item in results if item["rejected"]]
    report = {
        "status": "PILOT_COMPLETE",
        "input_graphs": len(results),
        "input_parent_indices_sha256": frozen.stable_hash(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(results) - len(rejected),
        "rejected_indices": rejected,
        "wall_seconds": time.perf_counter() - started,
        "workers": args.workers,
        "results": results,
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "survivor": "pilot non-rejection only",
            "arithmetic": "exact",
        },
    }
    if args.output:
        frozen.atomic_json(args.output.resolve(), report)
    print(json.dumps({key: report[key] for key in (
        "status", "input_graphs", "graphs_rejected", "graphs_surviving",
        "rejected_indices", "wall_seconds", "workers",
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
