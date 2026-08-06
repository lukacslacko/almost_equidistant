#!/usr/bin/env python3
"""Prototype the actual-support strengthening of the K6 Lorentz CSP.

This module is deliberately separate from the completed pure matching
milestone.  It imports only its transparent graph/seed construction helpers
and adds the next necessary condition.

For one candidate zero-factor set ``Z0`` and one coloring of the non-bipartite
components of ``L-Z0``, choose an *actual* support ``S_x subseteq D_x``.  A
zero-factor support has size at least three and is one shared choice in both
light-ray bins.  A nonzero light-ray support is merely nonempty.  Within each
orthonormal bin:

* no pair of supports may intersect in exactly one coordinate, because then
  the vectors' dot product consists of one nonzero summand and cannot vanish;
* the actual support family must admit a matching into the six coordinates,
  because an orthonormal family is linearly independent.

Intersection size at least two is accepted only as a necessary possibility:
the coefficients might cancel, but this screen does not claim they do.  There
is no constraint between opposite light-ray colors.  Candidate graph nonedges
remain unconstrained and allowed support coordinates may be unused.
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
    nonbipartite_components,
    support_matching,
    validate_graph,
    vertices,
)


@dataclass
class SupportSearchResult:
    """Result and work count for one fixed pair of light-ray bins."""

    feasible: bool
    dfs_nodes: int
    supports: dict[int, int] | None


@dataclass
class SupportSeedDecision:
    """Exhaustive actual-support decision for one fixed K6 seed."""

    feasible: bool
    z0_subsets_considered: int
    pure_colorings_considered: int
    pure_colorings_matchable: int
    support_searches: int
    support_dfs_nodes: int
    chosen_z0: list[int] | None
    chosen_odd_components: list[list[int]] | None
    chosen_component_colors: list[int] | None
    chosen_supports: dict[str, list[int]] | None
    witness: dict | None

    def jsonable(self) -> dict:
        return asdict(self)


@dataclass
class SupportGraphDecision:
    """One graph is rejected if any K6 seed exhausts all support choices."""

    applicable: bool
    rejected: bool
    seeds_checked: int
    impossible_seeds: int
    z0_subsets_considered: int
    pure_colorings_considered: int
    pure_colorings_matchable: int
    support_searches: int
    support_dfs_nodes: int
    feasible_z0_size_histogram: dict[str, int]
    first_witness: dict | None
    note: str | None = None

    def jsonable(self) -> dict:
        return asdict(self)


def support_domains(allowed: int, zero_factor: bool) -> tuple[int, ...]:
    """Enumerate every permitted nonempty actual support inside ``allowed``."""

    minimum = 3 if zero_factor else 1
    domains = []
    subset = allowed
    while subset:
        if subset.bit_count() >= minimum:
            domains.append(subset)
        subset = (subset - 1) & allowed
    # Small supports tend to expose contradictions quickly; lexical tie-break
    # makes the witness deterministic.
    domains.sort(key=lambda mask: (mask.bit_count(), mask))
    return tuple(domains)


def intersection_compatible(left: int, right: int) -> bool:
    """Two orthogonal nonzero vectors cannot share exactly one support slot."""

    return (left & right).bit_count() != 1


def _actual_matching(selected: int, supports: dict[int, int], n: int) -> bool:
    """Check structural full row rank for the assigned part of one bin."""

    masks = [0] * n
    assigned = 0
    for vertex, support in supports.items():
        bit = 1 << vertex
        if selected & bit:
            assigned |= bit
            masks[vertex] = support
    return support_matching(assigned, masks) is not None


def actual_supports_for_bins(
    instance: K6LorentzInstance, z0: int, bins: tuple[int, int]
) -> SupportSearchResult:
    """Search actual supports jointly, sharing every Z0 support across bins."""

    if (bins[0] & z0) != z0 or (bins[1] & z0) != z0:
        raise ValueError("both bins must contain Z0")
    if (bins[0] & bins[1]) != z0:
        raise ValueError("the two bins may overlap only in Z0")
    involved = bins[0] | bins[1]
    n = len(instance.outside)
    domains: dict[int, tuple[int, ...]] = {}
    peer_masks: dict[int, int] = {}
    for vertex in vertices(involved):
        domains[vertex] = support_domains(
            instance.defects[vertex], bool(z0 & (1 << vertex))
        )
        peers = 0
        for bin_mask in bins:
            if bin_mask & (1 << vertex):
                peers |= bin_mask & ~(1 << vertex)
        peer_masks[vertex] = peers

    assignments: dict[int, int] = {}
    nodes = 0

    def viable_domain(vertex: int) -> list[int]:
        answer = []
        assigned_peers = peer_masks[vertex]
        for support in domains[vertex]:
            if all(
                intersection_compatible(support, other_support)
                for other, other_support in assignments.items()
                if assigned_peers & (1 << other)
            ):
                answer.append(support)
        return answer

    def visit(unassigned: int) -> bool:
        nonlocal nodes
        nodes += 1
        if not unassigned:
            return True

        # Exact MRV ordering, with constrained degree as a deterministic tie.
        best_vertex = -1
        best_choices: list[int] | None = None
        best_key = None
        for vertex in vertices(unassigned):
            choices = viable_domain(vertex)
            key = (len(choices), -(peer_masks[vertex] & unassigned).bit_count(), vertex)
            if best_key is None or key < best_key:
                best_key = key
                best_vertex = vertex
                best_choices = choices
            if not choices:
                return False
        assert best_choices is not None

        bit = 1 << best_vertex
        for support in best_choices:
            assignments[best_vertex] = support
            if all(_actual_matching(bin_mask, assignments, n) for bin_mask in bins):
                if visit(unassigned ^ bit):
                    return True
            assignments.pop(best_vertex, None)
        return False

    feasible = visit(involved)
    return SupportSearchResult(
        feasible, nodes, assignments.copy() if feasible else None
    )


def _z0_subsets(eligible: int) -> Iterator[int]:
    bits = vertices(eligible)
    for size in range(min(COORDINATES, len(bits)) + 1):
        for chosen in combinations(bits, size):
            yield sum(1 << vertex for vertex in chosen)


def _pure_colorings(
    instance: K6LorentzInstance, z0: int, components: Sequence[int]
) -> Iterator[tuple[tuple[int, ...], tuple[int, int]]]:
    """Yield symmetry-reduced colorings passing the old allowed-mask matchings."""

    if not components:
        bins = (z0, z0)
        if support_matching(z0, instance.defects) is not None:
            yield (), bins
        return
    # Exchanging the names of the two light rays is a symmetry.  Fix the first
    # component in color zero and enumerate the remaining components.
    for encoded in range(1 << (len(components) - 1)):
        bins = [z0 | components[0], z0]
        colors = [0]
        for i, component in enumerate(components[1:], 1):
            color = (encoded >> (i - 1)) & 1
            colors.append(color)
            bins[color] |= component
        if support_matching(bins[0], instance.defects) is None:
            continue
        if support_matching(bins[1], instance.defects) is None:
            continue
        yield tuple(colors), (bins[0], bins[1])


def solve_support_seed(instance: K6LorentzInstance) -> SupportSeedDecision:
    """Exhaust every Z0, component coloring, and joint actual-support choice."""

    initial_odd = nonbipartite_components(instance, 0)
    initial_odd_union = 0
    for component in initial_odd:
        initial_odd_union |= component
    # Choosing a zero in an initially bipartite connected component can never
    # help: deletion preserves bipartiteness, while omitting that zero only
    # removes one vector from both support bins.  Tree attachments to an odd
    # component remain eligible because the whole initial component is kept.
    eligible = instance.eligible_z0_mask & initial_odd_union

    z0_count = 0
    colorings_count = 0
    matchable_count = 0
    searches = 0
    nodes = 0
    failed_support_examples = []
    for z0 in _z0_subsets(eligible):
        z0_count += 1
        if support_matching(z0, instance.defects) is None:
            continue
        components = nonbipartite_components(instance, z0)
        raw_colorings = 1 if not components else 1 << (len(components) - 1)
        colorings_count += raw_colorings
        for colors, bins in _pure_colorings(instance, z0, components):
            matchable_count += 1
            searches += 1
            result = actual_supports_for_bins(instance, z0, bins)
            nodes += result.dfs_nodes
            if not result.feasible:
                if len(failed_support_examples) < 16:
                    failed_support_examples.append(
                        {
                            "Z0": [
                                instance.outside[i] for i in vertices(z0)
                            ],
                            "odd_components": [
                                [
                                    instance.outside[i]
                                    for i in vertices(component)
                                ]
                                for component in components
                            ],
                            "component_colors": list(colors),
                            "bins": [
                                [
                                    instance.outside[i]
                                    for i in vertices(bin_mask)
                                ]
                                for bin_mask in bins
                            ],
                            "support_DFS_nodes": result.dfs_nodes,
                        }
                    )
                continue
            assert result.supports is not None
            return SupportSeedDecision(
                True,
                z0_count,
                colorings_count,
                matchable_count,
                searches,
                nodes,
                [instance.outside[i] for i in vertices(z0)],
                [
                    [instance.outside[i] for i in vertices(component)]
                    for component in components
                ],
                list(colors),
                {
                    str(instance.outside[i]): [
                        instance.seed[q] for q in vertices(support)
                    ]
                    for i, support in result.supports.items()
                },
                None,
            )

    witness = instance.jsonable()
    witness.update(
        {
            "failure_kind": "no_joint_actual_support_assignment",
            "eligible_Z0_after_bipartite_component_pruning": [
                instance.outside[i] for i in vertices(eligible)
            ],
            "Z0_subsets_considered": z0_count,
            "pure_colorings_considered": colorings_count,
            "pure_colorings_allowed_mask_matchable": matchable_count,
            "actual_support_searches": searches,
            "actual_support_DFS_nodes": nodes,
            "failed_support_examples": failed_support_examples,
        }
    )
    return SupportSeedDecision(
        False,
        z0_count,
        colorings_count,
        matchable_count,
        searches,
        nodes,
        None,
        None,
        None,
        None,
        witness,
    )


def evaluate_support_graph(adj: Sequence[int], scan_all: bool = True) -> SupportGraphDecision:
    """Apply the support strengthening to every K6 seed of a K6-only graph."""

    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        return SupportGraphDecision(
            False, False, 0, 0, 0, 0, 0, 0, 0, {}, None,
            note="contains a required K7; target population is K6-only",
        )

    totals = {
        "z0": 0,
        "colorings": 0,
        "matchable": 0,
        "searches": 0,
        "nodes": 0,
    }
    checked = 0
    impossible = 0
    sizes: dict[str, int] = {}
    first_witness = None
    for seed_mask in clique_masks(adj, 6):
        checked += 1
        decision = solve_support_seed(build_instance(adj, vertices(seed_mask)))
        totals["z0"] += decision.z0_subsets_considered
        totals["colorings"] += decision.pure_colorings_considered
        totals["matchable"] += decision.pure_colorings_matchable
        totals["searches"] += decision.support_searches
        totals["nodes"] += decision.support_dfs_nodes
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
    return SupportGraphDecision(
        checked > 0,
        impossible > 0,
        checked,
        impossible,
        totals["z0"],
        totals["colorings"],
        totals["matchable"],
        totals["searches"],
        totals["nodes"],
        sizes,
        first_witness,
        note=None if checked else "graph has no required K6 seed",
    )


def _evaluate_record(record: dict) -> dict:
    decision = evaluate_support_graph(record["adjacency"], scan_all=True)
    return {"index": record["index"], "decision": decision.jsonable()}


def evaluate_sample(sample: dict, workers: int) -> dict:
    started = time.perf_counter()
    if workers == 1:
        raw = [_evaluate_record(record) for record in sample["graphs"]]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            raw = list(pool.map(_evaluate_record, sample["graphs"], chunksize=2))
    wall = time.perf_counter() - started
    rejected = [item for item in raw if item["decision"]["rejected"]]

    count_keys = (
        "seeds_checked",
        "impossible_seeds",
        "z0_subsets_considered",
        "pure_colorings_considered",
        "pure_colorings_matchable",
        "support_searches",
        "support_dfs_nodes",
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
        decision["first_impossible_seed"] = None if witness is None else witness["seed"]
        compact.append({"index": item["index"], "decision": decision})

    return {
        "schema": 1,
        "filter": "K6_joint_actual_support_intersection_CSP",
        "status": "exact necessary-condition prototype",
        "sample_graphs": len(raw),
        "workers": workers,
        "wall_seconds": wall,
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(raw) - len(rejected),
        "rejection_rate": len(rejected) / len(raw),
        **totals,
        "rejected_indices": [item["index"] for item in rejected],
        "seed_witnesses": witnesses,
        "graph_results": compact,
        "soundness_scope": {
            "Z0_support_size": "at least 3",
            "nonzero_light_ray_support_size": "at least 1",
            "shared_Z0_supports": True,
            "within_bin_pair_intersection": "0 or at least 2",
            "actual_support_matching_per_bin": True,
            "opposite_color_constraints": "none",
            "candidate_nonedges": "unconstrained; D is only an allowed support",
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
    dependency = source.with_name("d6_k6_lorentz_reference.py")
    report["inputs"] = {
        "sample": {"file": args.sample.name, "sha256": sha256(args.sample)},
        "source": {"file": source.name, "sha256": sha256(source)},
        "pure_lorentz_reference_dependency": {
            "file": dependency.name,
            "sha256": sha256(dependency),
        },
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
        args.output.write_text(rendered, encoding="utf-8")
        print(
            f"wrote {args.output}: {report['graphs_rejected']}/"
            f"{report['sample_graphs']} rejected in {report['wall_seconds']:.3f}s"
        )
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
