#!/usr/bin/env python3
"""Exact nonlinear K7 special-H mask templates in rank six.

For a saturating six-clique C, write the six positive Sherman--Morrison
weights as coordinates of the basis.  This module recognizes two five-mask
templates whose full Schur equations have short positivity contradictions.
The proof is documented in ``d6_k7_special_h_templates.md``.

Candidate nonedges remain optional unit distances.  The 0/1 entries used here
are entries of the normalized K matrix only after fixing a K7 seed and a
zero-factor cover, exactly as in ``d6_k7_rank_reference.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as prior
import d6_k7_special_h_reference as href


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _outside_data(
    adj: Sequence[int], clique_mask: int
) -> tuple[list[int], dict[int, int], list[tuple[int, int]]]:
    clique = list(prior.bits(clique_mask))
    if len(clique) != 6:
        raise ValueError("special-H templates require a six-clique basis")
    if any(
        not (adj[first] & (1 << second))
        for first, second in combinations(clique, 2)
    ):
        raise ValueError("specified basis is not a required clique")
    remainder = [
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    ]
    masks = {
        vertex: sum(
            1 << coordinate
            for coordinate, basis in enumerate(clique)
            if adj[vertex] & (1 << basis)
        )
        for vertex in remainder
    }
    nonedges = [
        (first, second)
        for first, second in combinations(remainder, 2)
        if not (adj[first] & (1 << second))
    ]
    return remainder, masks, nonedges


def _type_one(
    remainder: Sequence[int], masks: dict[int, int], nonedges: Sequence[tuple[int, int]]
) -> dict | None:
    """Recognize the one-nonedge size-(3,5,3,3,3) template."""

    if len(remainder) != 5 or len(nonedges) != 1:
        return None
    first, second = nonedges[0]
    oriented = [
        (b_vertex, d_vertex)
        for b_vertex, d_vertex in ((first, second), (second, first))
        if masks[b_vertex].bit_count() == 5
        and masks[d_vertex].bit_count() == 3
    ]
    for b_vertex, d_vertex in oriented:
        triangle = [
            vertex for vertex in remainder
            if vertex not in (b_vertex, d_vertex)
        ]
        if any(masks[vertex].bit_count() != 3 for vertex in triangle):
            continue
        bmask, dmask = masks[b_vertex], masks[d_vertex]
        if any(masks[vertex] & ~bmask for vertex in triangle):
            continue
        common = masks[triangle[0]] & masks[triangle[1]] & masks[triangle[2]]
        union = masks[triangle[0]] | masks[triangle[1]] | masks[triangle[2]]
        if common.bit_count() != 1 or union.bit_count() != 4:
            continue
        if any(
            (masks[first_vertex] & masks[second_vertex]).bit_count() != 2
            for first_vertex, second_vertex in combinations(triangle, 2)
        ):
            continue
        if any((dmask & masks[vertex]) != common for vertex in triangle):
            continue
        if (dmask & bmask).bit_count() != 2:
            continue
        if not (dmask & bmask & common):
            continue
        if bmask != union | (dmask & bmask):
            continue
        return {
            "template": "type_I_one_nonedge",
            "B": b_vertex,
            "D": d_vertex,
            "triangle": triangle,
            "basis_masks": {str(vertex): masks[vertex] for vertex in remainder},
        }
    return None


def _type_two(
    remainder: Sequence[int], masks: dict[int, int], nonedges: Sequence[tuple[int, int]]
) -> dict | None:
    """Recognize the three-nonedge size-(4,5,4,2,4) template."""

    if len(remainder) != 5 or len(nonedges) != 3:
        return None
    nondegree = Counter(vertex for edge in nonedges for vertex in edge)
    centers = [vertex for vertex in remainder if nondegree[vertex] == 3]
    if len(centers) != 1:
        return None
    b_vertex = centers[0]
    if masks[b_vertex].bit_count() != 5:
        return None
    leaves = {vertex for edge in nonedges for vertex in edge if vertex != b_vertex}
    if len(leaves) != 3:
        return None
    d_candidates = [
        vertex for vertex in remainder
        if vertex != b_vertex and vertex not in leaves
    ]
    if len(d_candidates) != 1:
        return None
    d_vertex = d_candidates[0]
    triangle = sorted(leaves)
    if masks[d_vertex].bit_count() != 2:
        return None
    if any(masks[vertex].bit_count() != 4 for vertex in triangle):
        return None
    bmask, dmask = masks[b_vertex], masks[d_vertex]
    full = (1 << 6) - 1
    missing = full & ~bmask
    if missing.bit_count() != 1 or dmask & ~bmask:
        return None
    if any(not (masks[vertex] & missing) for vertex in triangle):
        return None
    common = masks[triangle[0]] & masks[triangle[1]] & masks[triangle[2]]
    union = masks[triangle[0]] | masks[triangle[1]] | masks[triangle[2]]
    if common.bit_count() != 2 or not (common & missing):
        return None
    central = common & bmask
    if central.bit_count() != 1 or union.bit_count() != 5:
        return None
    if any(
        (masks[first_vertex] & masks[second_vertex]).bit_count() != 3
        for first_vertex, second_vertex in combinations(triangle, 2)
    ):
        return None
    if any((dmask & masks[vertex]) != central for vertex in triangle):
        return None
    if bmask != (union & bmask) | dmask:
        return None
    return {
        "template": "type_II_three_nonedge_star",
        "B": b_vertex,
        "D": d_vertex,
        "triangle": triangle,
        "basis_masks": {str(vertex): masks[vertex] for vertex in remainder},
    }


def classify_template(adj: Sequence[int], clique_mask: int) -> dict | None:
    """Return a canonical exact obstruction witness, if either template fits."""

    if clique_mask.bit_count() != 6:
        return None
    remainder, masks, nonedges = _outside_data(adj, clique_mask)
    witness = _type_one(remainder, masks, nonedges)
    if witness is None:
        witness = _type_two(remainder, masks, nonedges)
    if witness is not None:
        witness = {
            **witness,
            "clique": list(prior.bits(clique_mask)),
            "remainder": remainder,
            "nonedges": [list(edge) for edge in nonedges],
        }
    return witness


@dataclass
class GraphTemplateResult:
    index: int | None
    seeds: int
    covers: int
    enhanced_passing_covers: int
    strict_h_passing_covers: int
    template_passing_covers: int
    saturating_cliques: int
    type_I_cliques: int
    type_II_cliques: int
    rejected: bool
    first_witness: dict | None = None
    elapsed_seconds: float = 0.0


def analyze_graph(graph: dict) -> GraphTemplateResult:
    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    seed_count = cover_count = enhanced_passes = strict_passes = template_passes = 0
    saturating_count = type_one_count = type_two_count = 0
    rejected = False
    first_witness = None
    for seed_mask in prior.clique_masks(adj, 7):
        seed_count += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(adj, seed_mask)
        total_term_rank = prior.matching_size(defects)
        seed_template_passes = 0
        covers = prior.eligible_covers(ladj, eligible)
        cover_count += len(covers)
        for zmask in covers:
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            enhanced_passes += 1
            nvertices = [
                outside[i] for i in range(len(outside))
                if not (zmask & (1 << i))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            cliques = list(prior.clique_masks(graph_n, analysis.k_rank_upper))
            saturating_count += len(cliques)
            if any(
                href.assess_saturating_clique(graph_n, clique).failed
                for clique in cliques
            ):
                continue
            strict_passes += 1
            witnesses = [
                witness
                for clique in cliques
                if (witness := classify_template(graph_n, clique)) is not None
            ]
            type_one_count += sum(
                witness["template"] == "type_I_one_nonedge"
                for witness in witnesses
            )
            type_two_count += sum(
                witness["template"] == "type_II_three_nonedge_star"
                for witness in witnesses
            )
            if witnesses:
                if first_witness is None:
                    first_witness = {
                        "seed": seed,
                        "zmask": zmask,
                        "nvertices": nvertices,
                        **witnesses[0],
                    }
                continue
            template_passes += 1
            seed_template_passes += 1
        if seed_template_passes == 0:
            rejected = True
    return GraphTemplateResult(
        index=graph.get("index"),
        seeds=seed_count,
        covers=cover_count,
        enhanced_passing_covers=enhanced_passes,
        strict_h_passing_covers=strict_passes,
        template_passing_covers=template_passes,
        saturating_cliques=saturating_count,
        type_I_cliques=type_one_count,
        type_II_cliques=type_two_count,
        rejected=rejected,
        first_witness=first_witness,
        elapsed_seconds=time.perf_counter() - started,
    )


def sample_report(sample: dict, indices: set[int], workers: int) -> dict:
    graphs = [
        graph for graph in sample["graphs"] if graph.get("index") in indices
    ]
    started = time.perf_counter()
    if workers == 1:
        results = [analyze_graph(graph) for graph in graphs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(analyze_graph, graphs, chunksize=1))
    rejected = sorted(
        result.index for result in results
        if result.index is not None and result.rejected
    )
    survivors = sorted(
        result.index for result in results
        if result.index is not None and not result.rejected
    )
    return {
        "schema": 1,
        "description": (
            "Exact nonlinear special-H template screen on the strict-H "
            "sample residue."
        ),
        "graphs": len(results),
        "workers": workers,
        "seeds": sum(result.seeds for result in results),
        "covers": sum(result.covers for result in results),
        "enhanced_passing_covers": sum(
            result.enhanced_passing_covers for result in results
        ),
        "strict_h_passing_covers": sum(
            result.strict_h_passing_covers for result in results
        ),
        "template_passing_covers": sum(
            result.template_passing_covers for result in results
        ),
        "saturating_cliques": sum(result.saturating_cliques for result in results),
        "template_clique_counts": {
            "type_I": sum(result.type_I_cliques for result in results),
            "type_II": sum(result.type_II_cliques for result in results),
        },
        "graph_counts": {
            "template_rejected": len(rejected),
            "survivors": len(survivors),
        },
        "graph_decisions": {
            "template_rejected": rejected,
            "survivors": survivors,
        },
        "runtime": {
            "wall_seconds": time.perf_counter() - started,
            "sum_graph_seconds": sum(result.elapsed_seconds for result in results),
            "maximum_graph_seconds": max(
                (result.elapsed_seconds for result in results), default=0.0
            ),
        },
        "per_graph": [asdict(result) for result in results],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("selection_report", type=Path)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sample = json.loads(args.sample.read_text(encoding="utf-8"))
    selection = json.loads(args.selection_report.read_text(encoding="utf-8"))
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("require 0 <= shard-index < shard-count")
    all_indices = sorted(selection["graph_decisions"]["combined_survivors"])
    indices = {
        index for position, index in enumerate(all_indices)
        if position % args.shard_count == args.shard_index
    }
    report = sample_report(sample, indices, args.workers)
    report["shard"] = {
        "count": args.shard_count,
        "index": args.shard_index,
        "selected_indices": sorted(indices),
    }
    report["input_sample"] = {
        "file": args.sample.name,
        "sha256": file_sha256(args.sample),
    }
    report["selection_report"] = {
        "file": args.selection_report.name,
        "sha256": file_sha256(args.selection_report),
    }
    report["sources"] = {
        Path(__file__).name: file_sha256(Path(__file__)),
        "d6_k7_special_h_reference.py": file_sha256(
            Path("d6_k7_special_h_reference.py")
        ),
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
