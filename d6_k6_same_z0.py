#!/usr/bin/env python3
"""Exact same-Z0 conjunction of the strongest current K6 necessary rules.

For each required K6 seed, a realization has one actual zero-Lorentz-factor
set Z0.  Earlier filters existentially selected Z0 independently.  This
layer preserves the quantifier and requires one Z0 to satisfy both:

* every bipartite Lorentz component passes the exact normal block-support
  inequalities in one generic orientation or its lightlike alternative;
* all non-bipartite components admit one two-light-ray coloring and one joint
  actual-support assignment, with the same supports for Z0 in both bins.

Candidate nonedges remain unconstrained and allowed coordinates may be zero.
All decisions are exact; floating point is used only for elapsed time.
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
from pathlib import Path
from typing import Sequence

from d6_k6_bipartite_rank_reference import _pure_colorings, lorentz_components
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
from d6_k6_normal_block_support import check_component, z0_subsets
from d6_k6_normal_inertia import InertiaCache
from d6_k6_support_reference import actual_supports_for_bins
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
SCHEMA = "d6-k6-same-z0-v1"
EXPECTED_INPUT = 990
EXPECTED_INDICES_SHA256 = (
    "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
)
EXPECTED = {
    "d6_current_residue_manifest.json": (
        "d06cdb23b03b1f8523fbe3c3ede5bc297ace21a844282830c5cf51a6cee3512d"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "d6_k6_support_reference.py": (
        "d6481137ef49d88154882285660dd755eb7cb652ab87280f05c3792de9742577"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_normal_inertia.py": (
        "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023"
    ),
    "d6_k6_normal_block_support.py": (
        "14faafc6b2f9c3c34c954bb9cbaa21e9a7f646cfb76fcdfb30dfe7af742cf9cf"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


COUNTER_FIELDS = (
    "z0_considered",
    "z0_matchable",
    "z0_block_passed",
    "z0_block_failed",
    "z0_nonbipartite_failed",
    "bipartite_components_checked",
    "bipartite_components_failed",
    "generic_orientation_cases_checked",
    "generic_orientation_support_new_failures",
    "block_subsets_checked",
    "pure_colorings_considered",
    "pure_colorings_matchable",
    "joint_support_searches",
    "joint_support_dfs_nodes",
)


def empty_counts() -> dict[str, int]:
    return {name: 0 for name in COUNTER_FIELDS}


def add_counts(target: dict[str, int], source: dict[str, int]) -> None:
    for name in COUNTER_FIELDS:
        target[name] += source[name]


def absolute_vertices(instance: K6LorentzInstance, mask: int) -> list[int]:
    return [instance.outside[local] for local in vertices(mask)]


def compact_block_failure(
    instance: K6LorentzInstance, component, decision
) -> dict:
    generic = []
    for case in decision.generic_cases:
        generic.append({
            "case": case["case"],
            "rank_lower_A": case["rank_lower_A"],
            "rank_lower_B": case["rank_lower_B"],
            "ambient_dimension_lower": case["ambient_dimension_lower"],
            "ambient_dimension_passed": case["ambient_dimension_passed"],
            "passed": case["passed"],
            "first_failure": case["first_failure"],
        })
    return {
        "component": absolute_vertices(instance, component.component),
        "side_A": absolute_vertices(instance, component.side_a),
        "side_B": absolute_vertices(instance, component.side_b),
        "generic_cases": generic,
        "lightlike_case": decision.lightlike_case,
    }


def rejected_seed_cross_classification(
    adj: Sequence[int], instance: K6LorentzInstance
) -> dict:
    """Recompute both systems independently for every Z0 of a rejected seed.

    The main search stops testing the nonbipartite system as soon as the block
    system fails.  This certificate-only audit deliberately does not: it
    records the full two-by-two truth table, making the quantifier obstruction
    directly checkable.  Its work is not included in production counters.
    """

    inertia = InertiaCache()
    rows = []
    histogram = {
        "both": 0,
        "block_only": 0,
        "nonbipartite_only": 0,
        "neither": 0,
    }
    representatives: dict[str, list[int]] = {}
    for z0 in z0_subsets(instance.eligible_z0_mask):
        z0_absolute = absolute_vertices(instance, z0)
        if support_matching(z0, instance.defects) is None:
            block_passed = False
            nonbipartite_passed = False
        else:
            local_components = lorentz_components(instance, z0)
            block_passed = all(
                check_component(adj, instance, z0, component, inertia).passed
                for component in local_components
                if component.bipartite
            )
            odd = [
                component.component
                for component in local_components
                if not component.bipartite
            ]
            nonbipartite_passed = False
            for _, bins in _pure_colorings(instance, z0, odd):
                if actual_supports_for_bins(instance, z0, bins).feasible:
                    nonbipartite_passed = True
                    break
        if block_passed and nonbipartite_passed:
            category = "both"
        elif block_passed:
            category = "block_only"
        elif nonbipartite_passed:
            category = "nonbipartite_only"
        else:
            category = "neither"
        histogram[category] += 1
        representatives.setdefault(category, z0_absolute)
        rows.append({
            "Z0": z0_absolute,
            "block_support_passed": block_passed,
            "nonbipartite_joint_support_passed": nonbipartite_passed,
            "category": category,
        })
    return {
        "histogram": histogram,
        "representative_Z0": representatives,
        "every_Z0_passes_exactly_one_system": (
            histogram["both"] == 0
            and histogram["neither"] == 0
            and histogram["block_only"] > 0
            and histogram["nonbipartite_only"] > 0
        ),
        "choices": rows,
    }


@dataclass
class SameZ0SeedDecision:
    feasible: bool
    counts: dict[str, int]
    chosen_z0: list[int] | None
    chosen_component_colors: list[int] | None
    chosen_actual_supports: dict[str, list[int]] | None
    certificate: dict | None


def solve_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    inertia: InertiaCache | None = None,
) -> SameZ0SeedDecision:
    """Exhaust the one common Z0 for one fixed required K6 seed."""

    if inertia is None:
        inertia = InertiaCache()
    counts = empty_counts()
    z0_failures = []

    for z0 in z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        z0_absolute = absolute_vertices(instance, z0)
        if support_matching(z0, instance.defects) is None:
            z0_failures.append({
                "Z0": z0_absolute,
                "reason": "Z0_allowed_mask_matching",
            })
            continue
        counts["z0_matchable"] += 1
        components = lorentz_components(instance, z0)

        failed_blocks = []
        for component in components:
            if not component.bipartite:
                continue
            decision = check_component(adj, instance, z0, component, inertia)
            counts["bipartite_components_checked"] += 1
            counts["generic_orientation_cases_checked"] += 2
            counts["generic_orientation_support_new_failures"] += sum(
                case["ambient_dimension_passed"] and not case["passed"]
                for case in decision.generic_cases
            )
            counts["block_subsets_checked"] += decision.subsets_checked
            if not decision.passed:
                counts["bipartite_components_failed"] += 1
                failed_blocks.append(
                    compact_block_failure(instance, component, decision)
                )
        if failed_blocks:
            counts["z0_block_failed"] += 1
            z0_failures.append({
                "Z0": z0_absolute,
                "reason": "bipartite_block_support",
                "failed_components": failed_blocks,
            })
            continue
        counts["z0_block_passed"] += 1

        odd = [
            component.component
            for component in components
            if not component.bipartite
        ]
        raw_colorings = 1 if not odd else 1 << (len(odd) - 1)
        counts["pure_colorings_considered"] += raw_colorings
        coloring_failures = []
        for colors, bins in _pure_colorings(instance, z0, odd):
            counts["pure_colorings_matchable"] += 1
            counts["joint_support_searches"] += 1
            support = actual_supports_for_bins(instance, z0, bins)
            counts["joint_support_dfs_nodes"] += support.dfs_nodes
            if support.feasible:
                assert support.supports is not None
                return SameZ0SeedDecision(
                    True,
                    counts,
                    z0_absolute,
                    list(colors),
                    {
                        str(instance.outside[local]): [
                            instance.seed[coordinate]
                            for coordinate in vertices(actual)
                        ]
                        for local, actual in support.supports.items()
                    },
                    None,
                )
            coloring_failures.append({
                "component_colors": list(colors),
                "bins": [absolute_vertices(instance, bin_mask) for bin_mask in bins],
                "support_DFS_nodes": support.dfs_nodes,
            })

        counts["z0_nonbipartite_failed"] += 1
        z0_failures.append({
            "Z0": z0_absolute,
            "reason": "nonbipartite_joint_actual_support",
            "odd_components": [
                absolute_vertices(instance, component) for component in odd
            ],
            "raw_colorings": raw_colorings,
            "allowed_mask_matchable_colorings": len(coloring_failures),
            "coloring_failures": coloring_failures,
        })

    certificate = {
        "failure_kind": "no_single_Z0_satisfies_both_exact_K6_systems",
        "seed": list(instance.seed),
        "outside": list(instance.outside),
        "eligible_Z0": absolute_vertices(instance, instance.eligible_z0_mask),
        "allowed_defects": {
            str(instance.outside[local]): [
                instance.seed[coordinate] for coordinate in vertices(mask)
            ]
            for local, mask in enumerate(instance.defects)
        },
        "counts": counts.copy(),
        "Z0_failures": z0_failures,
        "same_Z0_cross_classification": rejected_seed_cross_classification(
            adj, instance
        ),
    }
    return SameZ0SeedDecision(False, counts, None, None, None, certificate)


def evaluate_graph(adj: Sequence[int], stop_on_rejection: bool = True) -> dict:
    """Reject iff one K6 seed has no common Z0 satisfying both systems."""

    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("same-Z0 production target must be K6-only")
    inertia = InertiaCache()
    totals = empty_counts()
    seeds_checked = 0
    impossible_seeds = 0
    first_seed = None
    first_certificate = None
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_seed(adj, instance, inertia)
        add_counts(totals, decision.counts)
        if not decision.feasible:
            impossible_seeds += 1
            if first_seed is None:
                first_seed = list(instance.seed)
                first_certificate = decision.certificate
            if stop_on_rejection:
                break
    return {
        "rejected": impossible_seeds > 0,
        "seeds_checked": seeds_checked,
        "impossible_seeds": impossible_seeds,
        **totals,
        "inertia_cache_entries": len(inertia.values),
        "inertia_cache_hits": inertia.hits,
        "first_impossible_seed": first_seed,
        "certificate": first_certificate,
    }


def evaluate_record(record: dict) -> dict:
    return {
        "index": record["index"],
        "decision": evaluate_graph(record["adjacency"], stop_on_rejection=True),
    }


def verify_inputs() -> tuple[dict[str, str], list[dict], list[int]]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(
            f"same-Z0 source/input hash mismatch: observed={observed}"
        )
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest.json").read_text(encoding="utf-8")
    )
    indices = manifest["classes"]["K6_only"]["final_indices"]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INDICES_SHA256:
        raise ValueError("same-Z0 K6 residue index set changed")
    payload = json.loads(
        (ROOT / "d6_k6_bipartite_rank_input.json").read_text(encoding="utf-8")
    )
    by_index = {record["index"]: record for record in payload["graphs"]}
    if any(index not in by_index for index in indices):
        raise ValueError("same-Z0 residue index absent from adjacency input")
    return observed, [by_index[index] for index in indices], indices


def run(output: Path, workers: int) -> dict:
    observed, records, indices = verify_inputs()
    started = time.perf_counter()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    wall = time.perf_counter() - started
    rejected = [item for item in results if item["decision"]["rejected"]]
    aggregate_fields = (
        "seeds_checked",
        "impossible_seeds",
        *COUNTER_FIELDS,
        "inertia_cache_entries",
        "inertia_cache_hits",
    )
    totals = {
        name: sum(item["decision"][name] for item in results)
        for name in aggregate_fields
    }

    positive = evaluate_graph(lower_bound_18_graph(), stop_on_rejection=True)
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point positive control failed")

    report = {
        "schema": SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact conjunction requiring one common Z0 for the bipartite "
            "normal block-support and nonbipartite joint-support systems."
        ),
        "input_graphs": len(results),
        "input_indices_sha256": stable_hash(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(results) - len(rejected),
        "rejected_indices": [item["index"] for item in rejected],
        **totals,
        "graph_results": results,
        "positive_control": {
            "name": "known realizable 18-point construction",
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
            "impossible_seeds": positive["impossible_seeds"],
        },
        "sources": {
            **observed,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": wall,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "quantifiers": {
            "seed": "a graph is rejected if one required K6 seed is impossible",
            "Z0": (
                "for that seed, every eligible Z0 fails at least one system; "
                "both systems are tested on the same Z0"
            ),
            "generic_orientation": (
                "a bipartite component passes if either generic orientation "
                "or the separate lightlike alternative passes"
            ),
            "nonbipartite_coloring": (
                "a Z0 passes if some symmetry-reduced light-ray coloring has "
                "a joint actual-support assignment"
            ),
        },
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "survivor": "filter non-rejection only, not a realizability claim",
            "arithmetic": "exact rational/integer decisions; no floating tolerance",
        },
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_same_z0_report.json"
    )
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    report = run(args.output.resolve(), args.workers)
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "rejected_indices": report["rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
