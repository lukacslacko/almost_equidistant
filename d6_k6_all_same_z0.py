#!/usr/bin/env python3
"""Exact conjunction of all current K6 systems on one common Z0.

For every required K6 seed, a realization has one actual zero-Lorentz-factor
set Z0.  This layer requires that same Z0 to satisfy all three independently
proved necessary systems:

1. the ordinary zero-forcing bipartite-component rank bound;
2. the exact normal-inertia block-support bound (including lightlike cases);
3. the non-bipartite two-light-ray joint actual-support CSP.

The search is exhaustive over eligible Z0 subsets of size at most six.  All
mathematical decisions are exact.  Candidate nonedges remain unconstrained,
and allowed defect coordinates may be zero.
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
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    _pure_colorings,
    check_bipartite_components,
    lorentz_components,
)
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
REPORT_SCHEMA = "d6-k6-all-same-z0-v1"
CERTIFICATE_SCHEMA = "d6-k6-all-same-z0-certificates-v1"
EXPECTED_INPUT = 990
EXPECTED_INDICES_SHA256 = (
    "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
)
EXPECTED_PRIOR_REJECTIONS = [652_900, 2_842_523]
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
    "d6_k6_same_z0.py": (
        "6d6dd948a245d5d8245f6ee10ce58c725aba29b9d3a9b325c8be11b59d519b5f"
    ),
    "d6_k6_same_z0_report.json": (
        "cea11ab95f24aa5fe79f06cc4d1a719b775ac5647551e0fb737b85835c443ca1"
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


def absolute_vertices(instance: K6LorentzInstance, mask: int) -> list[int]:
    return [instance.outside[local] for local in vertices(mask)]


COUNTER_FIELDS = (
    "z0_considered",
    "z0_matchable",
    "z0_zero_forcing_passed",
    "z0_zero_forcing_failed",
    "z0_block_passed_after_zero_forcing",
    "z0_block_failed_after_zero_forcing",
    "z0_nonbipartite_failed",
    "zero_forcing_components_checked",
    "block_components_checked",
    "block_components_failed",
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


def compact_block_failure(instance, component, decision) -> dict:
    return {
        "component": absolute_vertices(instance, component.component),
        "side_A": absolute_vertices(instance, component.side_a),
        "side_B": absolute_vertices(instance, component.side_b),
        "generic_cases": [
            {
                "case": case["case"],
                "rank_lower_A": case["rank_lower_A"],
                "rank_lower_B": case["rank_lower_B"],
                "ambient_dimension_lower": case["ambient_dimension_lower"],
                "ambient_dimension_passed": case["ambient_dimension_passed"],
                "passed": case["passed"],
                "first_failure": case["first_failure"],
            }
            for case in decision.generic_cases
        ],
        "lightlike_case": decision.lightlike_case,
    }


@dataclass
class NonbipartiteDecision:
    passed: bool
    raw_colorings: int
    matchable_colorings: int
    support_searches: int
    support_dfs_nodes: int
    odd_components: list[list[int]]
    chosen_colors: list[int] | None
    chosen_supports: dict[str, list[int]] | None
    coloring_failures: list[dict]


def check_nonbipartite_system(
    instance: K6LorentzInstance, z0: int, components
) -> NonbipartiteDecision:
    odd = [
        component.component for component in components if not component.bipartite
    ]
    odd_absolute = [absolute_vertices(instance, component) for component in odd]
    raw = 1 if not odd else 1 << (len(odd) - 1)
    matchable = 0
    searches = 0
    nodes = 0
    failures = []
    for colors, bins in _pure_colorings(instance, z0, odd):
        matchable += 1
        searches += 1
        support = actual_supports_for_bins(instance, z0, bins)
        nodes += support.dfs_nodes
        if support.feasible:
            assert support.supports is not None
            chosen = {
                str(instance.outside[local]): [
                    instance.seed[coordinate]
                    for coordinate in vertices(actual)
                ]
                for local, actual in support.supports.items()
            }
            return NonbipartiteDecision(
                True,
                raw,
                matchable,
                searches,
                nodes,
                odd_absolute,
                list(colors),
                chosen,
                failures,
            )
        failures.append({
            "component_colors": list(colors),
            "bins": [absolute_vertices(instance, bin_mask) for bin_mask in bins],
            "support_DFS_nodes": support.dfs_nodes,
        })
    return NonbipartiteDecision(
        False,
        raw,
        matchable,
        searches,
        nodes,
        odd_absolute,
        None,
        None,
        failures,
    )


def block_system(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    components,
    inertia: InertiaCache,
) -> tuple[bool, list[dict], dict[str, int]]:
    failures = []
    counts = {
        "components": 0,
        "generic": 0,
        "new_failures": 0,
        "subsets": 0,
    }
    for component in components:
        if not component.bipartite:
            continue
        decision = check_component(adj, instance, z0, component, inertia)
        counts["components"] += 1
        counts["generic"] += 2
        counts["new_failures"] += sum(
            case["ambient_dimension_passed"] and not case["passed"]
            for case in decision.generic_cases
        )
        counts["subsets"] += decision.subsets_checked
        if not decision.passed:
            failures.append(compact_block_failure(instance, component, decision))
    return not failures, failures, counts


def passed_category(block: bool, zero_forcing: bool, nonbipartite: bool) -> str:
    names = []
    if block:
        names.append("block")
    if zero_forcing:
        names.append("zero_forcing")
    if nonbipartite:
        names.append("nonbipartite")
    return "+".join(names) if names else "none"


def exhaustive_seed_certificate(
    adj: Sequence[int], instance: K6LorentzInstance
) -> dict:
    """Evaluate all three systems independently on every eligible Z0."""

    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    choices = []
    histogram: dict[str, int] = {}
    representatives: dict[str, list[int]] = {}
    for z0 in z0_subsets(instance.eligible_z0_mask):
        z0_absolute = absolute_vertices(instance, z0)
        matchable = support_matching(z0, instance.defects) is not None
        if matchable:
            local_components = lorentz_components(instance, z0)
            rank_passed, rank_checks = check_bipartite_components(
                adj, instance, z0, local_components, zero_forcing
            )
            block_passed, block_failures, _ = block_system(
                adj, instance, z0, local_components, inertia
            )
            nonbipartite = check_nonbipartite_system(
                instance, z0, local_components
            )
        else:
            rank_passed = False
            rank_checks = []
            block_passed = False
            block_failures = []
            nonbipartite = NonbipartiteDecision(
                False, 0, 0, 0, 0, [], None, None, []
            )
        category = passed_category(
            block_passed, rank_passed, nonbipartite.passed
        )
        histogram[category] = histogram.get(category, 0) + 1
        representatives.setdefault(category, z0_absolute)
        choices.append({
            "Z0": z0_absolute,
            "Z0_allowed_mask_matchable": matchable,
            "block_support_passed": block_passed,
            "zero_forcing_rank_passed": rank_passed,
            "nonbipartite_joint_support_passed": nonbipartite.passed,
            "category": category,
            "failed_block_components": block_failures,
            "failed_zero_forcing_components": [
                check.jsonable() for check in rank_checks if not check.passed
            ],
            "nonbipartite": {
                "odd_components": nonbipartite.odd_components,
                "raw_colorings": nonbipartite.raw_colorings,
                "allowed_mask_matchable_colorings": (
                    nonbipartite.matchable_colorings
                ),
                "support_searches": nonbipartite.support_searches,
                "support_DFS_nodes": nonbipartite.support_dfs_nodes,
                "chosen_component_colors": nonbipartite.chosen_colors,
                "chosen_actual_supports": nonbipartite.chosen_supports,
                "coloring_failures": nonbipartite.coloring_failures,
            },
        })
    return {
        "seed": list(instance.seed),
        "outside": list(instance.outside),
        "eligible_Z0": absolute_vertices(instance, instance.eligible_z0_mask),
        "allowed_defects": {
            str(instance.outside[local]): [
                instance.seed[coordinate] for coordinate in vertices(mask)
            ]
            for local, mask in enumerate(instance.defects)
        },
        "histogram": histogram,
        "representative_Z0": representatives,
        "all_three_pass_count": histogram.get(
            "block+zero_forcing+nonbipartite", 0
        ),
        "choices": choices,
    }


@dataclass
class AllSameZ0SeedDecision:
    feasible: bool
    counts: dict[str, int]
    chosen_z0: list[int] | None
    certificate: dict | None


def solve_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> AllSameZ0SeedDecision:
    counts = empty_counts()
    for z0 in z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            continue
        counts["z0_matchable"] += 1
        local_components = lorentz_components(instance, z0)

        rank_passed, rank_checks = check_bipartite_components(
            adj, instance, z0, local_components, zero_forcing
        )
        counts["zero_forcing_components_checked"] += len(rank_checks)
        if not rank_passed:
            counts["z0_zero_forcing_failed"] += 1
            continue
        counts["z0_zero_forcing_passed"] += 1

        blocks_passed, block_failures, block_counts = block_system(
            adj, instance, z0, local_components, inertia
        )
        counts["block_components_checked"] += block_counts["components"]
        counts["block_components_failed"] += len(block_failures)
        counts["generic_orientation_cases_checked"] += block_counts["generic"]
        counts["generic_orientation_support_new_failures"] += block_counts[
            "new_failures"
        ]
        counts["block_subsets_checked"] += block_counts["subsets"]
        if not blocks_passed:
            counts["z0_block_failed_after_zero_forcing"] += 1
            continue
        counts["z0_block_passed_after_zero_forcing"] += 1

        nonbipartite = check_nonbipartite_system(instance, z0, local_components)
        counts["pure_colorings_considered"] += nonbipartite.raw_colorings
        counts["pure_colorings_matchable"] += nonbipartite.matchable_colorings
        counts["joint_support_searches"] += nonbipartite.support_searches
        counts["joint_support_dfs_nodes"] += nonbipartite.support_dfs_nodes
        if nonbipartite.passed:
            return AllSameZ0SeedDecision(
                True, counts, absolute_vertices(instance, z0), None
            )
        counts["z0_nonbipartite_failed"] += 1

    certificate = exhaustive_seed_certificate(adj, instance)
    if certificate["all_three_pass_count"]:
        raise AssertionError("main search and exhaustive certificate disagree")
    return AllSameZ0SeedDecision(False, counts, None, certificate)


def evaluate_graph(adj: Sequence[int], stop_on_rejection: bool = True) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("all-same-Z0 production target must be K6-only")
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    totals = empty_counts()
    seeds_checked = 0
    impossible_seeds = 0
    first_seed = None
    first_certificate = None
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_seed(adj, instance, inertia, zero_forcing)
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
        "zero_forcing_cache_entries": len(zero_forcing.cache),
        "zero_forcing_cache_hits": zero_forcing.cache_hits,
        "zero_forcing_initial_sets_checked": zero_forcing.initial_sets_checked,
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
        raise ValueError(f"all-same-Z0 dependency hash mismatch: {observed}")
    prior = json.loads(
        (ROOT / "d6_k6_same_z0_report.json").read_text(encoding="utf-8")
    )
    if prior.get("rejected_indices") != EXPECTED_PRIOR_REJECTIONS:
        raise ValueError("prior two-index same-Z0 result changed")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest.json").read_text(encoding="utf-8")
    )
    indices = manifest["classes"]["K6_only"]["final_indices"]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INDICES_SHA256:
        raise ValueError("all-same-Z0 K6 residue selection changed")
    payload = json.loads(
        (ROOT / "d6_k6_bipartite_rank_input.json").read_text(encoding="utf-8")
    )
    by_index = {record["index"]: record for record in payload["graphs"]}
    if any(index not in by_index for index in indices):
        raise ValueError("all-same-Z0 selected index absent from adjacency input")
    return observed, [by_index[index] for index in indices], indices


def run(report_path: Path, certificate_path: Path, workers: int) -> dict:
    observed, records, indices = verify_inputs()
    started = time.perf_counter()
    if workers == 1:
        raw_results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            raw_results = list(pool.map(evaluate_record, records, chunksize=1))
    wall = time.perf_counter() - started
    rejected = [item for item in raw_results if item["decision"]["rejected"]]

    certificate_archive = {
        "schema": CERTIFICATE_SCHEMA,
        "status": "COMPLETE",
        "input_indices_sha256": stable_hash(indices),
        "production_source_sha256": sha256(Path(__file__).resolve()),
        "rejected_indices": [item["index"] for item in rejected],
        "certificates": [
            {
                "index": item["index"],
                "certificate": item["decision"]["certificate"],
            }
            for item in rejected
        ],
    }
    atomic_json(certificate_path, certificate_archive)

    graph_results = []
    for item in raw_results:
        decision = item["decision"].copy()
        decision.pop("certificate")
        graph_results.append({"index": item["index"], "decision": decision})

    aggregate_fields = (
        "seeds_checked",
        "impossible_seeds",
        *COUNTER_FIELDS,
        "inertia_cache_entries",
        "inertia_cache_hits",
        "zero_forcing_cache_entries",
        "zero_forcing_cache_hits",
        "zero_forcing_initial_sets_checked",
    )
    totals = {
        name: sum(item["decision"][name] for item in raw_results)
        for name in aggregate_fields
    }

    positive = evaluate_graph(lower_bound_18_graph(), stop_on_rejection=True)
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point positive control failed")

    rejected_indices = [item["index"] for item in rejected]
    if not set(EXPECTED_PRIOR_REJECTIONS).issubset(rejected_indices):
        raise AssertionError("all-same-Z0 result lost a prior exact rejection")
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact same-Z0 conjunction of ordinary zero-forcing rank, normal "
            "block support, and nonbipartite joint actual support."
        ),
        "input_graphs": len(graph_results),
        "input_indices_sha256": stable_hash(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(graph_results) - len(rejected),
        "rejected_indices": rejected_indices,
        "prior_same_Z0_rejected_indices": EXPECTED_PRIOR_REJECTIONS,
        "new_rejected_indices": sorted(
            set(rejected_indices) - set(EXPECTED_PRIOR_REJECTIONS)
        ),
        **totals,
        "graph_results": graph_results,
        "certificate_archive": {
            "path": certificate_path.name,
            "sha256": sha256(certificate_path),
            "bytes": certificate_path.stat().st_size,
            "certificates": len(rejected),
        },
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
            "seed": "one impossible required K6 seed rejects a graph",
            "Z0": (
                "every eligible size-at-most-six Z0 fails at least one of "
                "three systems; all three are evaluated on the same Z0"
            ),
            "zero_forcing": (
                "every bipartite component must satisfy the ordinary "
                "zero-forcing span lower bound"
            ),
            "block": (
                "every bipartite component must pass a generic orientation "
                "or the separate lightlike alternative"
            ),
            "nonbipartite": (
                "some light-ray coloring must have one joint support assignment"
            ),
        },
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "survivor": "filter non-rejection only, not realizability",
            "arithmetic": "exact integer/rational decisions; no tolerance",
        },
    }
    atomic_json(report_path, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_all_same_z0_report.json"
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=ROOT / "d6_k6_all_same_z0_certificates.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    report = run(
        args.output.resolve(), args.certificates.resolve(), args.workers
    )
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "rejected_indices": report["rejected_indices"],
        "new_rejected_indices": report["new_rejected_indices"],
        "output": str(args.output.resolve()),
        "certificates": str(args.certificates.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
