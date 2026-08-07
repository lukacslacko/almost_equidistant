#!/usr/bin/env python3
"""Exact same-Z0 conjunction on the frozen 625-graph K6 residue.

A geometric realization has one actual zero-Lorentz-factor set Z0.  For each
required K6 seed this incremental layer requires the same enumerated Z0 to
pass both the frozen tight-Hall/singleton/repeated-arm bipartite system and
the proved nonbipartite two-light-ray actual-support CSP.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import d6_k6_repeated_two_support_arm as prior
import d6_k6_tight_hall_support as kernel
from d6_k6_bipartite_rank_reference import _pure_colorings
from d6_k6_support_reference import actual_supports_for_bins


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_repeated_two_support_arm_report.json"
PARENT_CERTIFICATES_GZ = (
    ROOT / "d6_k6_repeated_two_support_arm_certificates.json.gz"
)
PARENT_VERIFICATION = (
    ROOT / "d6_k6_repeated_two_support_arm_verification.json"
)
REPORT_SCHEMA = "d6-k6-tight-same-z0-v1"
CERTIFICATE_SCHEMA = "d6-k6-tight-same-z0-certificates-v1"
EXPECTED_INPUT = 625
EXPECTED_INPUT_SHA256 = (
    "04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1"
)
EXPECTED_PARENT_RAW_CERTIFICATE_SHA256 = (
    "fbaebe286acce8b0c3603f89f05a5316010b51dceabfeface7e73ff55ce90ed2"
)
DEFAULT_WORKERS = 8

EXPECTED_DEPENDENCIES = {
    "d6_k6_repeated_two_support_arm.py": (
        "cf4a0a534504e312c47fcf1974d3362f19beeec270e4a49f04eba3080d06e1d4"
    ),
    "d6_k6_repeated_two_support_arm_report.json": (
        "ea8d9438d062b92857a1057b950e76b8ae08bbe1bf97b4327130fd4058171ebe"
    ),
    "d6_k6_repeated_two_support_arm_certificates.json.gz": (
        "604ec86b943ec909d177acc1fd36ef623c3e3bff67b4a8a5e5ffa601fc202af5"
    ),
    "d6_k6_repeated_two_support_arm_verification.json": (
        "3abf48bab092671a9dff03c2dc47c0063a75515cd219e1d5a4bfa9a5fa797b5b"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_support_reference.py": (
        "d6481137ef49d88154882285660dd755eb7cb652ab87280f05c3792de9742577"
    ),
    "d6_k6_same_z0.py": (
        "6d6dd948a245d5d8245f6ee10ce58c725aba29b9d3a9b325c8be11b59d519b5f"
    ),
    "d6_k6_same_z0_report.json": (
        "cea11ab95f24aa5fe79f06cc4d1a719b775ac5647551e0fb737b85835c443ca1"
    ),
    "verify_d6_k6_same_z0.py": (
        "4db05ff06910a11395e4bbdd2fd13161d6c19f58f81cc99851b8c37558d7000f"
    ),
    "d6_k6_same_z0_verification.json": (
        "9d1f06c41827193f3a324bf2d08f450cf816567cbae727e21a6258afd68a9fbf"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gunzipped_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    prior.activate_kernel()


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"tight same-Z0 dependency boundary changed: {observed}")
    if gunzipped_sha256(PARENT_CERTIFICATES_GZ) != (
        EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
    ):
        raise ValueError("repeated-arm compressed certificate payload changed")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    if (
        report.get("schema") != "d6-k6-repeated-two-support-arm-v1"
        or report.get("status") != "COMPLETE"
        or report.get("graphs_rejected") != 9
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
        or verification.get("schema")
        != "d6-k6-repeated-two-support-arm-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 9
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen repeated-arm result boundary changed")
    all_records, _, _ = prior.load_input()
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("625-graph residue changed")
    return [by_index[index] for index in indices], indices, observed


def nonbipartite_support_witness(instance, z0: int, components) -> dict | None:
    odd = [component.component for component in components if not component.bipartite]
    for colors, bins in _pure_colorings(instance, z0, odd):
        support = actual_supports_for_bins(instance, z0, bins)
        if support.feasible:
            return {
                "component_colors": list(colors),
                "actual_supports": support.supports,
                "dfs_nodes": support.dfs_nodes,
            }
    return None


def bipartite_states(
    adjacency: Sequence[int], instance, z0: int,
    inertia, forcing, subset_rank,
) -> tuple[set[tuple[int, ...]], tuple, Counter[str]]:
    states: set[tuple[int, ...]] = {()}
    counts: Counter[str] = Counter()
    components = kernel.lorentz_components(instance, z0)
    for component in components:
        if not component.bipartite:
            counts["nonbipartite_components_deferred_to_conjunction"] += 1
            continue
        counts["bipartite_components"] += 1
        options, _, local = kernel.component_options(
            adjacency, instance, z0, component,
            inertia, forcing, subset_rank,
        )
        counts.update(local)
        states = kernel.extend_states(states, options, adjacency)
        if not states:
            break
    return states, components, counts


def full_cross_classification(
    adjacency: Sequence[int], instance,
) -> dict:
    """Evaluate both systems independently on every Z0 of a rejected seed."""

    inertia = kernel.InertiaCache()
    forcing = kernel.ZeroForcingSolver()
    subset_rank = kernel.arbitrary.PositiveSubsetRank(inertia, forcing)
    rows = []
    histogram: Counter[str] = Counter()
    for z0 in kernel.ranks.z0_subsets(instance.eligible_z0_mask):
        zvertices = [instance.outside[local] for local in kernel.vertices(z0)]
        if kernel.support_matching(z0, instance.defects) is None:
            bipartite_passed = False
            nonbipartite_passed = False
        else:
            states, components, _ = bipartite_states(
                adjacency, instance, z0, inertia, forcing, subset_rank
            )
            bipartite_passed = bool(states)
            nonbipartite_passed = (
                nonbipartite_support_witness(instance, z0, components) is not None
            )
        if bipartite_passed and nonbipartite_passed:
            category = "both"
        elif bipartite_passed:
            category = "bipartite_only"
        elif nonbipartite_passed:
            category = "nonbipartite_only"
        else:
            category = "neither"
        histogram[category] += 1
        rows.append({
            "Z0": zvertices,
            "tight_bipartite_passed": bipartite_passed,
            "nonbipartite_actual_support_passed": nonbipartite_passed,
            "category": category,
        })
    return {"histogram": dict(sorted(histogram.items())), "choices": rows}


def solve_seed(adjacency: Sequence[int], instance) -> dict:
    inertia = kernel.InertiaCache()
    forcing = kernel.ZeroForcingSolver()
    subset_rank = kernel.arbitrary.PositiveSubsetRank(inertia, forcing)
    counts: Counter[str] = Counter()
    failures = []
    for z0 in kernel.ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        zvertices = [instance.outside[local] for local in kernel.vertices(z0)]
        if kernel.support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            failures.append({"Z0": zvertices, "reason": "Z0_matching"})
            continue
        states, components, local = bipartite_states(
            adjacency, instance, z0, inertia, forcing, subset_rank
        )
        counts.update(local)
        if not states:
            counts["z0_bipartite_failed"] += 1
            failures.append({
                "Z0": zvertices, "reason": "tight_bipartite_system"
            })
            continue
        counts["z0_bipartite_passed"] += 1
        counts["nonbipartite_support_searches"] += 1
        nonbipartite = nonbipartite_support_witness(instance, z0, components)
        if nonbipartite is None:
            counts["z0_nonbipartite_failed"] += 1
            failures.append({
                "Z0": zvertices,
                "reason": "nonbipartite_two_light_ray_actual_support",
            })
            continue
        counts["z0_same_conjunction_passed"] += 1
        return {
            "feasible": True,
            "counts": dict(counts),
            "witness": {
                "Z0": zvertices,
                "empty_state": list(min(states, key=lambda state: (len(state), state))),
                "nonbipartite_component_colors": nonbipartite["component_colors"],
            },
            "failures": None,
        }
    return {
        "feasible": False,
        "counts": dict(counts),
        "witness": None,
        "failures": failures,
    }


def evaluate_record(record: dict) -> dict:
    prior.activate_kernel()
    adjacency = tuple(map(int, record["adjacency"]))
    total: Counter[str] = Counter()
    seeds_checked = 0
    for seed_mask in kernel.clique_masks(adjacency, kernel.COORDINATES):
        seeds_checked += 1
        instance = kernel.build_instance(adjacency, kernel.vertices(seed_mask))
        decision = solve_seed(adjacency, instance)
        total.update(decision["counts"])
        if not decision["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "first_impossible_seed_mask": seed_mask,
                "first_impossible_seed": list(instance.seed),
                "counts": dict(total),
                "certificate": {
                    "seed_mask": seed_mask,
                    "seed": list(instance.seed),
                    "Z0_failures": decision["failures"],
                    "same_Z0_cross_classification": full_cross_classification(
                        adjacency, instance
                    ),
                },
            }
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "first_impossible_seed_mask": None,
        "first_impossible_seed": None,
        "counts": dict(total),
        "certificate": None,
    }


def positive_control() -> dict:
    prior.activate_kernel()
    adjacency = kernel.lower_bound_18_graph()
    passed = 0
    for seed_mask in kernel.clique_masks(adjacency, kernel.COORDINATES):
        instance = kernel.build_instance(adjacency, kernel.vertices(seed_mask))
        if not solve_seed(adjacency, instance)["feasible"]:
            raise AssertionError("tight same-Z0 rejects positive 18 control")
        passed += 1
    if passed != 32:
        raise AssertionError("positive control K6 seed count changed")
    return {"passed": True, "K6_seeds": passed}


def run(workers: int, output: Path, certificates: Path) -> dict:
    records, indices, dependencies = load_input()
    started = time.monotonic()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [row["index"] for row in results if row["rejected"]]
    residue = [row["index"] for row in results if not row["rejected"]]
    counts: Counter[str] = Counter()
    for row in results:
        counts.update(row["counts"])
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "production_source_sha256": sha256(Path(__file__)),
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "rejected_graphs": [
            {"index": row["index"], **row["certificate"]}
            for row in results if row["rejected"]
        ],
    }
    atomic_json(certificates, archive)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": "same-Z0 conjunction of strongest bipartite and nonbipartite K6 systems",
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_actual_subsets_exhausted_in_nonbipartite_CSP",
            "same_Z0": "one_geometric_zero_factor_set_shared_by_both_systems",
            "arithmetic": "exact_integer_graph_bitmask_Hall_inertia_and_DFS",
        },
        "production_source_sha256": sha256(Path(__file__)),
        "dependencies": dependencies,
        "parent_raw_certificate_payload_sha256": (
            EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
        ),
        "ordered_input_indices": indices,
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "input_graphs": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices": residue,
        "ordered_residue_indices_sha256": stable_hash(residue),
        "totals": {
            "seeds_checked": sum(row["seeds_checked"] for row in results),
            **dict(counts),
        },
        "quantifiers": {
            "seed": (
                "a graph is rejected only if one required K6 seed has no "
                "surviving branch"
            ),
            "Z0": (
                "every eligible actual zero-factor set is enumerated, and "
                "both constituent systems must pass for that identical Z0"
            ),
            "constituent_witnesses": (
                "for one fixed Z0 the two necessary systems may choose "
                "different actual-support witnesses, a sound relaxation"
            ),
        },
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "rejected_graphs": len(archive["rejected_graphs"]),
        },
        "positive_18_control": positive_control(),
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "graph_results": [
            {key: value for key, value in row.items() if key != "certificate"}
            for row in results
        ],
        "nonclaims": [
            "a survivor is not a realization",
            "the two systems may still use different actual Z0 support subsets",
            "this layer does not settle the K6 or dimension-six problem",
        ],
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_tight_same_z0_report.json"
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_tight_same_z0_certificates.json",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    report = run(args.workers, args.output, args.certificates)
    print(json.dumps({
        "status": report["status"],
        "input": report["input_graphs"],
        "rejected": report["graphs_rejected"],
        "surviving": report["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
