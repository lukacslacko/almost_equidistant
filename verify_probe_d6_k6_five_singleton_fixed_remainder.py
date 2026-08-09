#!/usr/bin/env python3
"""Independent verifier for the five-singleton fixed-remainder census.

The checker imports neither the producer nor its helper functions.  It starts
from the frozen independent saturated-singleton verifier, replays its distinct
actual-support DFS, and evaluates every fixed root with exact SymPy radical
arithmetic.  The producer uses the derived symmetric-difference/sign rule;
this checker substitutes directly in the K6 pair polynomial.
"""

from __future__ import annotations

import argparse
import ast
import itertools
import json
import os
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import sympy as sp

import verify_d6_k6_saturated_singleton_basis as basev


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_saturated_singleton_basis_report.json"
EXPECTED_PARENT_REPORT_SHA256 = (
    "5f25193311ad769a402ec5245a6a690619cc41225e63a56546fe50a30d7a93e5"
)
EXPECTED_INPUT_SHA256 = (
    "f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31"
)
EXPECTED_REJECTED_INDICES = [
    428_414, 430_135, 1_688_505, 2_669_917, 2_802_776,
    3_655_397, 3_931_365, 3_950_363, 3_968_088,
]
EXPECTED_REJECTED_SHA256 = (
    "b7de35407737cc884b7c227e880654c8d738f75889678d3f0d842e5f3696a452"
)
EXPECTED_RESIDUE_SHA256 = (
    "7b9c2e26e2ac178ab4e179dc0623b20f746bbd8eba85afb5c9c2c2013df8956e"
)
EXPECTED_GRAPH_ROWS_SHA256 = (
    "51a3fb891230cf26808deaf634d2ac823f23c251f93ebec0560423944936719d"
)
EXPECTED_POSITIVE_SEED_ROWS_SHA256 = (
    "c340fd60403f18ef8f0f9a0353c4edfa342d637a7581d6466a310ec06ccebb9c"
)
EXPECTED_TOTALS = {
    "complete_parent_leaves_reached": 8770,
    "complete_parent_leaves_older_six_basis_failed": 278,
    "complete_parent_leaves_z0_nonempty": 33,
    "complete_parent_leaves_z0_empty": 8459,
    "complete_parent_leaves_with_five_basis": 1024,
    "five_basis_extensions_total": 1024,
    "five_basis_extensions_not_applicable": 0,
    "five_basis_extensions_passing": 0,
    "five_basis_extensions_failing": 1024,
    "fixed_sign_branches_total": 17304,
    "fixed_sign_branches_passing": 0,
    "fixed_sign_branches_failing": 17304,
    "fixed_sign_branch_fixed_support_not_allowed": 1816,
    "fixed_sign_branch_assigned_support_mismatch": 0,
    "fixed_sign_branch_required_pair_wrong_symmetric_difference": 15116,
    "fixed_sign_branch_required_pair_equal_sign": 372,
    "fixed_sign_branch_collision": 0,
    "fixed_sign_branch_relaxed_free_points_ignored": 0,
}
DEFAULT_PRODUCTION = Path(
    "/private/tmp/d6_k6_five_singleton_fixed_remainder_production.json"
)
DEFAULT_OUTPUT = Path(
    "/private/tmp/d6_k6_five_singleton_fixed_remainder_verification.json"
)
ORIGINAL_LEAF_RULE = basev.saturated_singleton_leaves_pass
ORIGINAL_SOLVE_SEED = basev.solve_seed
R = sp.sqrt(7)

COUNTER_KEYS = (
    "complete_parent_leaves_reached",
    "complete_parent_leaves_older_six_basis_failed",
    "complete_parent_leaves_z0_nonempty",
    "complete_parent_leaves_z0_empty",
    "complete_parent_leaves_with_five_basis",
    "five_basis_extensions_total",
    "five_basis_extensions_not_applicable",
    "five_basis_extensions_passing",
    "five_basis_extensions_failing",
    "fixed_sign_branches_total",
    "fixed_sign_branches_passing",
    "fixed_sign_branches_failing",
    "fixed_sign_branch_fixed_support_not_allowed",
    "fixed_sign_branch_assigned_support_mismatch",
    "fixed_sign_branch_required_pair_wrong_symmetric_difference",
    "fixed_sign_branch_required_pair_equal_sign",
    "fixed_sign_branch_collision",
    "fixed_sign_branch_relaxed_free_points_ignored",
)

LOCAL_COUNTS: Counter[str] = Counter()
LOCAL_EXTENSIONS: list[dict] = []
LOCAL_SEEDS: list[dict] = []


@dataclass(frozen=True)
class FixedPoint:
    local: int
    absolute: int
    neighbourhood: int
    allowed: int
    assigned: int | None
    domain_failure: str | None


def canonical(value):
    return sp.expand(sp.radsimp(sp.cancel(value)))


def fixed_root(sign_bit: bool):
    epsilon = 1 if sign_bit else -1
    return canonical((1 + epsilon * R) / 3)


def independent_branch_reason(
    fixed: Sequence[FixedPoint], signs: int, adjacency: Sequence[int],
) -> tuple[str, list[int] | None]:
    for point in fixed:
        if point.domain_failure is not None:
            return point.domain_failure, [point.absolute]
    for left_number, right_number in itertools.combinations(range(len(fixed)), 2):
        left = fixed[left_number]
        right = fixed[right_number]
        t = fixed_root(bool(signs & (1 << left_number)))
        s = fixed_root(bool(signs & (1 << right_number)))
        common = (left.neighbourhood & right.neighbourhood).bit_count()
        coefficient = (
            left.neighbourhood.bit_count()
            + right.neighbourhood.bit_count() - 3 - 2 * common
        )
        required = bool(adjacency[left.absolute] & (1 << right.absolute))
        if required and canonical(t + s + coefficient * t * s) != 0:
            # Classify the independently found nonzero polynomial only for
            # branch-account comparison with the producer.
            if (left.neighbourhood ^ right.neighbourhood).bit_count() != 4:
                return (
                    "required_pair_wrong_symmetric_difference",
                    [left.absolute, right.absolute],
                )
            return "required_pair_equal_sign", [left.absolute, right.absolute]
        if left.neighbourhood == right.neighbourhood and canonical(t - s) == 0:
            return "collision", [left.absolute, right.absolute]
    return "relaxed_free_points_ignored", None


def independent_extension(
    adjacency: Sequence[int], instance, basis_locals: Sequence[int],
    assignments: dict[int, int],
) -> tuple[bool, dict]:
    if len(basis_locals) != 5:
        raise ValueError("independent five-basis needs five vertices")
    coordinate_of = {}
    for local in basis_locals:
        support = int(assignments[local])
        if support.bit_count() != 1:
            return True, {"applicable": False, "reason": "basis_not_singleton"}
        coordinate = support.bit_length() - 1
        if coordinate in coordinate_of.values():
            return False, {
                "applicable": True,
                "reason": "duplicate_singleton_basis_coordinate",
                "branch_outcomes": [],
            }
        coordinate_of[local] = coordinate
    used = sum(1 << coordinate for coordinate in coordinate_of.values())
    missing = 0b111111 ^ used
    if used.bit_count() != 5 or missing.bit_count() != 1:
        return False, {
            "applicable": True,
            "reason": "basis_does_not_span_five_coordinates",
            "branch_outcomes": [],
        }
    fixed = []
    free_vertices = []
    for local in range(len(instance.outside)):
        if local in basis_locals:
            continue
        absolute = instance.outside[local]
        neighbourhood = 0
        for basis_local in basis_locals:
            if adjacency[absolute] & (1 << instance.outside[basis_local]):
                neighbourhood |= 1 << coordinate_of[basis_local]
        allowed = int(instance.defects[local])
        assigned = assignments.get(local)
        uncontrolled = allowed & used & ~neighbourhood
        if uncontrolled:
            return True, {
                "applicable": False,
                "reason": "basis_nonedge_coordinate_not_forced_zero",
                "branch_outcomes": [],
            }
        fixed_missing = not bool(allowed & missing) or bool(
            assigned is not None and not (assigned & missing)
        )
        if not fixed_missing:
            free_vertices.append(absolute)
            continue
        failure = None
        if neighbourhood & ~allowed:
            failure = "fixed_support_not_allowed"
        elif assigned is not None and assigned != neighbourhood:
            failure = "assigned_support_mismatch"
        fixed.append(FixedPoint(
            local, absolute, neighbourhood, allowed, assigned, failure
        ))

    counts = Counter()
    outcomes = []
    for signs in range(1 << len(fixed)):
        counts["fixed_sign_branches_total"] += 1
        reason, vertices = independent_branch_reason(fixed, signs, adjacency)
        counts[f"fixed_sign_branch_{reason}"] += 1
        passed = reason == "relaxed_free_points_ignored"
        counts[
            "fixed_sign_branches_passing" if passed
            else "fixed_sign_branches_failing"
        ] += 1
        outcomes.append({
            "fixed_sign_bits": signs,
            "passed": passed,
            "reason": reason,
            "vertices": vertices,
        })
    feasible = counts["fixed_sign_branches_passing"] > 0
    if (
        len(outcomes) != counts["fixed_sign_branches_total"]
        or counts["fixed_sign_branches_total"]
        != counts["fixed_sign_branches_passing"]
        + counts["fixed_sign_branches_failing"]
    ):
        raise AssertionError("independent fixed-sign accounting differs")
    return feasible, {
        "applicable": True,
        "reason": "some_fixed_sign_branch_passes" if feasible
        else "all_fixed_sign_branches_fail",
        "basis_vertices": [instance.outside[local] for local in basis_locals],
        "missing_coordinate": missing.bit_length() - 1,
        "fixed_points": [{
            "vertex": point.absolute,
            "basis_neighbourhood": list(basev.base.bits(point.neighbourhood)),
            "allowed_seed_defects": list(basev.base.bits(point.allowed)),
            "assigned_actual_support": None if point.assigned is None
            else list(basev.base.bits(point.assigned)),
            "domain_failure": point.domain_failure,
        } for point in fixed],
        "free_vertices_ignored": free_vertices,
        "counts": dict(counts),
        "branch_outcomes": outcomes,
    }


def strengthened_leaf_rule(adjacency, instance, z0, bins, assignments):
    LOCAL_COUNTS["complete_parent_leaves_reached"] += 1
    if not ORIGINAL_LEAF_RULE(adjacency, instance, z0, bins, assignments):
        LOCAL_COUNTS["complete_parent_leaves_older_six_basis_failed"] += 1
        return False
    if z0:
        LOCAL_COUNTS["complete_parent_leaves_z0_nonempty"] += 1
        return True
    LOCAL_COUNTS["complete_parent_leaves_z0_empty"] += 1
    leaf_has_basis = False
    leaf_passes = True
    for bin_mask in bins:
        singleton_locals = tuple(
            local for local in basev.base.bits(bin_mask)
            if assignments[local].bit_count() == 1
        )
        for basis_locals in itertools.combinations(singleton_locals, 5):
            leaf_has_basis = True
            feasible, detail = independent_extension(
                adjacency, instance, basis_locals, assignments
            )
            LOCAL_COUNTS["five_basis_extensions_total"] += 1
            if not detail.get("applicable"):
                LOCAL_COUNTS["five_basis_extensions_not_applicable"] += 1
            LOCAL_COUNTS[
                "five_basis_extensions_passing" if feasible
                else "five_basis_extensions_failing"
            ] += 1
            LOCAL_COUNTS.update(detail.get("counts", {}))
            LOCAL_EXTENSIONS.append(detail)
            leaf_passes &= feasible
    if leaf_has_basis:
        LOCAL_COUNTS["complete_parent_leaves_with_five_basis"] += 1
    return leaf_passes


def instrumented_solve_seed(adjacency, instance):
    before = Counter(LOCAL_COUNTS)
    extension_start = len(LOCAL_EXTENSIONS)
    decision = ORIGINAL_SOLVE_SEED(adjacency, instance)
    delta = Counter(LOCAL_COUNTS)
    delta.subtract(before)
    LOCAL_SEEDS.append({
        "seed": list(instance.seed),
        "seed_mask": sum(1 << vertex for vertex in instance.seed),
        "feasible": bool(decision["feasible"]),
        "counts": {key: delta[key] for key in COUNTER_KEYS},
        "extension_start": extension_start,
        "extension_stop": len(LOCAL_EXTENSIONS),
    })
    return decision


def worker_initializer():
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    basev.worker_initializer()
    basev.saturated_singleton_leaves_pass = strengthened_leaf_rule
    basev.solve_seed = instrumented_solve_seed


def evaluate_record(record):
    global LOCAL_COUNTS, LOCAL_EXTENSIONS, LOCAL_SEEDS
    LOCAL_COUNTS = Counter()
    LOCAL_EXTENSIONS = []
    LOCAL_SEEDS = []
    row = basev.evaluate_record(record)
    if len(LOCAL_SEEDS) != row["seeds_checked"]:
        raise AssertionError("independent seed instrumentation differs")
    row["fixed_remainder_counts"] = {
        key: LOCAL_COUNTS[key] for key in COUNTER_KEYS
    }
    row["fixed_remainder_seeds"] = list(LOCAL_SEEDS)
    row["fixed_remainder_extensions"] = list(LOCAL_EXTENSIONS)
    return row


def assert_import_independence():
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imports.update(
        node.module for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    )
    if "probe_d6_k6_five_singleton_fixed_remainder" in imports:
        raise AssertionError("independent checker imports producer")


def load_records():
    if basev.sha256(PARENT_REPORT) != EXPECTED_PARENT_REPORT_SHA256:
        raise ValueError("certified singleton-basis report boundary changed")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    indices = list(map(int, report["ordered_residue_indices"]))
    if basev.stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("certified 249-graph residue changed")
    records, parent_indices = basev.load_input()
    by_index = {int(record["index"]): record for record in records}
    if not set(indices) <= set(parent_indices):
        raise ValueError("249-graph residue is outside independent input")
    return [by_index[index] for index in indices], indices


def run_census(records, workers):
    if workers == 1:
        worker_initializer()
        return list(map(evaluate_record, records))
    with ProcessPoolExecutor(
        max_workers=workers, initializer=worker_initializer
    ) as pool:
        return list(pool.map(evaluate_record, records, chunksize=1))


def totals(results):
    aggregate = Counter()
    for row in results:
        aggregate.update(row["fixed_remainder_counts"])
    result = {key: aggregate[key] for key in COUNTER_KEYS}
    if (
        result["five_basis_extensions_total"]
        != result["five_basis_extensions_passing"]
        + result["five_basis_extensions_failing"]
        or result["fixed_sign_branches_total"]
        != result["fixed_sign_branches_passing"]
        + result["fixed_sign_branches_failing"]
    ):
        raise AssertionError("independent aggregate accounting differs")
    return result


def validate_production(production, indices):
    graph_rows = production.get("graph_rows")
    positive_rows = production.get("positive_seed_rows")
    if (
        production.get("schema")
        != "d6-k6-five-singleton-fixed-remainder-pilot-v1"
        or production.get("status") != "COMPLETE_DISCOVERY_ONLY"
        or production.get("input_graphs") != 249
        or production.get("ordered_input_indices_sha256")
        != EXPECTED_INPUT_SHA256
        or production.get("marginal_indices") != EXPECTED_REJECTED_INDICES
        or production.get("marginal_indices_sha256")
        != EXPECTED_REJECTED_SHA256
        or production.get("ordered_residue_indices_sha256")
        != EXPECTED_RESIDUE_SHA256
        or not isinstance(graph_rows, list)
        or len(graph_rows) != len(indices)
        or [row.get("index") for row in graph_rows] != indices
        or production.get("graph_rows_sha256") != basev.stable_hash(graph_rows)
        or not isinstance(positive_rows, list)
        or len(positive_rows) != 32
        or production.get("positive_seed_rows_sha256")
        != basev.stable_hash(positive_rows)
    ):
        raise ValueError("production census boundary is malformed")
    for row in graph_rows:
        extensions = row.get("extensions")
        seeds = row.get("seeds")
        if (
            not isinstance(extensions, list)
            or not isinstance(seeds, list)
            or len(seeds) != row.get("seeds_checked")
            or sum(
                len(extension.get("branch_outcomes", []))
                for extension in extensions if extension.get("applicable")
            ) != row.get("counts", {}).get("fixed_sign_branches_total")
        ):
            raise ValueError("production graph branch/seed coverage differs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--production", type=Path, default=DEFAULT_PRODUCTION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    assert_import_independence()
    records, indices = load_records()
    production = json.loads(args.production.read_text(encoding="utf-8"))
    validate_production(production, indices)
    started = time.monotonic()
    results = run_census(records, args.workers)
    rejected = [int(row["index"]) for row in results if row["rejected"]]
    residue = [int(row["index"]) for row in results if not row["rejected"]]
    if (
        rejected != EXPECTED_REJECTED_INDICES
        or basev.stable_hash(rejected) != EXPECTED_REJECTED_SHA256
        or basev.stable_hash(residue) != EXPECTED_RESIDUE_SHA256
        or rejected != production["marginal_indices"]
        or residue != production["ordered_residue_indices"]
    ):
        raise AssertionError("independent graph decisions differ")
    graph_rows = [{
        "index": int(row["index"]),
        "rejected": bool(row["rejected"]),
        "seeds_checked": int(row["seeds_checked"]),
        "seed_mask": row["seed_mask"],
        "counts": row["fixed_remainder_counts"],
        "seeds": row["fixed_remainder_seeds"],
        "extensions": row["fixed_remainder_extensions"],
    } for row in results]
    independent_totals = totals(results)
    if (
        independent_totals != EXPECTED_TOTALS
        or basev.stable_hash(graph_rows) != EXPECTED_GRAPH_ROWS_SHA256
    ):
        raise AssertionError("independent exhaustive transcript changed")
    production_by_index = {
        int(row["index"]): row for row in production["graph_rows"]
    }
    for row in graph_rows:
        other = production_by_index[row["index"]]
        if (
            row["rejected"] != other["rejected"]
            or row["seeds_checked"] != other["seeds_checked"]
        ):
            raise AssertionError("per-graph seed/decision coverage differs")

    global LOCAL_COUNTS, LOCAL_EXTENSIONS, LOCAL_SEEDS
    LOCAL_COUNTS = Counter()
    LOCAL_EXTENSIONS = []
    LOCAL_SEEDS = []
    basev.saturated_singleton_leaves_pass = strengthened_leaf_rule
    basev.solve_seed = instrumented_solve_seed
    positive = basev.positive_control()
    if positive != {"passed": True, "K6_seeds": 32} or len(LOCAL_SEEDS) != 32:
        raise AssertionError("independent positive control differs")
    positive_rows = list(LOCAL_SEEDS)
    if basev.stable_hash(positive_rows) != EXPECTED_POSITIVE_SEED_ROWS_SHA256:
        raise AssertionError("independent positive-seed transcript changed")

    checks = {
        "producer_import_absent": True,
        "production_branch_coverage_valid": True,
        "ordered_rejections_agree": True,
        "ordered_residue_agrees": True,
        "per_graph_seed_counts_agree": True,
        "all_32_positive_K6_seeds_pass": True,
    }
    output = {
        "schema": "d6-k6-five-singleton-fixed-remainder-verification-v1",
        "status": "PASS",
        "checks": checks,
        "input_graphs": len(records),
        "ordered_input_indices_sha256": basev.stable_hash(indices),
        "marginal_rejections": len(rejected),
        "marginal_indices": rejected,
        "marginal_indices_sha256": basev.stable_hash(rejected),
        "residue": len(residue),
        "ordered_residue_indices_sha256": basev.stable_hash(residue),
        "independent_totals": independent_totals,
        "production_totals": production["totals"],
        "independent_graph_rows_sha256": basev.stable_hash(graph_rows),
        "production_graph_rows_sha256": production["graph_rows_sha256"],
        "positive_18_control": positive,
        "independent_positive_seed_rows_sha256": basev.stable_hash(positive_rows),
        "production_positive_seed_rows_sha256": (
            production["positive_seed_rows_sha256"]
        ),
        "runtime": {
            "workers": args.workers,
            "wall_seconds": time.monotonic() - started,
        },
        "sources": {
            Path(__file__).name: basev.sha256(Path(__file__)),
            Path(basev.__file__).name: basev.sha256(Path(basev.__file__)),
            PARENT_REPORT.name: basev.sha256(PARENT_REPORT),
            args.production.name: basev.sha256(args.production),
        },
        "arithmetic": "exact independently canonicalized SymPy radicals",
    }
    basev.atomic_json(args.output, output)
    print(json.dumps({
        "status": "PASS",
        "marginal_rejections": len(rejected),
        "residue": len(residue),
        "independent_totals": output["independent_totals"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
