#!/usr/bin/env python3
"""Exact source census for the five-singleton K6 fixed-remainder rule.

At a complete actual-support leaf of the certified K6 search, choose five
same-light-ray points with distinct singleton supports.  Any remaining point
whose sixth defect coordinate is forced to zero has only the two exact roots
``t=(1+epsilon*sqrt(7))/3``.  This probe exhausts every such sign and checks
only required edges plus mandatory distinctness.  The other, free-coordinate
points are ignored, making the rule a sound relaxation.

Candidate nonedges receive no distance condition.  Allowed defect entries
may vanish.  Output is discovery-only and defaults to /private/tmp.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import d6_k6_saturated_singleton_basis as base


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_saturated_singleton_basis_report.json"
EXPECTED_PARENT_REPORT_SHA256 = (
    "5f25193311ad769a402ec5245a6a690619cc41225e63a56546fe50a30d7a93e5"
)
EXPECTED_INPUT_SHA256 = (
    "f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31"
)
EXPECTED_REJECTED_INDICES = (
    428_414, 430_135, 1_688_505, 2_669_917, 2_802_776,
    3_655_397, 3_931_365, 3_950_363, 3_968_088,
)
EXPECTED_REJECTED_SHA256 = (
    "b7de35407737cc884b7c227e880654c8d738f75889678d3f0d842e5f3696a452"
)
EXPECTED_RESIDUE_SHA256 = (
    "7b9c2e26e2ac178ab4e179dc0623b20f746bbd8eba85afb5c9c2c2013df8956e"
)
EXPECTED_GRAPH_ROWS_SHA256 = (
    "f5297d9136b20cd9e307dbbfb8778798cd6d09cf5513dc174398a3e05300df02"
)
EXPECTED_POSITIVE_SEED_ROWS_SHA256 = (
    "c340fd60403f18ef8f0f9a0353c4edfa342d637a7581d6466a310ec06ccebb9c"
)
EXPECTED_TOTALS = {
    "complete_parent_leaves_reached": 8785,
    "complete_parent_leaves_older_six_basis_failed": 278,
    "complete_parent_leaves_z0_nonempty": 33,
    "complete_parent_leaves_z0_empty": 8474,
    "complete_parent_leaves_with_five_basis": 1039,
    "five_basis_extensions_total": 1039,
    "five_basis_extensions_not_applicable": 0,
    "five_basis_extensions_passing": 0,
    "five_basis_extensions_failing": 1039,
    "fixed_sign_branches_total": 16760,
    "fixed_sign_branches_passing": 0,
    "fixed_sign_branches_failing": 16760,
    "fixed_sign_branch_fixed_support_not_allowed": 2064,
    "fixed_sign_branch_assigned_support_mismatch": 0,
    "fixed_sign_branch_required_pair_wrong_symmetric_difference": 14300,
    "fixed_sign_branch_required_pair_equal_sign": 396,
    "fixed_sign_branch_collision": 0,
    "fixed_sign_branch_relaxed_free_points_ignored": 0,
}
DEFAULT_OUTPUT = Path(
    "/private/tmp/d6_k6_five_singleton_fixed_remainder_production.json"
)
ORIGINAL_LEAF_RULE = base.singleton_basis_leaf_passes
ORIGINAL_SOLVE_SEED = base.solve_seed

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


def branch_reason(
    fixed: Sequence[FixedPoint], signs: int, adjacency: Sequence[int],
) -> tuple[str, list[int] | None]:
    """Return the first exact failure, or the relaxed passing reason."""

    for number, point in enumerate(fixed):
        if point.domain_failure is not None:
            return point.domain_failure, [point.absolute]
    for left_number, right_number in itertools.combinations(range(len(fixed)), 2):
        left = fixed[left_number]
        right = fixed[right_number]
        left_sign = bool(signs & (1 << left_number))
        right_sign = bool(signs & (1 << right_number))
        same_sign = left_sign == right_sign
        required = bool(
            adjacency[left.absolute] & (1 << right.absolute)
        )
        if required:
            symmetric_difference = (
                left.neighbourhood ^ right.neighbourhood
            ).bit_count()
            if symmetric_difference != 4:
                return (
                    "required_pair_wrong_symmetric_difference",
                    [left.absolute, right.absolute],
                )
            if same_sign:
                return "required_pair_equal_sign", [left.absolute, right.absolute]
        # Equal neighbourhood plus equal root reconstructs the same point.
        # This is a distinctness condition, not a nonedge distance condition.
        if left.neighbourhood == right.neighbourhood and same_sign:
            return "collision", [left.absolute, right.absolute]
    return "relaxed_free_points_ignored", None


def fixed_remainder_extension(
    adjacency: Sequence[int], instance, basis_locals: Sequence[int],
    assignments: dict[int, int],
) -> tuple[bool, dict]:
    """Exhaust all fixed-point roots for one five-singleton basis."""

    if len(basis_locals) != 5:
        raise ValueError("five singleton basis vertices are required")
    coordinate_of: dict[int, int] = {}
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

    basis_set = set(basis_locals)
    fixed: list[FixedPoint] = []
    free_vertices = []
    for local in range(len(instance.outside)):
        if local in basis_set:
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
            # This is absent in the alpha<=2 corpus.  Accept rather than use
            # zero coordinates not justified by the graph boundary.
            return True, {
                "applicable": False,
                "reason": "basis_nonedge_coordinate_not_forced_zero",
                "vertex": absolute,
                "uncontrolled_coordinates": (
                    base.parent.parent.kernel.vertices(uncontrolled)
                ),
                "branch_outcomes": [],
            }
        forced_zero_missing = not bool(allowed & missing) or bool(
            assigned is not None and not (assigned & missing)
        )
        if not forced_zero_missing:
            free_vertices.append(absolute)
            continue
        domain_failure = None
        if neighbourhood & ~allowed:
            # Fixed v=0 gives nonzero t, so every basis neighbour is in the
            # actual support and must be allowed.
            domain_failure = "fixed_support_not_allowed"
        elif assigned is not None and assigned != neighbourhood:
            domain_failure = "assigned_support_mismatch"
        fixed.append(FixedPoint(
            local=local,
            absolute=absolute,
            neighbourhood=neighbourhood,
            allowed=allowed,
            assigned=assigned,
            domain_failure=domain_failure,
        ))

    counts: Counter[str] = Counter()
    outcomes = []
    for signs in range(1 << len(fixed)):
        counts["fixed_sign_branches_total"] += 1
        reason, pair = branch_reason(fixed, signs, adjacency)
        counts[f"fixed_sign_branch_{reason}"] += 1
        passing = reason == "relaxed_free_points_ignored"
        counts[
            "fixed_sign_branches_passing" if passing
            else "fixed_sign_branches_failing"
        ] += 1
        outcomes.append({
            "fixed_sign_bits": signs,
            "passed": passing,
            "reason": reason,
            "vertices": pair,
        })
    feasible = counts["fixed_sign_branches_passing"] > 0
    if (
        counts["fixed_sign_branches_total"]
        != counts["fixed_sign_branches_passing"]
        + counts["fixed_sign_branches_failing"]
        or len(outcomes) != counts["fixed_sign_branches_total"]
    ):
        raise AssertionError("fixed sign accounting is incomplete")
    return feasible, {
        "applicable": True,
        "reason": "some_fixed_sign_branch_passes" if feasible
        else "all_fixed_sign_branches_fail",
        "basis_vertices": [instance.outside[local] for local in basis_locals],
        "basis_vertices_by_coordinate": {
            str(coordinate): instance.outside[local]
            for local, coordinate in coordinate_of.items()
        },
        "missing_coordinate": missing.bit_length() - 1,
        "fixed_points": [
            {
                "vertex": point.absolute,
                "basis_neighbourhood": (
                    base.parent.parent.kernel.vertices(point.neighbourhood)
                ),
                "allowed_seed_defects": (
                    base.parent.parent.kernel.vertices(point.allowed)
                ),
                "assigned_actual_support": None if point.assigned is None
                else base.parent.parent.kernel.vertices(point.assigned),
                "domain_failure": point.domain_failure,
            }
            for point in fixed
        ],
        "free_vertices_ignored": free_vertices,
        "counts": dict(counts),
        "branch_outcomes": outcomes,
    }


def strengthened_leaf_rule(
    adjacency: Sequence[int], instance, z0: int, bins: tuple[int, int],
    assignments: dict[int, int], counts: Counter[str], witnesses: list[dict],
) -> bool:
    """Frozen six-basis rule followed by all five-basis extensions."""

    LOCAL_COUNTS["complete_parent_leaves_reached"] += 1
    if not ORIGINAL_LEAF_RULE(
        adjacency, instance, z0, bins, assignments, counts, witnesses
    ):
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
            local for local in base.parent.parent.kernel.vertices(bin_mask)
            if assignments[local].bit_count() == 1
        )
        for basis_locals in itertools.combinations(singleton_locals, 5):
            leaf_has_basis = True
            feasible, detail = fixed_remainder_extension(
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


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    base.parent.parent.prior.activate_kernel()
    base.singleton_basis_leaf_passes = strengthened_leaf_rule
    base.solve_seed = instrumented_solve_seed


def evaluate_record(record: dict) -> dict:
    global LOCAL_COUNTS, LOCAL_EXTENSIONS, LOCAL_SEEDS
    LOCAL_COUNTS = Counter()
    LOCAL_EXTENSIONS = []
    LOCAL_SEEDS = []
    row = base.evaluate_record(record)
    if len(LOCAL_SEEDS) != row["seeds_checked"]:
        raise AssertionError("seed instrumentation coverage differs")
    row["fixed_remainder_counts"] = {
        key: LOCAL_COUNTS[key] for key in COUNTER_KEYS
    }
    row["fixed_remainder_seeds"] = list(LOCAL_SEEDS)
    row["fixed_remainder_extensions"] = list(LOCAL_EXTENSIONS)
    return row


def load_records() -> tuple[list[dict], list[int]]:
    if base.sha256(PARENT_REPORT) != EXPECTED_PARENT_REPORT_SHA256:
        raise ValueError("certified singleton-basis report boundary changed")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    indices = list(map(int, report["ordered_residue_indices"]))
    if base.stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("certified 249-graph residue changed")
    records, _, _ = base.load_input()
    by_index = {int(record["index"]): record for record in records}
    return [by_index[index] for index in indices], indices


def run_census(records: Sequence[dict], workers: int) -> list[dict]:
    if workers == 1:
        worker_initializer()
        return list(map(evaluate_record, records))
    with ProcessPoolExecutor(
        max_workers=workers, initializer=worker_initializer
    ) as pool:
        return list(pool.map(evaluate_record, records, chunksize=1))


def aggregate_counts(results: Sequence[dict]) -> dict[str, int]:
    total: Counter[str] = Counter()
    for row in results:
        total.update(row["fixed_remainder_counts"])
    rendered = {key: total[key] for key in COUNTER_KEYS}
    if (
        rendered["five_basis_extensions_total"]
        != rendered["five_basis_extensions_passing"]
        + rendered["five_basis_extensions_failing"]
        or rendered["fixed_sign_branches_total"]
        != rendered["fixed_sign_branches_passing"]
        + rendered["fixed_sign_branches_failing"]
    ):
        raise AssertionError("global extension/branch accounting differs")
    return rendered


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    records, indices = load_records()
    started = time.monotonic()
    results = run_census(records, args.workers)
    rejected = [int(row["index"]) for row in results if row["rejected"]]
    residue = [int(row["index"]) for row in results if not row["rejected"]]
    totals = aggregate_counts(results)
    if (
        rejected != list(EXPECTED_REJECTED_INDICES)
        or base.stable_hash(rejected) != EXPECTED_REJECTED_SHA256
        or base.stable_hash(residue) != EXPECTED_RESIDUE_SHA256
        or totals != EXPECTED_TOTALS
    ):
        raise AssertionError("fixed-remainder decisions changed")

    # Exercise every labeled K6 seed of the positive control under the exact
    # same patched theorem path, and retain its per-seed accounting.
    global LOCAL_COUNTS, LOCAL_EXTENSIONS, LOCAL_SEEDS
    LOCAL_COUNTS = Counter()
    LOCAL_EXTENSIONS = []
    LOCAL_SEEDS = []
    base.singleton_basis_leaf_passes = strengthened_leaf_rule
    base.solve_seed = instrumented_solve_seed
    positive = base.positive_control()
    if positive != {"passed": True, "K6_seeds": 32} or len(LOCAL_SEEDS) != 32:
        raise AssertionError("fixed-remainder positive control changed")
    positive_seed_rows = list(LOCAL_SEEDS)

    graph_rows = [{
        "index": int(row["index"]),
        "rejected": bool(row["rejected"]),
        "seeds_checked": int(row["seeds_checked"]),
        "first_impossible_seed_mask": row["first_impossible_seed_mask"],
        "counts": row["fixed_remainder_counts"],
        "seeds": row["fixed_remainder_seeds"],
        "extensions": row["fixed_remainder_extensions"],
    } for row in results]
    if (
        base.stable_hash(graph_rows) != EXPECTED_GRAPH_ROWS_SHA256
        or base.stable_hash(positive_seed_rows)
        != EXPECTED_POSITIVE_SEED_ROWS_SHA256
    ):
        raise AssertionError("fixed-remainder exhaustive transcript changed")
    output = {
        "schema": "d6-k6-five-singleton-fixed-remainder-pilot-v1",
        "status": "COMPLETE_DISCOVERY_ONLY",
        "input_graphs": len(records),
        "ordered_input_indices_sha256": base.stable_hash(indices),
        "marginal_rejections": len(rejected),
        "marginal_indices": rejected,
        "marginal_indices_sha256": base.stable_hash(rejected),
        "residue": len(residue),
        "ordered_residue_indices": residue,
        "ordered_residue_indices_sha256": base.stable_hash(residue),
        "totals": totals,
        "graph_rows": graph_rows,
        "graph_rows_sha256": base.stable_hash(graph_rows),
        "positive_18_control": positive,
        "positive_seed_rows": positive_seed_rows,
        "positive_seed_rows_sha256": base.stable_hash(positive_seed_rows),
        "runtime": {
            "workers": args.workers,
            "wall_seconds": time.monotonic() - started,
        },
        "sources": {
            Path(__file__).name: base.sha256(Path(__file__)),
            Path(base.__file__).name: base.sha256(Path(base.__file__)),
            PARENT_REPORT.name: base.sha256(PARENT_REPORT),
        },
        "semantics": {
            "candidate_nonedges": "no distance equation",
            "allowed_supports": "upper bounds; optional zeros retained",
            "free_missing_coordinate_points": "ignored relaxation",
            "arithmetic": "integer masks and exhaustive binary roots",
            "root": "t=(1+epsilon*sqrt(7))/3",
        },
    }
    base.atomic_json(args.output, output)
    print(json.dumps({
        "marginal_rejections": len(rejected),
        "residue": len(residue),
        "totals": totals,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
