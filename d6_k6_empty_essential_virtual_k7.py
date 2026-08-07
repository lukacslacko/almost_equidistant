#!/usr/bin/env python3
"""Exact virtual-K7 closure of the K6 empty-essential branches.

The independently checked empty-essential profile supplies every K6 seed for
which the current Hall relaxation has no all-nonempty branch, together with
the complete relaxed list of possible one- or two-empty states.  For every
vertex appearing in one of those states, this program branches on its K6
defect vector being zero, adds the six branch-forced unit edges to the K6,
and subjects the resulting actual K7 to the generic exact K7 cover/rank
quantifier.

An empty branch is eliminated only if every eligible zero-factor cover is
eliminated: covers of size 4--7 by the proved direct dimension caps, and
covers of size at most 3 by an exact support/rank obstruction.  Cover sizes
above 7 are impossible for the orthogonal zero-factor family.  No candidate
nonedge is made non-unit.  A graph is rejected only when one of its K6 seeds
has no surviving all-nonempty or empty branch.
"""

from __future__ import annotations

import argparse
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

import d6_k6_empty_essential_profile as profile_parent
import d6_k7_rank_reference as k7


ROOT = Path(__file__).resolve().parent
PROFILE_REPORT = ROOT / "d6_k6_empty_essential_profile.json"
PROFILE_VERIFICATION = ROOT / "d6_k6_empty_essential_profile_verification.json"
EXPECTED_PROFILE_REPORT_SHA256 = (
    "205c0877982c8de8ad6e5a685e33eff9ee05236c9907e60d488a593792d90781"
)
EXPECTED_PROFILE_VERIFICATION_SHA256 = (
    "c8f73e8a9bed6d7a102ab168eb3842312a916b23dbd577cceb6d089c8832f719"
)
EXPECTED_PROFILE_SOURCE_SHA256 = (
    "dab7603b11946635a29c3c37daf44a895db342b2cdddb17c3699f87c182f507a"
)
EXPECTED_K7_SOURCE_SHA256 = (
    "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
)
EXPECTED_INPUT = 805
EXPECTED_INPUT_SHA256 = (
    "53cad35eedf4f087d1e1dadbffcde54cbbbf04dbd5ee6ea6e10506329c958818"
)
EXPECTED_TARGET_GRAPHS = 49
EXPECTED_TARGET_SEEDS = 111
EXPECTED_PROMOTIONS = 291
REPORT_SCHEMA = "d6-k6-empty-essential-virtual-k7-v1"
DEFAULT_WORKERS = 11


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
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


def load_inputs() -> tuple[list[tuple[dict, dict]], list[int], dict[str, str]]:
    """Bind the full 805 boundary and select its 49 target graphs."""

    observed = {
        PROFILE_REPORT.name: sha256(PROFILE_REPORT),
        PROFILE_VERIFICATION.name: sha256(PROFILE_VERIFICATION),
        "d6_k6_empty_essential_profile.py": sha256(
            ROOT / "d6_k6_empty_essential_profile.py"
        ),
        "d6_k7_rank_reference.py": sha256(ROOT / "d6_k7_rank_reference.py"),
    }
    expected = {
        PROFILE_REPORT.name: EXPECTED_PROFILE_REPORT_SHA256,
        PROFILE_VERIFICATION.name: EXPECTED_PROFILE_VERIFICATION_SHA256,
        "d6_k6_empty_essential_profile.py": EXPECTED_PROFILE_SOURCE_SHA256,
        "d6_k7_rank_reference.py": EXPECTED_K7_SOURCE_SHA256,
    }
    if observed != expected:
        raise ValueError(f"virtual-K7 dependency boundary changed: {observed}")

    records, indices, _ = profile_parent.load_residue()
    report = json.loads(PROFILE_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PROFILE_VERIFICATION.read_text(encoding="utf-8"))
    aggregate = report.get("aggregate", {})
    if (
        report.get("schema") != "d6-k6-empty-essential-profile-v1"
        or report.get("status") != "COMPLETE"
        or report.get("ordered_input_indices") != indices
        or len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or aggregate.get("totals", {}).get("graphs_with_empty_essential_seed")
        != EXPECTED_TARGET_GRAPHS
        or aggregate.get("seed_classification_counts", {}).get("EMPTY_ESSENTIAL")
        != EXPECTED_TARGET_SEEDS
        or len(aggregate.get("virtual_K7_promotion_keys", []))
        != EXPECTED_PROMOTIONS
    ):
        raise ValueError("empty-essential profile contents changed")
    if (
        verification.get("schema")
        != "d6-k6-empty-essential-profile-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("report_sha256") != EXPECTED_PROFILE_REPORT_SHA256
        or verification.get("graphs_recomputed") != EXPECTED_INPUT
        or verification.get("K6_seeds_recomputed") != 25_354
        or verification.get("empty_essential_seeds") != EXPECTED_TARGET_SEEDS
        or verification.get("graphs_with_empty_essential_seed")
        != EXPECTED_TARGET_GRAPHS
        or verification.get("virtual_K7_promotion_keys") != EXPECTED_PROMOTIONS
    ):
        raise ValueError("empty-essential independent verification changed")

    by_index = {int(record["index"]): record for record in records}
    targets = []
    promotion_keys = []
    for graph_row in report["graph_results"]:
        essential = [
            row
            for row in graph_row["seed_results"]
            if row["classification"] == "EMPTY_ESSENTIAL"
        ]
        if not essential:
            continue
        index = int(graph_row["index"])
        targets.append((by_index[index], {"index": index, "seeds": essential}))
        for seed in essential:
            promotion_keys.extend(
                [index, int(seed["seed_mask"]), int(apex)]
                for apex in seed["promotion_apices"]
            )
    if (
        len(targets) != EXPECTED_TARGET_GRAPHS
        or sum(len(row["seeds"]) for _, row in targets) != EXPECTED_TARGET_SEEDS
        or len(promotion_keys) != EXPECTED_PROMOTIONS
        or promotion_keys != aggregate["virtual_K7_promotion_keys"]
    ):
        raise ValueError("target extraction differs from the checked profile")
    return targets, indices, observed


def add_required_edge(adjacency: list[int], first: int, second: int) -> None:
    adjacency[first] |= 1 << second
    adjacency[second] |= 1 << first


def augment_with_apex(
    adjacency: Sequence[int], seed: Sequence[int], apex: int
) -> tuple[tuple[int, ...], list[list[int]], int]:
    if len(seed) != 6 or len(set(seed)) != 6 or apex in seed:
        raise ValueError("virtual promotion requires a K6 and outside apex")
    augmented = list(map(int, adjacency))
    added = []
    for vertex in seed:
        if not (augmented[apex] & (1 << vertex)):
            add_required_edge(augmented, apex, vertex)
            added.append([min(apex, vertex), max(apex, vertex)])
    k7.validate_graph(augmented, require_alpha_two=True)
    seed_mask = sum(1 << vertex for vertex in (*seed, apex))
    intended = list(k7.bits(seed_mask))
    if len(intended) != 7 or any(
        not (augmented[first] & (1 << second))
        for i, first in enumerate(intended)
        for second in intended[:i]
    ):
        raise AssertionError("branch augmentation did not create intended K7")
    return tuple(augmented), added, seed_mask


def failure_names(analysis: k7.CoverAnalysis) -> list[str]:
    tests = (
        ("zero_factor_support", analysis.support_failed),
        ("subspace_K", analysis.subspace_k_failed),
        ("component_B", analysis.component_b_failed),
        ("PD_clique", analysis.pd_clique_failed),
        ("perpendicular_degree", analysis.perpendicular_degree_failed),
        ("basis_kernel", analysis.basis_kernel_failed),
        ("saturating_mask", analysis.saturating_mask_failed),
    )
    return [name for name, failed in tests if failed]


def compact_cover_certificate(
    zmask: int, outside: Sequence[int], analysis: k7.CoverAnalysis
) -> dict:
    return {
        "zmask": zmask,
        "Z_vertices": [outside[local] for local in k7.bits(zmask)],
        "failures": failure_names(analysis),
        "support_failed": analysis.support_failed,
        "subspace_K_failed": analysis.subspace_k_failed,
        "component_B_failed": analysis.component_b_failed,
        "PD_clique_failed": analysis.pd_clique_failed,
        "perpendicular_degree_failed": analysis.perpendicular_degree_failed,
        "basis_kernel_failed": analysis.basis_kernel_failed,
        "saturating_mask_failed": analysis.saturating_mask_failed,
        "N_size": analysis.n_size,
        "term_rank": analysis.term_rank,
        "total_term_rank": analysis.total_term_rank,
        "K_rank_upper": analysis.k_rank_upper,
        "K_clique_number": analysis.k_clique_number,
        "F_maximum_degree": analysis.f_maximum_degree,
        "K_zero_forcing_number": analysis.k_zero_forcing,
        "B_zero_forcing_number": analysis.b_zero_forcing,
        "B_component_nullity_cap": analysis.component_nullity_cap,
        "B_component_tier_failures": analysis.component_tier_failures,
        "basis_kernel_clique_failures": analysis.basis_kernel_clique_failures,
        "basis_kernel_failure_reasons": analysis.basis_kernel_failure_reasons,
        "saturating_mask_clique_failures": (
            analysis.saturating_mask_clique_failures
        ),
        "saturating_mask_failure_reasons": (
            analysis.saturating_mask_failure_reasons
        ),
    }


def screen_k7_seed(
    adjacency: Sequence[int],
    seed_mask: int,
    support_solver: k7.SupportSolver,
    zero_forcing: k7.ZeroForcingSolver,
    clique_solver: k7.CliqueStructureSolver,
) -> dict:
    """Exhaust every eligible cover in the complete K7 necessary system."""

    seed, outside, defects, ladj, eligible = k7.seed_instance(adjacency, seed_mask)
    covers = k7.eligible_covers(ladj, eligible, cap=7)
    total_term_rank = k7.matching_size(defects)
    size_histogram = Counter(str(cover.bit_count()) for cover in covers)
    direct = [cover for cover in covers if cover.bit_count() >= 4]
    small = [cover for cover in covers if cover.bit_count() <= 3]
    for cover in direct:
        if k7.direct_cover_cap_failure(cover.bit_count()) is None:
            raise AssertionError("large cover lacks its exact direct cap")
    small_rows = []
    passing_small = []
    failures: Counter[str] = Counter()
    for cover in small:
        analysis = k7.analyze_cover(
            adjacency,
            outside,
            defects,
            cover,
            support_solver,
            zero_forcing,
            total_term_rank,
            clique_solver,
        )
        row = compact_cover_certificate(cover, outside, analysis)
        small_rows.append(row)
        failures.update(row["failures"])
        if not analysis.enhanced_joint_failed:
            passing_small.append(cover)
        elif not row["failures"]:
            raise AssertionError("failed small cover has no exact failure name")
    eliminated = not passing_small
    return {
        "seed": seed,
        "outside": outside,
        "defect_masks": list(defects),
        "eligible_vertex_mask": eligible,
        "eligible_covers": len(covers),
        "eligible_cover_masks_sha256": stable_hash(covers),
        "cover_size_histogram": dict(size_histogram),
        "direct_cap_covers": len(direct),
        "small_covers": len(small),
        "small_cover_failure_counts": dict(failures),
        "small_cover_certificates": small_rows,
        "passing_small_covers": passing_small,
        "status": "ELIMINATED" if eliminated else "SURVIVOR",
        "elimination_quantifier": {
            "covers_above_7": "impossible_by_orthogonal_dimension",
            "covers_4_through_7": "eliminated_by_direct_cover_caps",
            "covers_0_through_3": "listed_and_exactly_screened",
            "all_eligible_covers_eliminated": eliminated,
        },
    }


def analyze_promotion(
    adjacency: Sequence[int],
    seed: Sequence[int],
    seed_mask: int,
    apex: int,
    support_solver: k7.SupportSolver,
    zero_forcing: k7.ZeroForcingSolver,
    clique_solver: k7.CliqueStructureSolver,
) -> dict:
    augmented, added_edges, promoted_mask = augment_with_apex(
        adjacency, seed, apex
    )
    if sum(1 << vertex for vertex in seed) != seed_mask:
        raise ValueError("profile seed labels do not equal its seed mask")
    result = screen_k7_seed(
        augmented, promoted_mask, support_solver, zero_forcing, clique_solver
    )
    return {
        "apex": apex,
        "added_edges": added_edges,
        "promoted_seed_mask": promoted_mask,
        **result,
    }


def evaluate_target(payload: tuple[dict, dict]) -> dict:
    record, target = payload
    adjacency = tuple(map(int, record["adjacency"]))
    support_solver = k7.SupportSolver()
    zero_forcing = k7.ZeroForcingSolver()
    clique_solver = k7.CliqueStructureSolver()
    seed_rows = []
    for profile_seed in target["seeds"]:
        seed = [int(vertex) for vertex in profile_seed["seed"]]
        seed_mask = int(profile_seed["seed_mask"])
        promotions = [
            analyze_promotion(
                adjacency,
                seed,
                seed_mask,
                int(apex),
                support_solver,
                zero_forcing,
                clique_solver,
            )
            for apex in profile_seed["promotion_apices"]
        ]
        killed = {
            row["apex"] for row in promotions if row["status"] == "ELIMINATED"
        }
        surviving_states = [
            state
            for state in profile_seed["possible_empty_states"]
            if not any(int(apex) in killed for apex in state)
        ]
        seed_rows.append(
            {
                "seed_mask": seed_mask,
                "seed": seed,
                "profile_possible_empty_states": profile_seed[
                    "possible_empty_states"
                ],
                "promotions": promotions,
                "surviving_empty_states": surviving_states,
                "status": "IMPOSSIBLE" if not surviving_states else "UNRESOLVED",
            }
        )
    return {
        "index": int(record["index"]),
        "target_seeds": len(seed_rows),
        "impossible_seed_masks": [
            row["seed_mask"] for row in seed_rows if row["status"] == "IMPOSSIBLE"
        ],
        "rejected": any(row["status"] == "IMPOSSIBLE" for row in seed_rows),
        "seed_results": seed_rows,
    }


def complete_graph(size: int) -> tuple[int, ...]:
    full = (1 << size) - 1
    return tuple(full & ~(1 << vertex) for vertex in range(size))


def controls() -> dict:
    support = k7.SupportSolver()
    forcing = k7.ZeroForcingSolver()
    clique = k7.CliqueStructureSolver()
    positive = screen_k7_seed(
        complete_graph(7), (1 << 7) - 1, support, forcing, clique
    )
    negative = screen_k7_seed(
        complete_graph(8), (1 << 7) - 1, support, forcing, clique
    )
    if positive["status"] != "SURVIVOR" or positive["passing_small_covers"] != [0]:
        raise AssertionError("realizable regular K7 positive control was eliminated")
    if negative["status"] != "ELIMINATED":
        raise AssertionError("K8 negative control unexpectedly survived in R6")
    return {
        "positive_regular_K7": {
            "passed": True,
            "eligible_covers": positive["eligible_covers"],
            "passing_small_covers": positive["passing_small_covers"],
        },
        "negative_K8": {
            "passed": True,
            "eligible_covers": negative["eligible_covers"],
            "small_cover_failure_counts": negative[
                "small_cover_failure_counts"
            ],
        },
    }


def aggregate(
    results: Sequence[dict], ordered_input: Sequence[int]
) -> dict:
    status_counts: Counter[str] = Counter()
    cover_sizes: Counter[str] = Counter()
    failure_counts: Counter[str] = Counter()
    totals: Counter[str] = Counter()
    impossible_keys = []
    promotion_keys = []
    no_cover_promotions = 0
    for graph in results:
        totals.update(
            target_graphs=1,
            target_seeds=graph["target_seeds"],
            impossible_seeds=len(graph["impossible_seed_masks"]),
            rejected_graphs=int(graph["rejected"]),
        )
        for seed in graph["seed_results"]:
            if seed["status"] == "IMPOSSIBLE":
                impossible_keys.append([graph["index"], seed["seed_mask"]])
            for promotion in seed["promotions"]:
                totals.update(
                    promotions=1,
                    eligible_covers=promotion["eligible_covers"],
                    direct_cap_covers=promotion["direct_cap_covers"],
                    small_covers=promotion["small_covers"],
                )
                promotion_keys.append(
                    [graph["index"], seed["seed_mask"], promotion["apex"]]
                )
                status_counts[promotion["status"]] += 1
                cover_sizes.update(promotion["cover_size_histogram"])
                failure_counts.update(promotion["small_cover_failure_counts"])
                no_cover_promotions += promotion["eligible_covers"] == 0
    rejected = [graph["index"] for graph in results if graph["rejected"]]
    rejected_set = set(rejected)
    residue = [index for index in ordered_input if index not in rejected_set]
    return {
        "totals": dict(totals),
        "promotion_status_counts": dict(status_counts),
        "promotions_with_no_eligible_cover": no_cover_promotions,
        "eligible_cover_size_histogram": dict(cover_sizes),
        "small_cover_failure_counts": dict(failure_counts),
        "promotion_keys": promotion_keys,
        "promotion_keys_sha256": stable_hash(promotion_keys),
        "impossible_seed_keys": impossible_keys,
        "impossible_seed_keys_sha256": stable_hash(impossible_keys),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices": residue,
        "ordered_residue_indices_sha256": stable_hash(residue),
    }


def run(output: Path, workers: int) -> dict:
    targets, ordered_input, dependencies = load_inputs()
    started = time.perf_counter()
    if workers == 1:
        results = [evaluate_target(payload) for payload in targets]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_target, targets, chunksize=1))
    if [row["index"] for row in results] != [row[1]["index"] for row in targets]:
        raise AssertionError("process pool changed target order")
    summary = aggregate(results, ordered_input)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact virtual-K7 elimination of every relaxed empty state at "
            "the empty-essential K6 seeds."
        ),
        "input_graphs": len(ordered_input),
        "input_indices_sha256": stable_hash(ordered_input),
        "target_graphs": len(results),
        "aggregate": summary,
        "graph_results": results,
        "controls": controls(),
        "sources": {
            **dependencies,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "semantics": {
            "branch_implication": (
                "u_Q(x)=0 forces all six x-Q pairs unit and makes Q union {x} "
                "an actual regular unit K7"
            ),
            "candidate_nonedges": "unconstrained and may also be unit",
            "added_edges": "only branch-forced apex-to-K6 unit pairs",
            "eligible_cover_quantifier": (
                "all eligible covers of sizes 0..7; sizes above 7 impossible"
            ),
            "empty_state_quantifier": (
                "complete necessary state list from the independently checked "
                "empty-essential profile"
            ),
            "double_apex_test": (
                "not needed when every constituent singleton promotion is eliminated"
            ),
            "survivor": "filter non-rejection only",
            "arithmetic": "exact integer, bit-mask, and rational logic only",
        },
        "runtime": {
            "command": list(sys.argv),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_empty_essential_virtual_k7.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > DEFAULT_WORKERS:
        parser.error(f"workers must lie in 1..{DEFAULT_WORKERS}")
    report = run(args.output.resolve(), args.workers)
    aggregate_row = report["aggregate"]
    print(
        json.dumps(
            {
                "status": report["status"],
                "promotions": aggregate_row["totals"]["promotions"],
                "impossible_seeds": aggregate_row["totals"]["impossible_seeds"],
                "rejected_graphs": aggregate_row["totals"]["rejected_graphs"],
                "residue_graphs": len(aggregate_row["ordered_residue_indices"]),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
