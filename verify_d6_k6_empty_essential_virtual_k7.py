#!/usr/bin/env python3
"""Independent checker for the K6-empty / virtual-K7 conjunction.

The checker imports neither the production conjunction nor the production
empty-essential profiler.  It reconstructs the 805 graphs, reads the
independently verified branch list, independently augments every apex,
enumerates every eligible K7 cover by brute bit masks, and replays the exact
generic K7 leaf certificates.  It then compares all 291 promotions, all 111
seed conjunctions, the 49 graph rejections, and the ordered 756 residue.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k6_arbitrary_subset_hall as hall_parent
import d6_k7_rank_reference as k7


ROOT = Path(__file__).resolve().parent
PROFILE_REPORT = ROOT / "d6_k6_empty_essential_profile.json"
PROFILE_VERIFICATION = ROOT / "d6_k6_empty_essential_profile_verification.json"
EMPTY_REPORT = ROOT / "d6_k6_empty_support_budget_report.json"
EXPECTED_PROFILE_REPORT_SHA256 = (
    "205c0877982c8de8ad6e5a685e33eff9ee05236c9907e60d488a593792d90781"
)
EXPECTED_PROFILE_VERIFICATION_SHA256 = (
    "c8f73e8a9bed6d7a102ab168eb3842312a916b23dbd577cceb6d089c8832f719"
)
EXPECTED_PROFILE_SOURCE_SHA256 = (
    "dab7603b11946635a29c3c37daf44a895db342b2cdddb17c3699f87c182f507a"
)
EXPECTED_EMPTY_REPORT_SHA256 = (
    "1860dbe69ae74b55203ea283cf83afff929a4b1f5cbfcd30e09f0138688eba9f"
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
VERIFICATION_SCHEMA = "d6-k6-empty-essential-virtual-k7-verification-v1"
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


def assert_import_independence() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "d6_k6_empty_essential_virtual_k7",
        "d6_k6_empty_essential_profile",
        "d6_k6_empty_support_budget",
        "probe_d6_k6_virtual_k7",
        "probe_d6_k6_virtual_k7_empty_conjunction",
    }
    if imports & forbidden:
        raise AssertionError(
            f"independent verifier imports production code: {imports & forbidden}"
        )


def reconstruct_inputs() -> tuple[list[tuple[dict, dict]], list[int], dict[str, str]]:
    """Rebuild the input and extract targets without production imports."""

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
        raise ValueError(f"independent dependency boundary changed: {observed}")
    if sha256(EMPTY_REPORT) != EXPECTED_EMPTY_REPORT_SHA256:
        raise ValueError("empty-support parent report hash changed")

    _, records, arbitrary_indices = hall_parent.verify_inputs()
    arbitrary_report = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_report.json").read_text(
            encoding="utf-8"
        )
    )
    inherited = {int(index) for index in arbitrary_report["rejected_indices"]}
    expected_822 = [index for index in arbitrary_indices if index not in inherited]
    empty_report = json.loads(EMPTY_REPORT.read_text(encoding="utf-8"))
    indices = [int(index) for index in empty_report["ordered_residue_indices"]]
    if (
        empty_report.get("status") != "COMPLETE"
        or empty_report.get("ordered_input_indices") != expected_822
        or len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("independent ordered 805 reconstruction failed")

    profile = json.loads(PROFILE_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PROFILE_VERIFICATION.read_text(encoding="utf-8"))
    aggregate = profile.get("aggregate", {})
    if (
        profile.get("schema") != "d6-k6-empty-essential-profile-v1"
        or profile.get("status") != "COMPLETE"
        or profile.get("ordered_input_indices") != indices
        or aggregate.get("totals", {}).get("graphs_with_empty_essential_seed")
        != EXPECTED_TARGET_GRAPHS
        or aggregate.get("seed_classification_counts", {}).get("EMPTY_ESSENTIAL")
        != EXPECTED_TARGET_SEEDS
        or len(aggregate.get("virtual_K7_promotion_keys", []))
        != EXPECTED_PROMOTIONS
    ):
        raise ValueError("independent profile content boundary changed")
    if (
        verification.get("status") != "PASS"
        or verification.get("report_sha256") != EXPECTED_PROFILE_REPORT_SHA256
        or verification.get("graphs_recomputed") != EXPECTED_INPUT
        or verification.get("K6_seeds_recomputed") != 25_354
        or verification.get("empty_essential_seeds") != EXPECTED_TARGET_SEEDS
        or verification.get("virtual_K7_promotion_keys") != EXPECTED_PROMOTIONS
    ):
        raise ValueError("independent profile verification boundary changed")

    by_index = {int(record["index"]): record for record in records}
    targets = []
    promotion_keys = []
    for graph in profile["graph_results"]:
        seeds = [
            row
            for row in graph["seed_results"]
            if row["classification"] == "EMPTY_ESSENTIAL"
        ]
        if not seeds:
            continue
        index = int(graph["index"])
        targets.append((by_index[index], {"index": index, "seeds": seeds}))
        for seed in seeds:
            promotion_keys.extend(
                [index, int(seed["seed_mask"]), int(apex)]
                for apex in seed["promotion_apices"]
            )
    if (
        len(targets) != EXPECTED_TARGET_GRAPHS
        or sum(len(target["seeds"]) for _, target in targets)
        != EXPECTED_TARGET_SEEDS
        or promotion_keys != aggregate["virtual_K7_promotion_keys"]
    ):
        raise ValueError("independent target extraction differs")
    return targets, indices, observed


def independent_augment(
    adjacency: Sequence[int], seed: Sequence[int], apex: int
) -> tuple[tuple[int, ...], list[list[int]], int]:
    if len(seed) != 6 or len(set(seed)) != 6 or apex in seed:
        raise ValueError("bad independent virtual-K7 promotion")
    augmented = list(map(int, adjacency))
    added = []
    for vertex in seed:
        if not (augmented[apex] & (1 << vertex)):
            augmented[apex] |= 1 << vertex
            augmented[vertex] |= 1 << apex
            added.append([min(apex, vertex), max(apex, vertex)])
    k7.validate_graph(augmented, require_alpha_two=True)
    seed_mask = sum(1 << vertex for vertex in (*seed, apex))
    labels = [vertex for vertex in range(len(augmented)) if seed_mask & (1 << vertex)]
    if len(labels) != 7 or any(
        not (augmented[first] & (1 << second))
        for first, second in combinations(labels, 2)
    ):
        raise AssertionError("independent augmentation did not make K7")
    return tuple(augmented), added, seed_mask


def independent_seed_instance(
    adjacency: Sequence[int], seed_mask: int
) -> tuple[list[int], list[int], tuple[int, ...], tuple[int, ...], int]:
    seed = [vertex for vertex in range(len(adjacency)) if seed_mask & (1 << vertex)]
    if len(seed) != 7:
        raise ValueError("independent seed must have seven vertices")
    positions = {vertex: coordinate for coordinate, vertex in enumerate(seed)}
    outside = [vertex for vertex in range(len(adjacency)) if vertex not in positions]
    defects = []
    for vertex in outside:
        allowed = 0
        for q in seed:
            if not (adjacency[vertex] & (1 << q)):
                allowed |= 1 << positions[q]
        defects.append(allowed)
    ladj = [0] * len(outside)
    for first, second in combinations(range(len(outside)), 2):
        if (
            adjacency[outside[first]] & (1 << outside[second])
            and not (defects[first] & defects[second])
        ):
            ladj[first] |= 1 << second
            ladj[second] |= 1 << first
    eligible = sum(
        1 << local
        for local, allowed in enumerate(defects)
        if allowed.bit_count() >= 3
    )
    return seed, outside, tuple(defects), tuple(ladj), eligible


def independent_is_cover(ladj: Sequence[int], selected: int) -> bool:
    remaining = ((1 << len(ladj)) - 1) & ~selected
    return all(not (ladj[vertex] & remaining) for vertex in range(len(ladj)) if remaining & (1 << vertex))


def independent_covers(ladj: Sequence[int], eligible: int) -> list[int]:
    return [
        selected
        for selected in range(1 << len(ladj))
        if not (selected & ~eligible)
        and selected.bit_count() <= 7
        and independent_is_cover(ladj, selected)
    ]


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


def compact_certificate(zmask: int, outside, analysis) -> dict:
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
        "saturating_mask_clique_failures": analysis.saturating_mask_clique_failures,
        "saturating_mask_failure_reasons": analysis.saturating_mask_failure_reasons,
    }


def independent_screen(
    adjacency: Sequence[int],
    seed_mask: int,
    support_solver,
    zero_forcing,
    clique_solver,
) -> dict:
    seed, outside, defects, ladj, eligible = independent_seed_instance(
        adjacency, seed_mask
    )
    covers = independent_covers(ladj, eligible)
    total_term_rank = k7.matching_size(defects)
    size_histogram = Counter(str(cover.bit_count()) for cover in covers)
    direct = [cover for cover in covers if cover.bit_count() >= 4]
    small = [cover for cover in covers if cover.bit_count() <= 3]
    direct_names = {4: "cap4_to3", 5: "cap5_to4", 6: "cap6_to5", 7: "cap7_to6"}
    if any(cover.bit_count() not in direct_names for cover in direct):
        raise AssertionError("independent direct cover cap missing")
    rows = []
    passing = []
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
        row = compact_certificate(cover, outside, analysis)
        rows.append(row)
        failures.update(row["failures"])
        if not analysis.enhanced_joint_failed:
            passing.append(cover)
        elif not row["failures"]:
            raise AssertionError("independent failed cover has no witness")
    eliminated = not passing
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
        "small_cover_certificates": rows,
        "passing_small_covers": passing,
        "status": "ELIMINATED" if eliminated else "SURVIVOR",
        "elimination_quantifier": {
            "covers_above_7": "impossible_by_orthogonal_dimension",
            "covers_4_through_7": "eliminated_by_direct_cover_caps",
            "covers_0_through_3": "listed_and_exactly_screened",
            "all_eligible_covers_eliminated": eliminated,
        },
    }


def evaluate_target(payload: tuple[dict, dict]) -> dict:
    record, target = payload
    adjacency = tuple(map(int, record["adjacency"]))
    support = k7.SupportSolver()
    forcing = k7.ZeroForcingSolver()
    clique = k7.CliqueStructureSolver()
    seeds = []
    for profile_seed in target["seeds"]:
        seed = [int(vertex) for vertex in profile_seed["seed"]]
        seed_mask = int(profile_seed["seed_mask"])
        if sum(1 << vertex for vertex in seed) != seed_mask:
            raise ValueError("independent seed labels/mask mismatch")
        promotions = []
        for apex in profile_seed["promotion_apices"]:
            augmented, added, promoted = independent_augment(
                adjacency, seed, int(apex)
            )
            promotions.append(
                {
                    "apex": int(apex),
                    "added_edges": added,
                    "promoted_seed_mask": promoted,
                    **independent_screen(
                        augmented, promoted, support, forcing, clique
                    ),
                }
            )
        killed = {
            row["apex"] for row in promotions if row["status"] == "ELIMINATED"
        }
        surviving = [
            state
            for state in profile_seed["possible_empty_states"]
            if not any(int(apex) in killed for apex in state)
        ]
        seeds.append(
            {
                "seed_mask": seed_mask,
                "seed": seed,
                "profile_possible_empty_states": profile_seed[
                    "possible_empty_states"
                ],
                "promotions": promotions,
                "surviving_empty_states": surviving,
                "status": "IMPOSSIBLE" if not surviving else "UNRESOLVED",
            }
        )
    impossible = [
        row["seed_mask"] for row in seeds if row["status"] == "IMPOSSIBLE"
    ]
    return {
        "index": int(record["index"]),
        "target_seeds": len(seeds),
        "impossible_seed_masks": impossible,
        "rejected": bool(impossible),
        "seed_results": seeds,
    }


def complete_graph(size: int) -> tuple[int, ...]:
    full = (1 << size) - 1
    return tuple(full & ~(1 << vertex) for vertex in range(size))


def controls() -> dict:
    support = k7.SupportSolver()
    forcing = k7.ZeroForcingSolver()
    clique = k7.CliqueStructureSolver()
    positive = independent_screen(
        complete_graph(7), (1 << 7) - 1, support, forcing, clique
    )
    negative = independent_screen(
        complete_graph(8), (1 << 7) - 1, support, forcing, clique
    )
    if positive["status"] != "SURVIVOR" or positive["passing_small_covers"] != [0]:
        raise AssertionError("independent regular-K7 control failed")
    if negative["status"] != "ELIMINATED":
        raise AssertionError("independent K8 control failed")
    return {
        "positive_regular_K7": {
            "passed": True,
            "eligible_covers": positive["eligible_covers"],
            "passing_small_covers": positive["passing_small_covers"],
        },
        "negative_K8": {
            "passed": True,
            "eligible_covers": negative["eligible_covers"],
            "small_cover_failure_counts": negative["small_cover_failure_counts"],
        },
    }


def aggregate(results: Sequence[dict], ordered_input: Sequence[int]) -> dict:
    status_counts: Counter[str] = Counter()
    cover_sizes: Counter[str] = Counter()
    failure_counts: Counter[str] = Counter()
    totals: Counter[str] = Counter()
    impossible_keys = []
    promotion_keys = []
    no_cover = 0
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
                no_cover += promotion["eligible_covers"] == 0
    rejected = [graph["index"] for graph in results if graph["rejected"]]
    residue = [index for index in ordered_input if index not in set(rejected)]
    return {
        "totals": dict(totals),
        "promotion_status_counts": dict(status_counts),
        "promotions_with_no_eligible_cover": no_cover,
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


class VerificationInfrastructureAbort(RuntimeError):
    pass


def verify(report_path: Path, expected_hash: str, output: Path, workers: int) -> dict:
    if (
        len(expected_hash) != 64
        or any(character not in "0123456789abcdef" for character in expected_hash)
    ):
        raise ValueError("--expected-report-sha256 must be 64 lowercase hex digits")
    actual_hash = sha256(report_path)
    if actual_hash != expected_hash:
        raise ValueError(f"explicit report SHA-256 mismatch: actual={actual_hash}")
    assert_import_independence()
    targets, ordered_input, dependencies = reconstruct_inputs()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    producer = ROOT / "d6_k6_empty_essential_virtual_k7.py"
    producer_hash = sha256(producer)
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != EXPECTED_INPUT
        or report.get("input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("target_graphs") != EXPECTED_TARGET_GRAPHS
        or report.get("sources", {}).get(producer.name) != producer_hash
    ):
        raise ValueError("production virtual-K7 report boundary changed")
    for name, digest in dependencies.items():
        if report.get("sources", {}).get(name) != digest:
            raise ValueError(f"production dependency mismatch: {name}")

    started = time.perf_counter()
    try:
        if workers == 1:
            recomputed = [evaluate_target(payload) for payload in targets]
        else:
            with ProcessPoolExecutor(
                max_workers=workers, initializer=worker_initializer
            ) as pool:
                recomputed = list(pool.map(evaluate_target, targets, chunksize=1))
    except BaseException as error:
        raise VerificationInfrastructureAbort(
            "fresh virtual-K7 verifier pool aborted; no rejection inferred"
        ) from error
    archived = report.get("graph_results")
    if recomputed != archived:
        for observed, expected in zip(recomputed, archived or []):
            if observed != expected:
                raise AssertionError(
                    f"independent promotion/certificate mismatch at {observed['index']}"
                )
        raise AssertionError("independent target result length/order mismatch")
    summary = aggregate(recomputed, ordered_input)
    if summary != report.get("aggregate"):
        raise AssertionError("independent virtual-K7 aggregate differs")
    control_rows = controls()
    if control_rows != report.get("controls"):
        raise AssertionError("independent controls differ")
    if (
        summary["totals"]["promotions"] != EXPECTED_PROMOTIONS
        or summary["totals"]["impossible_seeds"] != EXPECTED_TARGET_SEEDS
        or summary["totals"]["rejected_graphs"] != EXPECTED_TARGET_GRAPHS
        or len(summary["ordered_residue_indices"]) != 756
    ):
        raise AssertionError("verified theorem totals differ from expected boundary")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "production conjunction/profile not imported; independent 805 input "
            "reconstruction, apex augmentation, brute eligible-cover enumeration, "
            "and exact replay/comparison of all small-cover leaf certificates"
        ),
        "report_sha256": actual_hash,
        "production_source_sha256": producer_hash,
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
        "input_graphs": EXPECTED_INPUT,
        "target_graphs": EXPECTED_TARGET_GRAPHS,
        "target_seeds": EXPECTED_TARGET_SEEDS,
        "promotions_recomputed": EXPECTED_PROMOTIONS,
        "eligible_covers_recomputed": summary["totals"]["eligible_covers"],
        "small_cover_certificates_recomputed": summary["totals"]["small_covers"],
        "rejected_graphs": EXPECTED_TARGET_GRAPHS,
        "ordered_residue_graphs": len(summary["ordered_residue_indices"]),
        "rejected_indices": summary["rejected_indices"],
        "ordered_residue_indices_sha256": summary[
            "ordered_residue_indices_sha256"
        ],
        "controls": control_rows,
        "dependency_hashes": dependencies,
        "import_independence": {
            "production_conjunction_imported": False,
            "production_profile_imported": False,
            "empty_support_producer_imported": False,
            "checked_by_AST": True,
        },
        "runtime": {
            "command": list(sys.argv),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k6_empty_essential_virtual_k7.json",
    )
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_empty_essential_virtual_k7_verification.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > DEFAULT_WORKERS:
        parser.error(f"workers must lie in 1..{DEFAULT_WORKERS}")
    try:
        result = verify(
            args.report.resolve(),
            args.expected_report_sha256,
            args.output.resolve(),
            args.workers,
        )
    except VerificationInfrastructureAbort as error:
        atomic_json(
            args.output.resolve(),
            {
                "schema": VERIFICATION_SCHEMA,
                "status": "INFRA_ABORT",
                "meaning": "no mathematical rejection is inferred",
                "message": str(error),
                "report_sha256": sha256(args.report.resolve()),
                "verifier_source_sha256": sha256(Path(__file__).resolve()),
            },
        )
        print(str(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": result["status"],
                "promotions_recomputed": result["promotions_recomputed"],
                "small_cover_certificates_recomputed": result[
                    "small_cover_certificates_recomputed"
                ],
                "rejected_graphs": result["rejected_graphs"],
                "ordered_residue_graphs": result["ordered_residue_graphs"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
