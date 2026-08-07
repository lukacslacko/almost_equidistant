#!/usr/bin/env python3
"""Bounded conjunction of virtual-K7 branches with the K6 empty Hall budget.

This is a pilot only.  It reads the frozen 20-graph virtual-K7 pilot, reruns a
stronger double-apex screen on its surviving candidate pairs, then restricts
the empty-singleton alternatives in the current K6 empty-support evaluator.
The all-nonempty choice is always retained.  No full-822 campaign is exposed.
"""

from __future__ import annotations

import argparse
import ctypes
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
import d6_k6_empty_support_budget as empty_budget
import d6_k6_psd_z_hereditary as hall
import d6_k6_psd_zmatrix as ranks
import probe_d6_k6_virtual_k7 as virtual_k7
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    add_edge,
    build_instance,
    clique_masks,
    support_matching,
    vertices,
)
from d6_k6_normal_inertia import InertiaCache


ROOT = Path(__file__).resolve().parent
VIRTUAL_REPORT = ROOT / "d6_k6_virtual_k7_pilot_report.json"
EXPECTED_VIRTUAL_REPORT_SHA256 = (
    "9f43436fddd5b1aa3511318251cdd992c0621d76807200166ffd49e7fb872148"
)
EXPECTED_VIRTUAL_SOURCE_SHA256 = (
    "531f4d123f249de03438e411f03fb7c70a591f3398b2a162b92158a31fd655b0"
)
EXPECTED_EMPTY_SOURCE_SHA256 = (
    "284b8c3bda4d2a43581665ade3d52d9295b29d460418454e3a760897cabc5fac"
)
EXPECTED_SELECTION_SHA256 = (
    "a45245b9a7fe3b73de131698f1f42d6f68c712e68078a2055c10d2b0d24fd37e"
)
EXPECTED_SELECTION_SIZE = 20
QOS_CLASS_USER_INITIATED = 0x19


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


def worker_init() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    if sys.platform != "darwin":
        return
    libc = ctypes.CDLL(None)
    setter = libc.pthread_set_qos_class_self_np
    setter.argtypes = (ctypes.c_uint, ctypes.c_int)
    setter.restype = ctypes.c_int
    error = int(setter(QOS_CLASS_USER_INITIATED, 0))
    if error:
        raise OSError(error, "pthread_set_qos_class_self_np failed")


def load_inputs() -> tuple[list[dict], dict[int, dict]]:
    if sha256(VIRTUAL_REPORT) != EXPECTED_VIRTUAL_REPORT_SHA256:
        raise ValueError("virtual-K7 pilot report hash changed")
    if sha256(ROOT / "probe_d6_k6_virtual_k7.py") != EXPECTED_VIRTUAL_SOURCE_SHA256:
        raise ValueError("virtual-K7 pilot source hash changed")
    if sha256(ROOT / "d6_k6_empty_support_budget.py") != EXPECTED_EMPTY_SOURCE_SHA256:
        raise ValueError("current empty-support evaluator source hash changed")

    report = json.loads(VIRTUAL_REPORT.read_text(encoding="utf-8"))
    selection = [int(index) for index in report["selection"]]
    if (
        report.get("status") != "PILOT_COMPLETE"
        or len(selection) != EXPECTED_SELECTION_SIZE
        or stable_hash(selection) != EXPECTED_SELECTION_SHA256
        or report.get("semantics", {}).get("original_graph_rejections_claimed") != 0
    ):
        raise ValueError("virtual-K7 pilot selection/semantics changed")

    _, residue, _, _ = empty_budget.verify_parent_input()
    by_index = {int(record["index"]): record for record in residue}
    selected = [by_index[index] for index in selection]
    virtual_by_index = {
        int(row["index"]): row for row in report["graph_results"]
    }
    if list(virtual_by_index) != selection:
        raise ValueError("virtual report graph order changed")
    return selected, virtual_by_index


def pair_key(first: int, second: int) -> tuple[int, int]:
    return (first, second) if first < second else (second, first)


def augment_both_apices(
    adjacency: Sequence[int], seed: Sequence[int], first: int, second: int
) -> tuple[int, ...]:
    if len(seed) != 6 or first in seed or second in seed or first == second:
        raise ValueError("bad double-apex promotion")
    if adjacency[first] & (1 << second):
        raise ValueError("opposite empty apices must be a candidate nonedge")
    augmented = list(map(int, adjacency))
    for apex in (first, second):
        for q in seed:
            if not (augmented[apex] & (1 << q)):
                add_edge(augmented, apex, q)
    # Deliberately leave the apex pair unchanged as a candidate nonedge.
    if augmented[first] & (1 << second):
        raise AssertionError("double promotion accidentally added the apex pair")
    return tuple(augmented)


def screen_double_apex_pairs(
    adjacency: Sequence[int], seed: Sequence[int], candidate_pairs: Sequence[Sequence[int]]
) -> tuple[set[tuple[int, int]], list[dict]]:
    """Run both promoted K7 seeds after simultaneously adding both apices."""

    support_solver = virtual_k7.k7.SupportSolver()
    zero_forcing = virtual_k7.k7.ZeroForcingSolver()
    clique_solver = virtual_k7.k7.CliqueStructureSolver()
    survivors: set[tuple[int, int]] = set()
    rows = []
    for raw_pair in candidate_pairs:
        if len(raw_pair) != 2:
            raise ValueError("bad virtual-pair record")
        first, second = pair_key(int(raw_pair[0]), int(raw_pair[1]))
        augmented = augment_both_apices(adjacency, seed, first, second)
        first_result = virtual_k7.analyze_virtual_seed(
            augmented,
            seed,
            first,
            support_solver,
            zero_forcing,
            clique_solver,
        )
        second_result = virtual_k7.analyze_virtual_seed(
            augmented,
            seed,
            second,
            support_solver,
            zero_forcing,
            clique_solver,
        )
        passes = (
            first_result["status"] == virtual_k7.STATUS_SURVIVOR
            and second_result["status"] == virtual_k7.STATUS_SURVIVOR
        )
        if passes:
            survivors.add((first, second))
        rows.append(
            {
                "pair": [first, second],
                "first_seed_status": first_result["status"],
                "second_seed_status": second_result["status"],
                "survives": passes,
                "first_seed": virtual_k7.compact_augmentation(first_result),
                "second_seed": virtual_k7.compact_augmentation(second_result),
            }
        )
    return survivors, rows


def orientation_empty_options(
    instance: K6LorentzInstance,
    z0: int,
    side_a,
    side_b,
    subset_rank: hall_parent.PositiveSubsetRank,
    negative_side: str,
    viable_apices: set[int],
) -> dict:
    """Enumerate every Hall-passing none/singleton choice on one orientation."""

    groups_a, _ = hall_parent.make_side_groups(instance, "A", side_a, subset_rank)
    groups_b, _ = hall_parent.make_side_groups(instance, "B", side_b, subset_rank)
    if negative_side == "A":
        negative_groups, negative_rank = groups_a, side_a
    elif negative_side == "B":
        negative_groups, negative_rank = groups_b, side_b
    else:
        raise ValueError("negative_side must be A or B")
    zero_capable = empty_budget.zero_capable_singletons(
        instance, negative_groups, negative_rank
    )
    absolute_by_group = {
        item.group_index: item.absolute_vertex for item in zero_capable
    }
    z0_groups = tuple(
        hall.singleton_group(
            f"Z0:{instance.outside[local]}",
            instance.outside[local],
            instance.defects[local],
        )
        for local in vertices(z0)
    )
    allowed_groups = tuple(
        item.group_index
        for item in zero_capable
        if item.absolute_vertex in viable_apices
    )
    passing: list[int | None] = []
    rows = []
    for chosen in (None, *allowed_groups):
        if negative_side == "A":
            test_a = empty_budget.force_singleton_status(
                groups_a, zero_capable, chosen
            )
            test_b = groups_b
        else:
            test_a = groups_a
            test_b = empty_budget.force_singleton_status(
                groups_b, zero_capable, chosen
            )
        hall_result = hall.check_groups(test_a + test_b + z0_groups)
        absolute = None if chosen is None else absolute_by_group[chosen]
        rows.append({"empty": absolute, "Hall_passed": hall_result.passed})
        if hall_result.passed:
            passing.append(absolute)
    return {
        "negative_singletons": [item.absolute_vertex for item in zero_capable],
        "virtual_viable_negative_singletons": [
            item.absolute_vertex
            for item in zero_capable
            if item.absolute_vertex in viable_apices
        ],
        "passing_empty_choices": passing,
        "Hall_rows": rows,
    }


def component_empty_options(
    adjacency: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
    subset_rank: hall_parent.PositiveSubsetRank,
    viable_apices: set[int],
) -> dict:
    """Return every none/singleton choice surviving this Lorentz component."""

    graph_a, absolute_a = induced_required_graph(adjacency, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adjacency, instance, component.side_b)
    alternatives = []
    choices: set[int | None] = set()
    for name, sign_a, sign_b, negative_side in (
        ("A_positive", "positive", "negative", "B"),
        ("A_negative", "negative", "positive", "A"),
    ):
        side_a = ranks.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = ranks.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        row = orientation_empty_options(
            instance,
            z0,
            side_a,
            side_b,
            subset_rank,
            negative_side,
            viable_apices,
        )
        alternatives.append({"kind": name, **row})
        choices.update(row["passing_empty_choices"])

    light = z0 | component.component
    light_passed = (
        light.bit_count() <= COORDINATES
        and support_matching(light, instance.defects) is not None
    )
    alternatives.append(
        {
            "kind": "lightlike",
            "passed": light_passed,
            "orthonormal_vectors": light.bit_count(),
        }
    )
    if light_passed:
        choices.add(None)
    ordered = sorted((choice for choice in choices if choice is not None))
    if None in choices:
        ordered = [None, *ordered]
    return {"choices": ordered, "alternatives": alternatives}


def extend_empty_states(
    states: set[tuple[int, ...]],
    options: Sequence[int | None],
    surviving_pairs: set[tuple[int, int]],
) -> set[tuple[int, ...]]:
    """Conjoin component choices with the global opposite-apex constraint."""

    following: set[tuple[int, ...]] = set()
    for state in states:
        for option in options:
            if option is None:
                following.add(state)
                continue
            if option in state:
                raise AssertionError("one vertex appeared in two Lorentz components")
            chosen = tuple(sorted((*state, option)))
            if len(chosen) > 2:
                continue
            if len(chosen) == 2 and pair_key(*chosen) not in surviving_pairs:
                continue
            following.add(chosen)
    return following


def solve_seed_conjunction(
    adjacency: Sequence[int],
    instance: K6LorentzInstance,
    viable_apices: set[int],
    surviving_pairs: set[tuple[int, int]],
) -> dict:
    """Exhaust the common Z0, Hall choices, B_Q, and double-apex pairs."""

    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    subset_rank = hall_parent.PositiveSubsetRank(inertia, zero_forcing)
    rows = []
    counts: Counter[str] = Counter()
    for z0 in ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0"] += 1
        if support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            rows.append(
                {
                    "Z0": hall_parent.frozen.absolute_vertices(instance, z0),
                    "matchable": False,
                    "terminal": "Z0_ALLOWED_SUPPORT_UNMATCHABLE",
                }
            )
            continue

        states: set[tuple[int, ...]] = {()}
        component_rows = []
        terminal = None
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                counts["nonbipartite_components_ignored"] += 1
                continue
            counts["bipartite_components"] += 1
            decision = component_empty_options(
                adjacency,
                instance,
                z0,
                component,
                inertia,
                zero_forcing,
                subset_rank,
                viable_apices,
            )
            counts["orientation_Hall_systems"] += 2
            counts["component_empty_options"] += len(decision["choices"])
            before = sorted(map(list, states))
            states = extend_empty_states(
                states, decision["choices"], surviving_pairs
            )
            component_rows.append(
                {
                    "vertices": hall_parent.frozen.absolute_vertices(
                        instance, component.component
                    ),
                    "empty_options": decision["choices"],
                    "states_before": before,
                    "states_after": sorted(map(list, states)),
                    "alternatives": decision["alternatives"],
                }
            )
            if not states:
                terminal = "NO_GLOBALLY_COMPATIBLE_EMPTY_ASSIGNMENT"
                break
        if states:
            chosen = min(states, key=lambda state: (len(state), state))
            return {
                "feasible": True,
                "counts": dict(counts),
                "chosen_Z0": hall_parent.frozen.absolute_vertices(instance, z0),
                "chosen_empty_vertices": list(chosen),
                "certificate_rows": None,
            }
        rows.append(
            {
                "Z0": hall_parent.frozen.absolute_vertices(instance, z0),
                "matchable": True,
                "terminal": terminal,
                "components": component_rows,
            }
        )
    return {
        "feasible": False,
        "counts": dict(counts),
        "chosen_Z0": None,
        "chosen_empty_vertices": None,
        "certificate_rows": rows,
    }


def evaluate_record(payload: tuple[dict, dict]) -> dict:
    record, virtual_row = payload
    started = time.perf_counter()
    adjacency = tuple(map(int, record["adjacency"]))
    summaries = {
        int(row["seed_mask"]): row for row in virtual_row["seed_summaries"]
    }
    seed_masks = tuple(clique_masks(adjacency, COORDINATES))
    if set(summaries) != set(seed_masks):
        raise ValueError("virtual seed summaries do not match K6 enumeration")

    pair_rows = []
    pair_survivors_by_seed: dict[int, set[tuple[int, int]]] = {}
    for seed_mask in seed_masks:
        seed = tuple(vertices(seed_mask))
        summary = summaries[seed_mask]
        candidate_pairs = summary["possible_opposite_empty_pairs"]
        survivors, rows = screen_double_apex_pairs(
            adjacency, seed, candidate_pairs
        )
        pair_survivors_by_seed[seed_mask] = survivors
        pair_rows.extend(
            {"seed_mask": seed_mask, **row} for row in rows
        )

    base_impossible = []
    conjunction_impossible = []
    seed_rows = []
    first_certificate = None
    total_counts: Counter[str] = Counter()
    for seed_mask in seed_masks:
        seed = tuple(vertices(seed_mask))
        instance = build_instance(adjacency, seed)
        base = empty_budget.solve_seed(adjacency, instance)
        viable = {int(vertex) for vertex in summaries[seed_mask]["viable_apices"]}
        conjunction = solve_seed_conjunction(
            adjacency,
            instance,
            viable,
            pair_survivors_by_seed[seed_mask],
        )
        total_counts.update(conjunction["counts"])
        if not base["feasible"]:
            base_impossible.append(seed_mask)
        if not conjunction["feasible"]:
            conjunction_impossible.append(seed_mask)
            if first_certificate is None:
                first_certificate = {
                    "seed": list(seed),
                    "seed_mask": seed_mask,
                    "viable_apices": sorted(viable),
                    "surviving_pairs": [
                        list(pair) for pair in sorted(pair_survivors_by_seed[seed_mask])
                    ],
                    "choices": conjunction["certificate_rows"],
                }
        if conjunction["feasible"] and not base["feasible"]:
            raise AssertionError("restricted conjunction passed a base-impossible seed")
        seed_rows.append(
            {
                "seed_mask": seed_mask,
                "base_feasible": base["feasible"],
                "conjunction_feasible": conjunction["feasible"],
                "viable_apex_count": len(viable),
                "double_apex_pair_count": len(pair_survivors_by_seed[seed_mask]),
                "chosen_Z0": conjunction["chosen_Z0"],
                "chosen_empty_vertices": conjunction["chosen_empty_vertices"],
            }
        )

    base_rejected = bool(base_impossible)
    conjunction_rejected = bool(conjunction_impossible)
    return {
        "index": int(record["index"]),
        "k6_seeds": len(seed_masks),
        "base_rejected": base_rejected,
        "conjunction_rejected": conjunction_rejected,
        "marginal_rejected": conjunction_rejected and not base_rejected,
        "base_impossible_seed_masks": base_impossible,
        "conjunction_impossible_seed_masks": conjunction_impossible,
        "marginal_impossible_seed_masks": [
            mask for mask in conjunction_impossible if mask not in set(base_impossible)
        ],
        "singleton_viable_apices": sum(
            len(summary["viable_apices"]) for summary in summaries.values()
        ),
        "double_apex_candidate_pairs": len(pair_rows),
        "double_apex_surviving_pairs": sum(
            row["survives"] for row in pair_rows
        ),
        "pair_rows": pair_rows,
        "seed_rows": seed_rows,
        "conjunction_counts": dict(total_counts),
        "first_conjunction_certificate": first_certificate,
        "elapsed_seconds": time.perf_counter() - started,
    }


def aggregate(results: Sequence[dict]) -> dict:
    totals: Counter[str] = Counter()
    pair_statuses: Counter[str] = Counter()
    conjunction_counts: Counter[str] = Counter()
    for result in results:
        totals.update(
            graphs=1,
            k6_seeds=result["k6_seeds"],
            base_graph_rejections=int(result["base_rejected"]),
            conjunction_graph_rejections=int(result["conjunction_rejected"]),
            marginal_graph_rejections=int(result["marginal_rejected"]),
            base_impossible_seeds=len(result["base_impossible_seed_masks"]),
            conjunction_impossible_seeds=len(
                result["conjunction_impossible_seed_masks"]
            ),
            marginal_impossible_seeds=len(result["marginal_impossible_seed_masks"]),
            singleton_viable_apices=result["singleton_viable_apices"],
            double_apex_candidate_pairs=result["double_apex_candidate_pairs"],
            double_apex_surviving_pairs=result["double_apex_surviving_pairs"],
        )
        conjunction_counts.update(result["conjunction_counts"])
        for row in result["pair_rows"]:
            pair_statuses[f"first_{row['first_seed_status']}"] += 1
            pair_statuses[f"second_{row['second_seed_status']}"] += 1
            pair_statuses["pair_survives" if row["survives"] else "pair_eliminated"] += 1
    return {
        "totals": dict(totals),
        "double_apex_status_counts": dict(pair_statuses),
        "conjunction_effort": dict(conjunction_counts),
        "base_rejected_indices": [r["index"] for r in results if r["base_rejected"]],
        "conjunction_rejected_indices": [
            r["index"] for r in results if r["conjunction_rejected"]
        ],
        "marginal_rejected_indices": [
            r["index"] for r in results if r["marginal_rejected"]
        ],
    }


def run(output: Path, workers: int) -> dict:
    records, virtual_by_index = load_inputs()
    payloads = [(record, virtual_by_index[int(record["index"])]) for record in records]
    started = time.perf_counter()
    if workers == 1:
        results = [evaluate_record(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=worker_init) as pool:
            results = list(pool.map(evaluate_record, payloads, chunksize=1))
    summary = aggregate(results)
    decision_rows = [
        {
            "index": result["index"],
            "base_rejected": result["base_rejected"],
            "conjunction_rejected": result["conjunction_rejected"],
            "marginal_rejected": result["marginal_rejected"],
            "conjunction_impossible_seed_masks": result[
                "conjunction_impossible_seed_masks"
            ],
        }
        for result in results
    ]
    report = {
        "schema": "d6-k6-virtual-k7-empty-conjunction-pilot-v1",
        "status": "PILOT_COMPLETE",
        "selection": [int(record["index"]) for record in records],
        "selection_sha256": stable_hash([int(record["index"]) for record in records]),
        "aggregate": summary,
        "decision_rows_sha256": stable_hash(decision_rows),
        "graph_results": results,
        "runtime": {
            "command": (
                f"{sys.executable} {Path(__file__).name} --workers {workers} "
                f"--output {output}"
            ),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "sources": {
            Path(__file__).name: sha256(Path(__file__)),
            VIRTUAL_REPORT.name: sha256(VIRTUAL_REPORT),
            "probe_d6_k6_virtual_k7.py": sha256(
                ROOT / "probe_d6_k6_virtual_k7.py"
            ),
            "d6_k6_empty_support_budget.py": sha256(
                ROOT / "d6_k6_empty_support_budget.py"
            ),
        },
        "semantics": {
            "all_nonempty_branch_retained": True,
            "empty_singletons_restricted_to_B_Q": True,
            "maximum_empty_vertices": 2,
            "two_empty_vertices": (
                "candidate nonedge and both promoted K7 seeds survive after "
                "simultaneously adding both apex-to-Q stars"
            ),
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defect_coordinates": "upper bounds and may be zero",
            "nonbipartite_components": "ignored conservatively as in base layer",
            "pilot_only_no_theorem_credit": True,
            "survivor": "filter non-rejection only",
        },
    }
    atomic_json(output, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_virtual_k7_empty_conjunction_pilot_report.json",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    report = run(args.output, args.workers)
    print(json.dumps(report["aggregate"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
