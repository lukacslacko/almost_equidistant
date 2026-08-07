#!/usr/bin/env python3
"""Independent verifier for the graph-3919831 odd-cycle/full-pin increment."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import verify_d6_k7_double_pin_conjunction as base_independent
import verify_d6_k7_sparse_value_full as sparse_independent
import verify_d6_k7_support_full as support_independent


ROOT = Path(__file__).resolve().parent
TARGET = 3_919_831
BASE_REPORT = ROOT / "d6_k7_double_pin_conjunction_report.json"
BASE_VERIFICATION = ROOT / "d6_k7_double_pin_conjunction_verification.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def coordinate_graph(required: Sequence[int], masks: Sequence[int], coordinate: int):
    bit = 1 << coordinate
    answer = [0] * len(required)
    for first in range(len(required)):
        for second in range(first):
            if required[first] & (1 << second) and masks[first] & masks[second] == bit:
                answer[first] |= 1 << second
                answer[second] |= 1 << first
    return tuple(answer)


def nonbipartite_components(graph: Sequence[int]):
    unseen = set(range(len(graph)))
    output = []
    while unseen:
        root = min(unseen)
        colors = {root: 0}
        queue = [root]
        unseen.remove(root)
        conflict = None
        for vertex in queue:
            for neighbour in support_independent.bit_positions(graph[vertex]):
                if neighbour not in colors:
                    colors[neighbour] = colors[vertex] ^ 1
                    unseen.discard(neighbour)
                    queue.append(neighbour)
                elif colors[neighbour] == colors[vertex] and conflict is None:
                    conflict = (min(vertex, neighbour), max(vertex, neighbour))
        if conflict is not None:
            output.append((tuple(sorted(colors)), conflict))
    return tuple(output)


def independent_full_pin(graph: Sequence[int], masks: Sequence[int]) -> dict | None:
    pinned = [0] * len(graph)
    data = {}
    for coordinate in range(7):
        local = coordinate_graph(graph, masks, coordinate)
        for component, conflict in nonbipartite_components(local):
            for vertex in component:
                pinned[vertex] |= 1 << coordinate
                data[(vertex, coordinate)] = (component, conflict)
    for central, mask in enumerate(masks):
        if pinned[central] & mask == mask:
            return {
                "central": central,
                "central_mask": int(mask),
                "coordinate_components": [
                    {
                        "coordinate": coordinate,
                        "component": list(data[(central, coordinate)][0]),
                        "conflict_edge": list(data[(central, coordinate)][1]),
                    }
                    for coordinate in support_independent.bit_positions(mask)
                ],
            }
    return None


def classify_cover(graph_n, z_allowed, n_allowed):
    counts: Counter = Counter()
    certificates = []
    passing = None
    for fixed in support_independent.independent_labeled_support_families(z_allowed):
        counts["labeled_z_families"] += 1
        masks, _ = sparse_independent.propagated_masks(n_allowed, fixed)
        failure = sparse_independent.simple_propagation_failure(graph_n, masks)
        if failure is not None:
            counts[f"propagation_failure:{failure}"] += 1
            continue
        sparse = sparse_independent.independent_small_support_check(graph_n, masks)
        if not sparse.feasible:
            counts["sparse_value_failure"] += 1
            continue
        counts["pre_full_pin_families"] += 1
        certificate = independent_full_pin(graph_n, masks)
        if certificate is None:
            counts["full_pin_passing_families"] += 1
            if passing is None:
                passing = {
                    "z_supports": list(fixed),
                    "propagated_masks": list(masks),
                }
        else:
            counts["full_pin_infeasible_families"] += 1
            certificates.append({
                "z_supports": list(fixed),
                "propagated_masks": list(masks),
                "certificate": certificate,
            })
    return {
        "status": "PASSING" if passing is not None else "INFEASIBLE",
        "family_counts": dict(sorted(counts.items())),
        "full_pin_certificates": certificates,
        "first_passing_witness": passing,
    }


def independently_evaluate_target(graph, prior, tetrad, expected_current):
    adj = tuple(map(int, graph["adjacency"]))
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    totals: Counter = Counter()
    seed_records = []
    observed_current = 0
    for seed_mask in reference.clique_masks(adj, 7):
        totals["seeds"] += 1
        outside, defects, ladj, eligible = support_independent.independent_seed_instance(
            adj, seed_mask
        )
        outside = tuple(map(int, outside))
        defects = tuple(map(int, defects))
        seed = tuple(support_independent.bit_positions(seed_mask))
        covers = base_independent.independently_all_eligible_covers(ladj, eligible)
        total_term_rank = reference.matching_size(defects)
        cover_records = []
        for zmask in covers:
            totals["eligible_covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask, support_solver, zero_forcing,
                total_term_rank, clique_solver,
            )
            status = base_independent.inherited_cover_status(
                adj, outside, defects, seed, zmask, baseline, prior, tetrad
            )
            totals[f"inherited_cover_{status}"] += 1
            if status != "passing":
                continue
            observed_current += 1
            totals["current_covers"] += 1
            zvertices = tuple(support_independent.bit_positions(zmask))
            nlocal = tuple(
                position for position in range(len(outside))
                if not zmask & (1 << position)
            )
            graph_n = support_independent.independent_induced_graph(
                adj, tuple(outside[position] for position in nlocal)
            )
            classified = classify_cover(
                graph_n,
                tuple(defects[position] for position in zvertices),
                tuple(defects[position] for position in nlocal),
            )
            totals.update(classified["family_counts"])
            totals[f"full_pin_{classified['status'].lower()}_covers"] += 1
            cover_records.append({"seed": list(seed), "zmask": zmask, **classified})
        seed_status = (
            "PASSING" if any(cover["status"] == "PASSING" for cover in cover_records)
            else "INFEASIBLE"
        )
        seed_records.append({
            "seed": list(seed),
            "seed_mask": seed_mask,
            "status": seed_status,
            "current_covers": cover_records,
        })
    if observed_current != expected_current:
        raise ValueError("independent current-cover count mismatch")
    rejecting = next(
        (seed for seed in seed_records if seed["status"] == "INFEASIBLE"), None
    )
    return {
        "index": TARGET,
        "decision": "REJECTED" if rejecting is not None else "SURVIVOR",
        "counts": dict(sorted(totals.items())),
        "first_rejecting_seed_mask": (
            int(rejecting["seed_mask"]) if rejecting is not None else 0
        ),
        "seeds": seed_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_full_pin_increment_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--base-report-sha256", required=True)
    parser.add_argument("--base-verification-sha256", required=True)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_full_pin_increment_verification.json",
    )
    args = parser.parse_args()
    if sha256(args.report) != args.report_sha256:
        raise ValueError("increment report differs from explicit hash pin")
    if sha256(BASE_REPORT) != args.base_report_sha256:
        raise ValueError("base report differs from explicit hash pin")
    if sha256(BASE_VERIFICATION) != args.base_verification_sha256:
        raise ValueError("base verification differs from explicit hash pin")
    report = json.loads(args.report.read_text(encoding="utf-8"))
    base_report = json.loads(BASE_REPORT.read_text(encoding="utf-8"))
    base_verification = json.loads(BASE_VERIFICATION.read_text(encoding="utf-8"))
    expected_base = {
        "d6_k7_double_pin_conjunction_report.json": args.base_report_sha256,
        "d6_k7_double_pin_conjunction_verification.json": args.base_verification_sha256,
    }
    if (
        report.get("kind") != "d6_k7_full_pin_odd_cycle_increment"
        or report.get("status") != "COMPLETE"
        or report.get("base_sha256") != dict(sorted(expected_base.items()))
        or base_verification.get("status") != "PASS"
        or base_verification.get("report", {}).get("sha256")
        != args.base_report_sha256
    ):
        raise ValueError("increment/base report gates failed")
    for name, expected in report["source_sha256"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"reported production source changed: {name}")
    residue = json.loads(base_independent.RESIDUE.read_text(encoding="utf-8"))
    graph = next(
        item for item in residue["classes"]["K7"]["graphs"]
        if int(item["index"]) == TARGET
    )
    selected = {TARGET}
    prior = base_independent.read_failure_keys(
        base_independent.PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
    ).get(TARGET, set())
    tetrad = base_independent.read_failure_keys(
        base_independent.TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
    ).get(TARGET, set())
    tetrad_row = base_independent.read_tetrad_rows(selected)[TARGET]
    started = time.monotonic()
    record = independently_evaluate_target(
        graph, prior, tetrad, int(tetrad_row["tetrad_passing_covers"])
    )
    if record != report["record"] or record["decision"] != "REJECTED":
        raise ValueError("independent target quantifier differs from report")
    base_survivors = list(base_report["summary"]["survivor_indices"])
    survivors = [index for index in base_survivors if index != TARGET]
    if (
        report["base_survivors"] != base_survivors
        or report["base_survivors_sha256"] != stable_hash(base_survivors)
        or report["exact_survivors"] != survivors
        or report["exact_survivors_sha256"] != stable_hash(survivors)
        or len(survivors) != 24
    ):
        raise ValueError("independent 25-to-24 residue update differs")
    output = {
        "schema": 1,
        "kind": "d6_k7_full_pin_odd_cycle_increment_verification",
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": args.report_sha256},
        "base_sha256": dict(sorted(expected_base.items())),
        "checked": {
            "target": TARGET,
            "seeds": record["counts"]["seeds"],
            "eligible_covers": record["counts"]["eligible_covers"],
            "current_covers": record["counts"]["current_covers"],
            "labeled_z_families": record["counts"]["labeled_z_families"],
            "pre_full_pin_families": record["counts"]["pre_full_pin_families"],
            "full_pin_certificates": record["counts"]["full_pin_infeasible_families"],
            "exact_survivors": len(survivors),
        },
        "exact_survivors_sha256": stable_hash(survivors),
        "source_sha256": {Path(__file__).name: sha256(Path(__file__))},
        "execution": {
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
    }
    atomic_json(args.output, output)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "status": "PASS",
        "checked": output["checked"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
