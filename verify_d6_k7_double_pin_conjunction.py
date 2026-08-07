#!/usr/bin/env python3
"""Independent verifier for the K7 two-defect double-pin conjunction.

This checker imports neither the builder nor the double-pin locator.  It uses
the pre-existing independent support and sparse-value implementations, a
separately transcribed triangle locator, and replays every graph/seed/cover/
support quantifier on the embedded 155-graph v5 K7 corpus.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import d6_k7_special_h_reference as strict_h
import verify_d6_k7_sparse_value_full as sparse_independent
import verify_d6_k7_support_full as support_independent


ROOT = Path(__file__).resolve().parent
RESIDUE = ROOT / "d6_current_residue_manifest_v5.json"
RESIDUE_VERIFICATION = ROOT / "d6_current_residue_manifest_v5_verification.json"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"

EXPECTED_HASHES = {
    "d6_current_residue_manifest_v5.json": (
        "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5"
    ),
    "d6_current_residue_manifest_v5_verification.json": (
        "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46"
    ),
    "d6_k7_positive_dual_full_certificates.jsonl.gz": (
        "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
    ),
    "d6_k7_rankone_tetrad_full_certificates.jsonl.gz": (
        "f719dbb6492fc20fb3103cb79079567aa2cad163835f97e750412c3d27897c50"
    ),
    "d6_k7_rankone_tetrad_full_decisions.tsv.gz": (
        "2552ce91f52727d9beda60d99d534c1d5a0be3f93686ab321f223d8e0d11dfc1"
    ),
}

EXPECTED_SOURCE_HASHES = {
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_special_h_reference.py": (
        "e649633d7262941c9c4608ea694ae0a5d6f29c146d47060c3dae71a036a9d190"
    ),
    "verify_d6_k7_support_full.py": (
        "c6b9dce47618bd9a4247eb5a9f093682265447729e1868c4dc0677a24f51c628"
    ),
    "verify_d6_k7_sparse_value_full.py": (
        "6c0012e66d0b89e280f3cebed7a61f24cc1c6e128e726663a7b246e7da0edbbc"
    ),
}


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
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def read_failure_keys(path: Path, selected: set[int], field: str):
    output: dict[int, set[tuple[tuple[int, ...], int]]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            item = json.loads(line)
            index = int(item["index"])
            if index not in selected:
                continue
            for witness in item[field]:
                key = (
                    tuple(map(int, witness["seed"])), int(witness["zmask"])
                )
                if key in output[index]:
                    raise ValueError("duplicate inherited witness key")
                output[index].add(key)
    return output


def read_tetrad_rows(selected: set[int]) -> dict[int, dict]:
    rows = {}
    with gzip.open(TETRAD_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            if index in selected:
                rows[index] = row
    if set(rows) != selected:
        raise ValueError("v5 graph missing from tetrad decisions")
    return rows


def independently_all_eligible_covers(
    ladj: Sequence[int], eligible: int
) -> tuple[int, ...]:
    full = (1 << len(ladj)) - 1
    answer = []
    for zmask in range(1 << len(ladj)):
        if zmask & ~eligible or zmask.bit_count() > 7:
            continue
        remaining = full & ~zmask
        if all(
            not (ladj[vertex] & remaining)
            for vertex in support_independent.bit_positions(remaining)
        ):
            answer.append(zmask)
    return tuple(answer)


def inherited_cover_status(
    adj: tuple[int, ...],
    outside: tuple[int, ...],
    defects: tuple[int, ...],
    seed: tuple[int, ...],
    zmask: int,
    baseline: reference.CoverAnalysis,
    prior: set[tuple[tuple[int, ...], int]],
    tetrad: set[tuple[tuple[int, ...], int]],
) -> str:
    if baseline.enhanced_joint_failed:
        return "baseline"
    nvertices = tuple(
        outside[position] for position in range(len(outside))
        if not zmask & (1 << position)
    )
    graph_n = support_independent.independent_induced_graph(adj, nvertices)
    saturating = tuple(reference.clique_masks(
        graph_n, baseline.k_rank_upper
    )) if baseline.k_rank_upper else ()
    if any(
        strict_h.assess_saturating_clique(graph_n, clique).failed
        for clique in saturating
    ):
        return "strict_h"
    key = (seed, zmask)
    if key in prior:
        return "prior_dual"
    if key in tetrad:
        return "tetrad"
    return "passing"


def independent_double_pin(
    graph: Sequence[int], masks: Sequence[int]
) -> dict | None:
    """Separate direct enumeration of the certificate pattern."""

    adjacency = tuple(map(int, graph))
    supports = tuple(map(int, masks))
    if len(adjacency) != len(supports):
        raise ValueError("graph/mask size mismatch")
    pins: dict[tuple[int, int], list[tuple[int, int, int]]] = defaultdict(list)
    for triangle in combinations(range(len(adjacency)), 3):
        first, second, third = triangle
        if not (
            adjacency[first] & (1 << second)
            and adjacency[first] & (1 << third)
            and adjacency[second] & (1 << third)
        ):
            continue
        intersections = (
            supports[first] & supports[second],
            supports[first] & supports[third],
            supports[second] & supports[third],
        )
        if (
            intersections[0]
            and intersections[0] == intersections[1] == intersections[2]
            and intersections[0].bit_count() == 1
        ):
            coordinate = intersections[0].bit_length() - 1
            for vertex in triangle:
                pins[(vertex, coordinate)].append(triangle)
    for central, mask in enumerate(supports):
        if mask.bit_count() != 2:
            continue
        coordinates = tuple(support_independent.bit_positions(mask))
        left = pins.get((central, coordinates[0]), ())
        right = pins.get((central, coordinates[1]), ())
        if left and right:
            return {
                "central": central,
                "central_mask": mask,
                "first_coordinate": coordinates[0],
                "first_triangle": list(left[0]),
                "second_coordinate": coordinates[1],
                "second_triangle": list(right[0]),
            }
    return None


def classify_cover(
    graph_n: tuple[int, ...],
    z_allowed: tuple[int, ...],
    n_allowed: tuple[int, ...],
) -> dict:
    counts: Counter = Counter()
    certificates = []
    passing_witness = None
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
        counts["pre_double_pin_families"] += 1
        pin = independent_double_pin(graph_n, masks)
        if pin is not None:
            counts["double_pin_infeasible_families"] += 1
            certificates.append({
                "z_supports": list(fixed),
                "propagated_masks": list(masks),
                "certificate": pin,
            })
        else:
            counts["double_pin_passing_families"] += 1
            if passing_witness is None:
                passing_witness = {
                    "z_supports": list(fixed),
                    "propagated_masks": list(masks),
                }
    return {
        "status": "PASSING" if passing_witness is not None else "INFEASIBLE",
        "family_counts": dict(sorted(counts.items())),
        "double_pin_certificates": certificates,
        "first_passing_witness": passing_witness,
    }


def verify_graph(payload: tuple) -> dict:
    graph, serialized_prior, serialized_tetrad, expected_current = payload
    prior = {(tuple(seed), int(zmask)) for seed, zmask in serialized_prior}
    tetrad = {(tuple(seed), int(zmask)) for seed, zmask in serialized_tetrad}
    index = int(graph["index"])
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
        outside, defects, ladj, eligible = (
            support_independent.independent_seed_instance(adj, seed_mask)
        )
        outside = tuple(map(int, outside))
        defects = tuple(map(int, defects))
        seed = tuple(support_independent.bit_positions(seed_mask))
        covers = independently_all_eligible_covers(ladj, eligible)
        if covers != tuple(reference.eligible_covers(ladj, eligible)):
            raise ValueError("independent eligible-cover enumeration differs")
        total_term_rank = reference.matching_size(defects)
        cover_records = []
        for zmask in covers:
            totals["eligible_covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask, support_solver, zero_forcing,
                total_term_rank, clique_solver,
            )
            status = inherited_cover_status(
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
            totals[f"double_pin_{classified['status'].lower()}_covers"] += 1
            totals.update(classified["family_counts"])
            cover_records.append({
                "seed": list(seed),
                "zmask": zmask,
                **classified,
            })
        if not cover_records:
            raise ValueError("v5 graph has no inherited-passing cover at a seed")
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
    if observed_current != int(expected_current):
        raise ValueError("independent inherited-cover count mismatch")
    rejecting = next(
        (seed for seed in seed_records if seed["status"] == "INFEASIBLE"), None
    )
    return {
        "index": index,
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
        default=ROOT / "d6_k7_double_pin_conjunction_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_double_pin_conjunction_verification.json",
    )
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    report_hash = sha256(args.report)
    if report_hash != args.report_sha256:
        raise ValueError("report differs from explicit SHA-256 pin")
    report = json.loads(args.report.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind") != "d6_k7_two_defect_double_pin_seed_conjunction"
        or report.get("status") != "COMPLETE"
        or report.get("upstream_sha256") != dict(sorted(EXPECTED_HASHES.items()))
    ):
        raise ValueError("unexpected report schema or upstream roots")
    for name, expected in {**EXPECTED_HASHES, **EXPECTED_SOURCE_HASHES}.items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"root hash mismatch: {name}")
    for name, expected in report["source_sha256"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"reported production source changed: {name}")
    verification = json.loads(RESIDUE_VERIFICATION.read_text(encoding="utf-8"))
    if verification.get("status") != "PASS":
        raise ValueError("v5 residue verification is not PASS")
    residue = json.loads(RESIDUE.read_text(encoding="utf-8"))
    graphs = residue["classes"]["K7"]["graphs"]
    indices = [int(graph["index"]) for graph in graphs]
    if len(graphs) != 155 or indices != residue["classes"]["K7"]["indices"]:
        raise ValueError("v5 K7 corpus population/order mismatch")
    selected = set(indices)
    prior = read_failure_keys(PRIOR_CERTIFICATES, selected, "dual_failure_witnesses")
    tetrad = read_failure_keys(TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses")
    tetrad_rows = read_tetrad_rows(selected)
    payloads = [(
        graph,
        [[list(seed), zmask] for seed, zmask in sorted(prior.get(index, set()))],
        [[list(seed), zmask] for seed, zmask in sorted(tetrad.get(index, set()))],
        int(tetrad_rows[index]["tetrad_passing_covers"]),
    ) for graph, index in zip(graphs, indices, strict=True)]

    started = time.monotonic()
    if args.workers == 1:
        records = list(map(verify_graph, payloads))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            records = list(executor.map(verify_graph, payloads, chunksize=1))
    if records != report["records"]:
        raise ValueError("independent full quantifier records differ from report")
    rejected = [record["index"] for record in records if record["decision"] == "REJECTED"]
    survivors = [record["index"] for record in records if record["decision"] == "SURVIVOR"]
    totals: Counter = Counter()
    for record in records:
        totals.update(record["counts"])
    summary = report["summary"]
    if (
        len(rejected) != 130
        or len(survivors) != 25
        or summary["rejected_indices"] != rejected
        or summary["rejected_indices_sha256"] != stable_hash(rejected)
        or summary["survivor_indices"] != survivors
        or summary["survivor_indices_sha256"] != stable_hash(survivors)
        or summary["totals"] != dict(sorted(totals.items()))
    ):
        raise ValueError("independent aggregate partition differs from report")
    output = {
        "schema": 1,
        "kind": "d6_k7_two_defect_double_pin_seed_conjunction_verification",
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": report_hash},
        "checked": {
            "graphs": len(graphs),
            "K7_seeds": totals["seeds"],
            "eligible_covers": totals["eligible_covers"],
            "current_covers": totals["current_covers"],
            "labeled_z_families": totals["labeled_z_families"],
            "pre_double_pin_families": totals["pre_double_pin_families"],
            "double_pin_certificates": totals["double_pin_infeasible_families"],
            "rejected_graphs": len(rejected),
            "survivors": len(survivors),
        },
        "ordered_rejected_indices_sha256": stable_hash(rejected),
        "ordered_survivor_indices_sha256": stable_hash(survivors),
        "root_sha256": dict(sorted(EXPECTED_HASHES.items())),
        "source_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            **EXPECTED_SOURCE_HASHES,
        },
        "execution": {
            "workers": args.workers,
            "logical_cpus": os.cpu_count(),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
        "claim": (
            "A separate support enumerator, propagation/value checker, and "
            "double-pin locator reproduced all 155 records and the exact "
            "130/25 graph partition."
        ),
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
