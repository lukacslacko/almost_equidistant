#!/usr/bin/env python3
"""Build the exact odd-cycle/full-pin increment on graph 3919831."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from pathlib import Path

import build_d6_k7_double_pin_conjunction as base
import d6_k7_full_pin_odd_cycle as full_pin
import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation
import run_d6_k7_support_capacity_pilot as inherited


ROOT = Path(__file__).resolve().parent
BASE_REPORT = ROOT / "d6_k7_double_pin_conjunction_report.json"
BASE_VERIFICATION = ROOT / "d6_k7_double_pin_conjunction_verification.json"
TARGET = 3_919_831


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


def evaluate_target(graph: dict, prior: set, tetrad: set, expected_current: int) -> dict:
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
        seed, outside, defects, ladj, eligible = reference.seed_instance(adj, seed_mask)
        total_term_rank = reference.matching_size(defects)
        cover_records = []
        for zmask in reference.eligible_covers(ladj, eligible):
            totals["eligible_covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask, support_solver, zero_forcing,
                total_term_rank, clique_solver,
            )
            status = inherited.current_cover_status(
                adj, outside, defects, zmask, baseline,
                (tuple(seed), zmask), prior, tetrad,
            )
            totals[f"inherited_cover_{status}"] += 1
            if status != "passing":
                continue
            observed_current += 1
            totals["current_covers"] += 1
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex for vertex in range(len(outside)) if not zmask & (1 << vertex)
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nvertices]
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            counts: Counter = Counter()
            certificates = []
            passing = None
            for z_supports in propagation.labeled_support_families(z_allowed):
                counts["labeled_z_families"] += 1
                propagated = propagation.analyze_support_assignment(
                    graph_n, z_supports, n_allowed, zero_forcing, clique_solver
                )
                if propagated.failure is not None:
                    counts[f"propagation_failure:{propagated.failure}"] += 1
                    continue
                sparse = sparse_value.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                if not sparse.feasible:
                    counts["sparse_value_failure"] += 1
                    continue
                counts["pre_full_pin_families"] += 1
                certificate = full_pin.find_full_pin(
                    graph_n, propagated.propagated_masks
                )
                if certificate is None:
                    counts["full_pin_passing_families"] += 1
                    if passing is None:
                        passing = {
                            "z_supports": list(z_supports),
                            "propagated_masks": list(propagated.propagated_masks),
                        }
                else:
                    full_pin.verify_certificate(
                        graph_n, propagated.propagated_masks, certificate
                    )
                    counts["full_pin_infeasible_families"] += 1
                    certificates.append({
                        "z_supports": list(z_supports),
                        "propagated_masks": list(propagated.propagated_masks),
                        "certificate": full_pin.certificate_json(certificate),
                    })
            totals.update(counts)
            cover_status = "PASSING" if passing is not None else "INFEASIBLE"
            totals[f"full_pin_{cover_status.lower()}_covers"] += 1
            cover_records.append({
                "seed": seed,
                "zmask": zmask,
                "status": cover_status,
                "family_counts": dict(sorted(counts.items())),
                "full_pin_certificates": certificates,
                "first_passing_witness": passing,
            })
        seed_status = (
            "PASSING" if any(cover["status"] == "PASSING" for cover in cover_records)
            else "INFEASIBLE"
        )
        seed_records.append({
            "seed": seed,
            "seed_mask": seed_mask,
            "status": seed_status,
            "current_covers": cover_records,
        })
    if observed_current != expected_current:
        raise ValueError("inherited current-cover count mismatch")
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
        "--output", type=Path,
        default=ROOT / "d6_k7_full_pin_increment_report.json",
    )
    parser.add_argument("--base-report-sha256", required=True)
    parser.add_argument("--base-verification-sha256", required=True)
    args = parser.parse_args()
    expected_base_hashes = {
        "d6_k7_double_pin_conjunction_report.json": args.base_report_sha256,
        "d6_k7_double_pin_conjunction_verification.json": (
            args.base_verification_sha256
        ),
    }
    for name, expected in expected_base_hashes.items():
        if sha256(ROOT / name) != expected:
            raise SystemExit(f"base artifact hash mismatch: {name}")
    base_report = json.loads(BASE_REPORT.read_text(encoding="utf-8"))
    base_verification = json.loads(BASE_VERIFICATION.read_text(encoding="utf-8"))
    if (
        base_report.get("kind")
        != "d6_k7_two_defect_double_pin_seed_conjunction"
        or base_report.get("status") != "COMPLETE"
        or base_verification.get("status") != "PASS"
        or base_verification.get("report", {}).get("sha256")
        != args.base_report_sha256
        or TARGET not in base_report["summary"]["survivor_indices"]
    ):
        raise SystemExit("verified double-pin base does not contain target survivor")
    residue = json.loads(base.RESIDUE.read_text(encoding="utf-8"))
    graph = next(
        item for item in residue["classes"]["K7"]["graphs"]
        if int(item["index"]) == TARGET
    )
    selected = {TARGET}
    prior = base.read_failure_keys(
        base.PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
    ).get(TARGET, set())
    tetrad = base.read_failure_keys(
        base.TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
    ).get(TARGET, set())
    tetrad_row = base.read_tetrad_rows(selected)[TARGET]
    started = time.monotonic()
    record = evaluate_target(
        graph, prior, tetrad, int(tetrad_row["tetrad_passing_covers"])
    )
    if record["decision"] != "REJECTED":
        raise SystemExit("full-pin target unexpectedly survived")
    base_survivors = list(base_report["summary"]["survivor_indices"])
    survivors = [index for index in base_survivors if index != TARGET]
    report = {
        "schema": 1,
        "kind": "d6_k7_full_pin_odd_cycle_increment",
        "status": "COMPLETE",
        "claim": (
            "The generalized odd-cycle/full-support pin obstruction rejects "
            "double-pin survivor 3919831, reducing the exact K7 residue from "
            "25 to 24 graphs."
        ),
        "semantics": {
            "candidate_nonedges_optional": True,
            "only_required_edges_enter_coordinate_graphs": True,
            "floating_point_enters_rejection": False,
        },
        "base_sha256": dict(sorted(expected_base_hashes.items())),
        "source_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            "d6_k7_full_pin_odd_cycle.py": sha256(
                ROOT / "d6_k7_full_pin_odd_cycle.py"
            ),
            "build_d6_k7_double_pin_conjunction.py": sha256(
                ROOT / "build_d6_k7_double_pin_conjunction.py"
            ),
        },
        "incremental_rejected_indices": [TARGET],
        "incremental_rejected_indices_sha256": stable_hash([TARGET]),
        "base_survivors": base_survivors,
        "base_survivors_sha256": stable_hash(base_survivors),
        "exact_survivors": survivors,
        "exact_survivors_sha256": stable_hash(survivors),
        "record": record,
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
    }
    atomic_json(args.output, report)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "incremental_rejected": [TARGET],
        "exact_survivors": len(survivors),
        "exact_survivors_sha256": stable_hash(survivors),
        "counts": record["counts"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
