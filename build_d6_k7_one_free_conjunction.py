#!/usr/bin/env python3
"""Build the exact K7 one-free-edge conjunction on the frozen 24 residue."""

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
from datetime import UTC, datetime
from pathlib import Path

import build_d6_k7_double_pin_conjunction as base
import d6_k7_correlated_one_free_edge as correlated
import d6_k7_full_pin_odd_cycle as full_pin
import d6_k7_one_free_edge as singleton
import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation
import d6_k7_two_defect_double_pin as double_pin
import run_d6_k7_support_capacity_pilot as inherited


ROOT = Path(__file__).resolve().parent
BASE_REPORT = ROOT / "d6_k7_double_pin_conjunction_report.json"
BASE_VERIFICATION = ROOT / "d6_k7_double_pin_conjunction_verification.json"
FULL_PIN_REPORT = ROOT / "d6_k7_full_pin_increment_report.json"
FULL_PIN_VERIFICATION = ROOT / "d6_k7_full_pin_increment_verification.json"
PIPELINES = ("old", "singleton", "correlated", "combined")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("ascii")
    ).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def evaluate_graph(payload: tuple) -> dict:
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
    observed_current_keys = set()

    for seed_mask in reference.clique_masks(adj, 7):
        totals["seeds"] += 1
        seed, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        cover_records = []
        for zmask in reference.eligible_covers(ladj, eligible):
            totals["eligible_covers"] += 1
            baseline = reference.analyze_cover(
                adj,
                outside,
                defects,
                zmask,
                support_solver,
                zero_forcing,
                total_term_rank,
                clique_solver,
            )
            key = (tuple(seed), zmask)
            status = inherited.current_cover_status(
                adj,
                outside,
                defects,
                zmask,
                baseline,
                key,
                prior,
                tetrad,
            )
            totals[f"inherited_cover_{status}"] += 1
            if status != "passing":
                continue
            observed_current_keys.add(key)
            totals["current_covers"] += 1
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex
                for vertex in range(len(outside))
                if not zmask & (1 << vertex)
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nvertices]
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            counts: Counter = Counter()
            first_witness = {pipeline: None for pipeline in PIPELINES}
            singleton_certificates = []
            correlated_certificates = []

            for z_supports in propagation.labeled_support_families(z_allowed):
                counts["labeled_z_families"] += 1
                propagated = propagation.analyze_support_assignment(
                    graph_n,
                    z_supports,
                    n_allowed,
                    zero_forcing,
                    clique_solver,
                )
                if propagated.failure is not None:
                    counts[f"propagation_failure:{propagated.failure}"] += 1
                    continue
                masks = propagated.propagated_masks
                sparse = sparse_value.check_small_support_masks(graph_n, masks)
                if not sparse.feasible:
                    counts["sparse_value_failure"] += 1
                    continue
                counts["pre_pinning_families"] += 1

                old_double = double_pin.find_double_pin(graph_n, masks)
                old_full = None
                if old_double is not None:
                    double_pin.verify_certificate(graph_n, masks, old_double)
                    counts["old_double_pin_failure"] += 1
                else:
                    old_full = full_pin.find_full_pin(graph_n, masks)
                    if old_full is not None:
                        full_pin.verify_certificate(graph_n, masks, old_full)
                        counts["old_full_pin_failure"] += 1
                old_failed = old_double is not None or old_full is not None
                if old_failed:
                    continue
                counts["pre_new_families"] += 1

                singleton_pin = singleton.find_one_free_edge(graph_n, masks)
                correlated_pin = correlated.find_correlated_one_free_edge(
                    graph_n, masks
                )
                if singleton_pin is not None:
                    singleton.verify_certificate(graph_n, masks, singleton_pin)
                    counts["singleton_failure"] += 1
                    singleton_certificates.append(
                        {
                            "z_supports": list(z_supports),
                            "propagated_masks": list(masks),
                            "certificate": singleton.certificate_json(
                                singleton_pin
                            ),
                        }
                    )
                if correlated_pin is not None:
                    correlated.verify_certificate(graph_n, masks, correlated_pin)
                    counts["correlated_failure"] += 1
                    correlated_certificates.append(
                        {
                            "z_supports": list(z_supports),
                            "propagated_masks": list(masks),
                            "certificate": correlated.certificate_json(
                                correlated_pin
                            ),
                        }
                    )
                witness = {
                    "z_supports": list(z_supports),
                    "propagated_masks": list(masks),
                }
                survives = {
                    "old": True,
                    "singleton": singleton_pin is None,
                    "correlated": correlated_pin is None,
                    "combined": (
                        singleton_pin is None and correlated_pin is None
                    ),
                }
                for pipeline, passing in survives.items():
                    if passing:
                        counts[f"{pipeline}_passing_families"] += 1
                        if first_witness[pipeline] is None:
                            first_witness[pipeline] = witness
                    else:
                        counts[f"{pipeline}_infeasible_families"] += 1

            status_by_pipeline = {
                pipeline: (
                    "PASSING"
                    if first_witness[pipeline] is not None
                    else "INFEASIBLE"
                )
                for pipeline in PIPELINES
            }
            for pipeline, pipeline_status in status_by_pipeline.items():
                totals[
                    f"{pipeline}_{pipeline_status.lower()}_covers"
                ] += 1
            totals.update(counts)
            cover_records.append(
                {
                    "seed": seed,
                    "zmask": zmask,
                    "status_by_pipeline": status_by_pipeline,
                    "family_counts": dict(sorted(counts.items())),
                    "singleton_certificates": singleton_certificates,
                    "correlated_certificates": correlated_certificates,
                    "first_passing_witness_by_pipeline": first_witness,
                }
            )

        if not cover_records:
            raise ValueError("frozen survivor has no inherited-passing cover")
        status_by_pipeline = {
            pipeline: (
                "PASSING"
                if any(
                    cover["status_by_pipeline"][pipeline] == "PASSING"
                    for cover in cover_records
                )
                else "INFEASIBLE"
            )
            for pipeline in PIPELINES
        }
        seed_records.append(
            {
                "seed": seed,
                "seed_mask": seed_mask,
                "status_by_pipeline": status_by_pipeline,
                "current_covers": cover_records,
            }
        )

    if len(observed_current_keys) != int(expected_current):
        raise ValueError("inherited current-cover reconstruction mismatch")
    decision_by_pipeline = {}
    first_rejecting_seed_mask_by_pipeline = {}
    for pipeline in PIPELINES:
        rejecting = next(
            (
                seed_record
                for seed_record in seed_records
                if seed_record["status_by_pipeline"][pipeline] == "INFEASIBLE"
            ),
            None,
        )
        decision_by_pipeline[pipeline] = (
            "REJECTED" if rejecting is not None else "SURVIVOR"
        )
        first_rejecting_seed_mask_by_pipeline[pipeline] = (
            int(rejecting["seed_mask"]) if rejecting is not None else 0
        )
    return {
        "index": index,
        "decision_by_pipeline": decision_by_pipeline,
        "first_rejecting_seed_mask_by_pipeline": (
            first_rejecting_seed_mask_by_pipeline
        ),
        "counts": dict(sorted(totals.items())),
        "seeds": seed_records,
    }


def validate_upstream(hashes: dict[str, str]) -> tuple[dict, dict, dict, dict]:
    paths = {
        BASE_REPORT.name: BASE_REPORT,
        BASE_VERIFICATION.name: BASE_VERIFICATION,
        FULL_PIN_REPORT.name: FULL_PIN_REPORT,
        FULL_PIN_VERIFICATION.name: FULL_PIN_VERIFICATION,
    }
    for name, path in paths.items():
        if sha256(path) != hashes[name]:
            raise ValueError(f"upstream artifact hash mismatch: {name}")
    base_report = json.loads(BASE_REPORT.read_text(encoding="utf-8"))
    base_verification = json.loads(BASE_VERIFICATION.read_text(encoding="utf-8"))
    full_report = json.loads(FULL_PIN_REPORT.read_text(encoding="utf-8"))
    full_verification = json.loads(
        FULL_PIN_VERIFICATION.read_text(encoding="utf-8")
    )
    if (
        base_report.get("kind")
        != "d6_k7_two_defect_double_pin_seed_conjunction"
        or base_report.get("status") != "COMPLETE"
        or base_verification.get("status") != "PASS"
        or base_verification.get("report", {}).get("sha256")
        != hashes[BASE_REPORT.name]
        or full_report.get("kind") != "d6_k7_full_pin_odd_cycle_increment"
        or full_report.get("status") != "COMPLETE"
        or full_verification.get("status") != "PASS"
        or full_verification.get("report", {}).get("sha256")
        != hashes[FULL_PIN_REPORT.name]
        or full_report.get("base_sha256")
        != {
            BASE_REPORT.name: hashes[BASE_REPORT.name],
            BASE_VERIFICATION.name: hashes[BASE_VERIFICATION.name],
        }
        or full_verification.get("base_sha256")
        != {
            BASE_REPORT.name: hashes[BASE_REPORT.name],
            BASE_VERIFICATION.name: hashes[BASE_VERIFICATION.name],
        }
    ):
        raise ValueError("upstream semantic or hash-binding gate failed")
    return base_report, base_verification, full_report, full_verification


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-report-sha256", required=True)
    parser.add_argument("--base-verification-sha256", required=True)
    parser.add_argument("--full-pin-report-sha256", required=True)
    parser.add_argument("--full-pin-verification-sha256", required=True)
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_one_free_conjunction_report.json",
    )
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    artifact_hashes = {
        BASE_REPORT.name: args.base_report_sha256,
        BASE_VERIFICATION.name: args.base_verification_sha256,
        FULL_PIN_REPORT.name: args.full_pin_report_sha256,
        FULL_PIN_VERIFICATION.name: args.full_pin_verification_sha256,
    }
    _, _, full_report, _ = validate_upstream(artifact_hashes)
    for name, expected in {
        **base.EXPECTED_HASHES,
        **base.EXPECTED_SOURCE_HASHES,
    }.items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"frozen root hash mismatch: {name}")

    residue = json.loads(base.RESIDUE.read_text(encoding="utf-8"))
    all_graphs = residue["classes"]["K7"]["graphs"]
    selected_indices = list(map(int, full_report["exact_survivors"]))
    if (
        len(selected_indices) != 24
        or stable_hash(selected_indices)
        != full_report["exact_survivors_sha256"]
    ):
        raise ValueError("upstream full-pin residue is not the frozen 24-list")
    graph_by_index = {int(graph["index"]): graph for graph in all_graphs}
    if any(index not in graph_by_index for index in selected_indices):
        raise ValueError("upstream survivor missing from v5 K7 corpus")
    graphs = [graph_by_index[index] for index in selected_indices]
    selected = set(selected_indices)
    prior = base.read_failure_keys(
        base.PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
    )
    tetrad = base.read_failure_keys(
        base.TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
    )
    tetrad_rows = base.read_tetrad_rows(selected)
    payloads = [
        (
            graph,
            [
                [list(seed), zmask]
                for seed, zmask in sorted(prior.get(index, set()))
            ],
            [
                [list(seed), zmask]
                for seed, zmask in sorted(tetrad.get(index, set()))
            ],
            int(tetrad_rows[index]["tetrad_passing_covers"]),
        )
        for graph, index in zip(graphs, selected_indices, strict=True)
    ]

    started_utc = datetime.now(UTC).isoformat()
    started = time.monotonic()
    if args.workers == 1:
        records = list(map(evaluate_graph, payloads))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            records = list(executor.map(evaluate_graph, payloads, chunksize=1))
    elapsed = time.monotonic() - started
    if [record["index"] for record in records] != selected_indices:
        raise ValueError("worker output order changed")

    rejected_by_pipeline = {
        pipeline: [
            record["index"]
            for record in records
            if record["decision_by_pipeline"][pipeline] == "REJECTED"
        ]
        for pipeline in PIPELINES
    }
    survivors_by_pipeline = {
        pipeline: [
            record["index"]
            for record in records
            if record["decision_by_pipeline"][pipeline] == "SURVIVOR"
        ]
        for pipeline in PIPELINES
    }
    if rejected_by_pipeline["old"] or survivors_by_pipeline["old"] != selected_indices:
        raise ValueError("old pinning replay disagrees with frozen 24 residue")
    simple = set(rejected_by_pipeline["singleton"])
    corr = set(rejected_by_pipeline["correlated"])
    combined = set(rejected_by_pipeline["combined"])
    if not simple <= combined or not corr <= combined:
        raise ValueError("combined conjunction lost an individual rejection")
    totals: Counter = Counter()
    for record in records:
        totals.update(record["counts"])

    source_names = (
        Path(__file__).name,
        "d6_k7_one_free_edge.py",
        "d6_k7_correlated_one_free_edge.py",
        "d6_k7_two_defect_double_pin.py",
        "d6_k7_full_pin_odd_cycle.py",
        "d6_k7_support_propagation.py",
        "d6_k7_small_support_value.py",
        "d6_k7_rank_reference.py",
        "run_d6_k7_support_capacity_pilot.py",
    )
    report = {
        "schema": 1,
        "kind": "d6_k7_one_free_edge_seed_conjunction",
        "status": "COMPLETE",
        "claim": (
            "The singleton and correlated pinned-product one-free-edge "
            "obstructions were conjoined with the frozen exact pinning "
            "layers over every seed, inherited cover, and propagated family "
            "of the frozen 24-graph K7 residue."
        ),
        "semantics": {
            "candidate_nonedges_optional": True,
            "only_required_edges_enter_rejection": True,
            "propagated_masks_are_support_supersets": True,
            "floating_point_enters_rejection": False,
            "graph_rejected_if_any_k7_seed_is_infeasible": True,
        },
        "upstream_artifact_sha256": dict(sorted(artifact_hashes.items())),
        "source_sha256": {
            name: sha256(ROOT / name) for name in source_names
        },
        "input_indices": selected_indices,
        "input_indices_sha256": stable_hash(selected_indices),
        "summary": {
            "input_graphs": len(selected_indices),
            "rejected_by_pipeline": rejected_by_pipeline,
            "rejected_by_pipeline_sha256": {
                pipeline: stable_hash(indices)
                for pipeline, indices in rejected_by_pipeline.items()
            },
            "survivors_by_pipeline": survivors_by_pipeline,
            "survivors_by_pipeline_sha256": {
                pipeline: stable_hash(indices)
                for pipeline, indices in survivors_by_pipeline.items()
            },
            "singleton_correlated_graph_overlap": sorted(simple & corr),
            "combined_only_synergy_rejections": sorted(
                combined - simple - corr
            ),
            "totals": dict(sorted(totals.items())),
        },
        "records": records,
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "workers": args.workers,
            "started_utc": started_utc,
            "finished_utc": datetime.now(UTC).isoformat(),
            "elapsed_seconds": elapsed,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "git": base.git_provenance(),
        },
    }
    atomic_json(args.output, report)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output),
                "input_graphs": len(selected_indices),
                "rejected_by_pipeline": rejected_by_pipeline,
                "exact_survivors": survivors_by_pipeline["combined"],
                "exact_survivors_sha256": stable_hash(
                    survivors_by_pipeline["combined"]
                ),
                "elapsed_seconds": elapsed,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
