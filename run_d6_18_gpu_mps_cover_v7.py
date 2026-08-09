#!/usr/bin/env python3
"""Heuristic all-K6-gauge MPS triage of the 181-class d=6 n=18 cover.

This program deliberately makes *no* mathematical rejection or realizability
claim.  It loads the exact production selection frozen by
``run_d6_interval_18_cover_v7.py``, enumerates every unit K6 in every selected
18-vertex graph, and applies the existing MPS/LM discovery kernel in every K6
gauge.  Candidate nonedges are never fitted to a non-unit distance; the only
all-pair term is the explicitly heuristic short-distance collision bias.

Chunks are atomic and resumable.  Random streams use a stable variant id
derived from the deletion-manifest class ordinal and K6-gauge ordinal, so they
do not depend on chunking.  The default uses MPS and one sequential CPU LM
worker, which is suitable while the interval campaign occupies the other CPU
cores.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import platform
import shlex
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Sequence

import run_d6_18_gpu_mps_screen as screen
import run_d6_interval_18_cover_v7 as cover18


ROOT = Path(__file__).resolve().parent
EXPECTED_SCREEN_SOURCE_SHA256 = (
    "cd57a6193fba1a9b05b88be566598413bd7ff79a13cdc7b5530df90a14293627"
)
EXPECTED_COVER_SOURCE_SHA256 = (
    "5eeca35b328d930c94fa33d96c714ef756f4a777fcce8c8bffccbe472d744019"
)
EXPECTED_CLASS_COUNT = 181
EXPECTED_GAUGE_COUNT = 4_099
EXPECTED_GAUGE_MANIFEST_SHA256 = (
    "571273aff8152f0d8b3afa7c03e2c230ac09c018381b86cae5b7893ce641287f"
)
EXPECTED_VARIANT_IDS_SHA256 = (
    "9d83883b7241d8db034ce43caa80a156ffe5b617052e49cba20f03dce8dce5b7"
)
EXPECTED_GAUGE_POPULATIONS = {"K7": 500, "K6": 3_599}
EXPECTED_POSITIVE_CLASS_INDEX = cover18.EXPECTED_KNOWN_POSITIVE_CLASS_INDEX
EXPECTED_POSITIVE_GAUGES = 9
VARIANT_STRIDE = 128

DEFAULT_CHECKPOINT_ROOT = ROOT / ".runs/d6_18_gpu_mps_cover_v7"
DEFAULT_REPORT = ROOT / ".runs/d6_18_gpu_mps_cover_v7_report.json"
DEFAULT_DETAILS = ROOT / ".runs/d6_18_gpu_mps_cover_v7_details.json.gz"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def dependency_boundary() -> dict:
    """Pin both imported producers before trusting their selection/semantics."""

    paths = {
        "mps_kernel": Path(screen.__file__).resolve(),
        "cover_selector": Path(cover18.__file__).resolve(),
    }
    hashes = {name: screen.file_sha256(path) for name, path in paths.items()}
    require(
        hashes["mps_kernel"] == EXPECTED_SCREEN_SOURCE_SHA256,
        "imported MPS kernel source hash drift",
    )
    require(
        hashes["cover_selector"] == EXPECTED_COVER_SOURCE_SHA256,
        "imported 181-class selector source hash drift",
    )
    return {
        name: {"path": str(paths[name]), "sha256": hashes[name]}
        for name in sorted(paths)
    }


def git_head() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def load_campaign() -> tuple[list[dict], dict, dict]:
    """Rebuild and audit the committed 181-class production selection."""

    parents, parent_provenance = cover18.load_v7_parents(
        ROOT / "d6_current_residue_manifest_v7.json",
        ROOT / "d6_current_residue_manifest_v7_verification.json",
    )
    deletions, deletion_provenance = cover18.load_deletion_boundary(
        ROOT / "d6_residue_18_deletions_v2.json",
        ROOT / "d6_residue_18_deletions_v2_verification.json",
    )
    base_graphs, dense_cover = cover18.build_dense_deletion_cover(parents, deletions)
    production_graphs, augmentation = cover18.build_actionable_campaign(
        base_graphs, dense_cover, deletions
    )
    require(len(production_graphs) == EXPECTED_CLASS_COUNT, "campaign class-count drift")
    require(
        cover18.stable_hash([graph["index"] for graph in production_graphs])
        == cover18.EXPECTED_PRODUCTION_CLASS_INDICES_SHA256,
        "campaign class-index root drift",
    )

    records = deletions["unique_deletions"]
    campaign = []
    for graph in production_graphs:
        class_index = int(graph["index"])
        raw = records[class_index]
        require(
            tuple(map(int, raw["adjacency"])) == tuple(graph["adjacency"]),
            "campaign/deletion adjacency mismatch",
        )
        campaign.append(
            {
                "class_index": class_index,
                "population": str(graph["population"]),
                "raw": raw,
            }
        )
    positives = [
        item for item in campaign if item["raw"].get("standard18_compatible")
    ]
    require(
        len(positives) == 1
        and positives[0]["class_index"] == EXPECTED_POSITIVE_CLASS_INDEX,
        "known standard-positive control drift",
    )
    provenance = {
        "parent_boundary": parent_provenance,
        "deletion_boundary": deletion_provenance,
        "production_class_count": len(campaign),
        "production_class_indices": [item["class_index"] for item in campaign],
        "production_class_indices_sha256": cover18.stable_hash(
            [item["class_index"] for item in campaign]
        ),
        "dense_cover_selected_count": dense_cover["selected_count"],
        "actionable_backup_count": augmentation["backup_count"],
        "known_positive_class_index": EXPECTED_POSITIVE_CLASS_INDEX,
        "candidate_nonedges": "unconstrained and may also have distance one",
    }
    return campaign, provenance, deletions


def gauge_metadata(campaign: Sequence[dict]) -> list[dict]:
    """Enumerate every K6 gauge in deterministic class/lexicographic order."""

    answer = []
    seen_variant_ids: set[int] = set()
    for item in campaign:
        raw = item["raw"]
        rows = screen.validate_rows(raw["adjacency"], raw["class_id"])
        seeds = screen.enumerate_cliques(rows, 6)
        require(bool(seeds), f"class {item['class_index']} has no K6 gauge")
        require(len(seeds) < VARIANT_STRIDE, "variant-id stride exhausted")
        for gauge_ordinal, seed in enumerate(seeds):
            variant_id = int(item["class_index"]) * VARIANT_STRIDE + gauge_ordinal
            require(variant_id not in seen_variant_ids, "duplicate gauge variant id")
            seen_variant_ids.add(variant_id)
            answer.append(
                {
                    "variant_id": variant_id,
                    "class_index": int(item["class_index"]),
                    "population": item["population"],
                    "gauge_ordinal": gauge_ordinal,
                    "seed": list(seed),
                }
            )
    return answer


def audit_full_gauge_manifest(metadata: Sequence[dict]) -> dict:
    require(len(metadata) == EXPECTED_GAUGE_COUNT, "all-K6 gauge-count drift")
    require(
        screen.stable_hash(list(metadata)) == EXPECTED_GAUGE_MANIFEST_SHA256,
        "all-K6 gauge manifest drift",
    )
    require(
        screen.stable_hash([item["variant_id"] for item in metadata])
        == EXPECTED_VARIANT_IDS_SHA256,
        "variant-id root drift",
    )
    populations = dict(Counter(item["population"] for item in metadata))
    require(populations == EXPECTED_GAUGE_POPULATIONS, "gauge population drift")
    positive = [
        item for item in metadata
        if item["class_index"] == EXPECTED_POSITIVE_CLASS_INDEX
    ]
    require(len(positive) == EXPECTED_POSITIVE_GAUGES, "positive gauge-count drift")
    return {
        "gauge_count": len(metadata),
        "gauge_manifest_sha256": screen.stable_hash(list(metadata)),
        "variant_ids_sha256": screen.stable_hash(
            [item["variant_id"] for item in metadata]
        ),
        "population_counts": populations,
        "known_positive_gauge_count": len(positive),
    }


def select_campaign(campaign: Sequence[dict], indices: Sequence[int]) -> list[dict]:
    if not indices:
        return list(campaign)
    requested = set(map(int, indices))
    available = {int(item["class_index"]) for item in campaign}
    require(not requested - available, f"unknown class indices: {sorted(requested - available)}")
    return [item for item in campaign if item["class_index"] in requested]


def make_variants(items: Sequence[dict]) -> tuple[list[dict], list[dict]]:
    variants = []
    metadata = gauge_metadata(items)
    by_class = {int(item["class_index"]): item for item in items}
    for info in metadata:
        item = by_class[info["class_index"]]
        variant = screen.reorder_record_with_seed(
            info["variant_id"], item["raw"], info["seed"]
        )
        require(variant["gauge_size"] == 6, "non-K6 variant entered campaign")
        variants.append(variant)
    return variants, metadata


def configuration(args: argparse.Namespace, selected: Sequence[dict], source_hash: str) -> dict:
    return {
        "schema": "d6-18-gpu-mps-cover-v7-config-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS_NO_REALIZABILITY_CONCLUSIONS",
        "orchestrator_source_sha256": source_hash,
        "dependency_source_sha256": {
            "mps_kernel": EXPECTED_SCREEN_SOURCE_SHA256,
            "cover_selector": EXPECTED_COVER_SOURCE_SHA256,
        },
        "production_class_indices_sha256": (
            cover18.EXPECTED_PRODUCTION_CLASS_INDICES_SHA256
        ),
        "selected_class_indices": [item["class_index"] for item in selected],
        "selected_class_indices_sha256": screen.stable_hash(
            [item["class_index"] for item in selected]
        ),
        "dimension": 6,
        "vertices": 18,
        "gauge": "every unit K6, lexicographically enumerated within each class",
        "variant_id_formula": "128 * deletion_manifest_class_index + gauge_ordinal",
        "total_restarts": args.total_restarts,
        "restarts_per_wave": args.restarts_per_wave,
        "optimizer": "PyTorch Adam with cosine learning-rate decay",
        "optimizer_steps": args.steps,
        "learning_rate": args.learning_rate,
        "initialization_method": "seed_sphere",
        "initialization_scale": args.initialization_scale,
        "seed_standard_controls": args.seed_standard_controls,
        "seedbase": args.seedbase,
        "graph_seed_formula": (
            "(seedbase*1000003 + variant_id*7919 + wave*104729) mod 2^64"
        ),
        "collision_distance": args.collision_distance,
        "collision_weight": args.collision_weight,
        "collision_semantics": (
            "all-pair short-distance repulsion is a heuristic distinctness bias, "
            "not a nonedge distance constraint"
        ),
        "distinctness_threshold": args.distinctness_threshold,
        "gpu_candidate_rms": args.gpu_candidate_rms,
        "cpu_refine_max_nfev": args.cpu_refine_max_nfev,
        "lm_top_k": args.lm_top_k,
        "lm_candidate_rms": args.lm_candidate_rms,
        "lm_collision_distance": args.lm_collision_distance,
        "lm_collision_weight": args.lm_collision_weight,
        "cpu_workers": args.cpu_workers,
        "class_chunk_size": args.class_chunk_size,
        "device": args.device,
        "dtype": "MPS/CPU float32 optimization and SciPy float64 LM refinement",
        "candidate_nonedges": (
            "absent from the unit-edge residual; unconstrained and may be unit"
        ),
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }


def decorate_result(result: dict, metadata: Sequence[dict]) -> dict:
    """Attach stable original-class and K6-gauge identities to kernel output."""

    lookup = {item["variant_id"]: item for item in metadata}
    require(len(lookup) == len(metadata), "duplicate metadata variant id")

    def decorate(items: Sequence[dict]) -> list[dict]:
        answer = []
        for original in items:
            item = dict(original)
            variant_id = int(item.pop("index"))
            require(variant_id in lookup, "kernel returned an unknown variant id")
            info = lookup[variant_id]
            item.update(
                {
                    "variant_id": variant_id,
                    "class_index": info["class_index"],
                    "population": info["population"],
                    "gauge_ordinal": info["gauge_ordinal"],
                    "gauge_seed_old_vertices": info["seed"],
                }
            )
            answer.append(item)
        return answer

    decorated = {
        key: decorate(result[key])
        for key in ("rows", "witnesses", "retained_endpoints", "lm_attempts")
    }
    decorated["wave_profiles"] = result["wave_profiles"]
    decorated["lm_profile"] = result["lm_profile"]
    return decorated


def checkpoint_path(run_directory: Path, start: int, end: int) -> Path:
    return run_directory / f"classes_{start:04d}_{end:04d}.json.gz"


def load_checkpoint(
    path: Path,
    config_hash: str,
    start: int,
    end: int,
    class_indices: Sequence[int],
) -> dict:
    with gzip.open(path, "rt", encoding="ascii") as stream:
        payload = json.load(stream)
    claimed = payload.pop("checkpoint_sha256", None)
    require(claimed == screen.stable_hash(payload), f"{path}: payload hash")
    payload["checkpoint_sha256"] = claimed
    require(payload.get("schema") == "d6-18-gpu-mps-cover-v7-checkpoint-v1", f"{path}: schema")
    require(payload.get("config_sha256") == config_hash, f"{path}: config")
    require(payload.get("class_range") == [start, end], f"{path}: range")
    require(payload.get("class_indices") == list(class_indices), f"{path}: classes")
    rows = payload.get("rows", [])
    metadata = payload.get("gauge_metadata", [])
    require(len(rows) == len(metadata) == payload.get("gauge_count"), f"{path}: gauges")
    require(
        [row["variant_id"] for row in rows]
        == [item["variant_id"] for item in metadata],
        f"{path}: gauge ordering",
    )
    return payload


def best_row(rows: Sequence[dict]) -> dict:
    def key(row: dict) -> tuple[float, float, int]:
        return (
            float(row["best_distinct_lm_edge_rms"]),
            float(row["best_distinct_edge_rms"]),
            int(row["gauge_ordinal"]),
        )

    return min(rows, key=key)


def summarize_classes(rows: Sequence[dict]) -> list[dict]:
    grouped: dict[int, list[dict]] = {}
    for row in rows:
        grouped.setdefault(int(row["class_index"]), []).append(row)
    answer = []
    for class_index in sorted(grouped):
        selected = sorted(grouped[class_index], key=lambda row: row["gauge_ordinal"])
        best = best_row(selected)
        lm_candidates = [row for row in selected if row["lm_candidate"]]
        gpu_candidates = [row for row in selected if row["gpu_candidate"]]
        answer.append(
            {
                "class_index": class_index,
                "class_id": selected[0]["class_id"],
                "population": selected[0]["population"],
                "edges": selected[0]["edges"],
                "standard18_compatible": bool(selected[0]["standard_compatible"]),
                "gauge_count": len(selected),
                "optimizer_starts": sum(row["optimizer_starts"] for row in selected),
                "cpu_lm_attempts": sum(row["cpu_lm_attempts"] for row in selected),
                "gpu_candidate_gauges": len(gpu_candidates),
                "lm_candidate_gauges": len(lm_candidates),
                "heuristic_status": (
                    "NUMERICAL_CANDIDATE_REQUIRES_CERTIFICATION"
                    if lm_candidates or gpu_candidates
                    else "NO_CANDIDATE_FOUND_NO_CONCLUSION"
                ),
                "best_variant_id": best["variant_id"],
                "best_gauge_ordinal": best["gauge_ordinal"],
                "best_gauge_seed_old_vertices": best["gauge_seed_old_vertices"],
                "best_distinct_lm_edge_rms": best["best_distinct_lm_edge_rms"],
                "best_distinct_lm_max_edge_error": best[
                    "best_distinct_lm_max_edge_error"
                ],
                "best_distinct_lm_minimum_distance": best[
                    "best_distinct_lm_minimum_distance"
                ],
                "mathematical_rejection": False,
                "mathematical_realizability_conclusion": False,
            }
        )
    return answer


def positive_control_audit(rows: Sequence[dict], config: dict) -> dict:
    selected = [
        row for row in rows if row["class_index"] == EXPECTED_POSITIVE_CLASS_INDEX
    ]
    if not selected:
        return {
            "included": False,
            "status": "NOT_IN_SELECTED_SUBSET",
            "mathematical_conclusion": False,
        }
    exact_initial = [
        row
        for row in selected
        if math.isfinite(row["initial_standard_control_edge_rms"])
        and row["initial_standard_control_edge_rms"] <= 2e-6
        and row["initial_standard_control_minimum_distance"]
        >= config["distinctness_threshold"]
    ]
    lm_recovered = [row for row in selected if row["lm_candidate"]]
    if not config["seed_standard_controls"]:
        status = "UNSEEDED_CALIBRATION_REPORTED_NOT_REQUIRED_TO_RECOVER"
    elif len(exact_initial) != len(selected):
        status = "FAIL_SEEDED_STANDARD_INITIALIZATION"
    elif not lm_recovered:
        status = "FAIL_STANDARD_CONTROL_NOT_RECOVERED"
    else:
        status = "PASS"
    return {
        "included": True,
        "class_index": EXPECTED_POSITIVE_CLASS_INDEX,
        "gauge_count": len(selected),
        "seeded": config["seed_standard_controls"],
        "exact_initial_gauges": len(exact_initial),
        "lm_recovered_gauges": len(lm_recovered),
        "best_lm_edge_rms": min(
            row["best_distinct_lm_edge_rms"] for row in selected
        ),
        "status": status,
        "mathematical_conclusion": False,
    }


def assemble(
    checkpoints: Sequence[dict],
    config: dict,
    config_hash: str,
    provenance: dict,
    dependencies: dict,
    args: argparse.Namespace,
    invocation_wall_seconds: float,
    run_state: dict,
) -> dict:
    keys = ("rows", "witnesses", "retained_endpoints", "lm_attempts")
    merged = {
        key: [item for checkpoint in checkpoints for item in checkpoint[key]]
        for key in keys
    }
    merged["rows"].sort(key=lambda row: (row["class_index"], row["gauge_ordinal"]))
    merged["witnesses"].sort(
        key=lambda row: (row["class_index"], row["gauge_ordinal"])
    )
    merged["retained_endpoints"].sort(
        key=lambda row: (row["class_index"], row["gauge_ordinal"])
    )
    merged["lm_attempts"].sort(
        key=lambda row: (row["class_index"], row["gauge_ordinal"])
    )
    class_summaries = summarize_classes(merged["rows"])
    control = positive_control_audit(merged["rows"], config)
    details = {
        "schema": "d6-18-gpu-mps-cover-v7-details-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS_NO_REALIZABILITY_CONCLUSIONS",
        "config_sha256": config_hash,
        "rows": merged["rows"],
        "witnesses": merged["witnesses"],
        "retained_endpoints": merged["retained_endpoints"],
        "lm_attempts": merged["lm_attempts"],
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }
    screen.atomic_gzip_json(args.details, details)
    report = {
        "schema": "d6-18-gpu-mps-cover-v7-report-v1",
        "status": (
            "COMPLETE_HEURISTIC_CONTROL_FAILURE"
            if control["status"].startswith("FAIL")
            else "COMPLETE_HEURISTIC_ONLY"
        ),
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS_NO_REALIZABILITY_CONCLUSIONS",
        "warning": (
            "Failure to find a candidate in any or all gauges proves nothing. "
            "Every candidate requires exact or rigorous interval certification."
        ),
        "config": config,
        "config_sha256": config_hash,
        "producer_source": {
            "path": str(Path(__file__).resolve()),
            "sha256": screen.file_sha256(Path(__file__).resolve()),
        },
        "dependency_boundary": dependencies,
        "git_head_at_launch": run_state["git_head_at_launch"],
        "launch_command": run_state["launch_command"],
        "assembly_command": shlex.join(sys.argv),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
        },
        "input_provenance": provenance,
        "invocation_wall_seconds": invocation_wall_seconds,
        "checkpoint_compute_wall_seconds": sum(
            checkpoint["wall_seconds"] for checkpoint in checkpoints
        ),
        "class_count": len(class_summaries),
        "gauge_count": len(merged["rows"]),
        "optimizer_starts": sum(row["optimizer_starts"] for row in merged["rows"]),
        "cpu_lm_attempts": sum(row["cpu_lm_attempts"] for row in merged["rows"]),
        "class_candidate_count": sum(
            summary["heuristic_status"].startswith("NUMERICAL_CANDIDATE")
            for summary in class_summaries
        ),
        "gauge_lm_candidate_count": sum(row["lm_candidate"] for row in merged["rows"]),
        "class_summaries": class_summaries,
        "positive_control_audit": control,
        "optional_nonedge_semantics_audit": {
            "status": "PASS_BY_PINNED_INPUT_AND_OBJECTIVE_CONSTRUCTION",
            "unit_fit_mask": "adjacency entries equal to one only",
            "candidate_nonedges": "not fitted; unconstrained and may also be unit",
            "all_pair_term": (
                "soft short-distance collision bias only; discovery aid, not a "
                "unit/non-unit constraint"
            ),
        },
        "details": {
            "path": str(args.details),
            "sha256": screen.file_sha256(args.details),
            "rows": len(merged["rows"]),
            "witnesses_with_coordinates": len(merged["witnesses"]),
            "retained_endpoints_with_coordinates": len(merged["retained_endpoints"]),
        },
        "checkpoints": {
            "count": len(checkpoints),
            "ordered_sha256": screen.stable_hash(
                [checkpoint["checkpoint_sha256"] for checkpoint in checkpoints]
            ),
        },
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }
    screen.atomic_json(args.report, report)
    return report


def run(
    selected: Sequence[dict],
    config: dict,
    config_hash: str,
    provenance: dict,
    dependencies: dict,
    args: argparse.Namespace,
) -> dict | None:
    run_directory = args.checkpoint_root / config_hash[:16]
    run_directory.mkdir(parents=True, exist_ok=True)
    state_path = run_directory / "run_state.json"
    expected_indices = [item["class_index"] for item in selected]
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        require(
            state.get("schema") == "d6-18-gpu-mps-cover-v7-run-state-v1"
            and state.get("config_sha256") == config_hash
            and state.get("class_indices") == expected_indices,
            "existing run-state boundary mismatch",
        )
    else:
        state = {
            "schema": "d6-18-gpu-mps-cover-v7-run-state-v1",
            "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
            "config_sha256": config_hash,
            "class_count": len(selected),
            "class_indices": expected_indices,
            "checkpoint_directory": str(run_directory),
            "started_at_utc": screen.utc_now(),
            "git_head_at_launch": git_head(),
            "launch_command": shlex.join(sys.argv),
            "completed_class_chunks": 0,
            "invocations": [],
            "mathematical_rejections": 0,
            "mathematical_realizability_conclusions": 0,
        }
    state.setdefault("invocations", []).append(
        {
            "started_at_utc": screen.utc_now(),
            "git_head": git_head(),
            "command": shlex.join(sys.argv),
        }
    )
    screen.atomic_json(state_path, state)
    checkpoints = []
    launched = 0
    started = time.perf_counter()
    executor = screen.make_executor(config)
    try:
        for start in range(0, len(selected), config["class_chunk_size"]):
            end = min(start + config["class_chunk_size"], len(selected))
            class_indices = [item["class_index"] for item in selected[start:end]]
            path = checkpoint_path(run_directory, start, end)
            if path.exists():
                payload = load_checkpoint(path, config_hash, start, end, class_indices)
            else:
                if args.max_new_chunks is not None and launched >= args.max_new_chunks:
                    break
                variants, metadata = make_variants(selected[start:end])
                chunk_started = time.perf_counter()
                result = screen.process_batch(variants, config, executor)
                chunk_wall_seconds = time.perf_counter() - chunk_started
                decorated = decorate_result(result, metadata)
                payload = {
                    "schema": "d6-18-gpu-mps-cover-v7-checkpoint-v1",
                    "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
                    "config_sha256": config_hash,
                    "class_range": [start, end],
                    "class_indices": class_indices,
                    "gauge_count": len(metadata),
                    "gauge_metadata": metadata,
                    "completed_at_utc": screen.utc_now(),
                    "wall_seconds": chunk_wall_seconds,
                    **decorated,
                    "mathematical_rejections": 0,
                    "mathematical_realizability_conclusions": 0,
                }
                payload["checkpoint_sha256"] = screen.stable_hash(payload)
                screen.atomic_gzip_json(path, payload)
                payload = load_checkpoint(path, config_hash, start, end, class_indices)
                launched += 1
                print(
                    f"classes {start}:{end} gauges={payload['gauge_count']} "
                    f"lm_candidates={sum(row['lm_candidate'] for row in payload['rows'])}",
                    flush=True,
                )
            checkpoints.append(payload)
            state["completed_class_chunks"] = len(checkpoints)
            state["completed_classes"] = sum(len(item["class_indices"]) for item in checkpoints)
            state["completed_gauges"] = sum(item["gauge_count"] for item in checkpoints)
            screen.atomic_json(state_path, state)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    expected_chunks = math.ceil(len(selected) / config["class_chunk_size"])
    if len(checkpoints) != expected_chunks:
        print(
            f"partial checkpoint-safe run: {len(checkpoints)}/{expected_chunks} chunks",
            flush=True,
        )
        return None
    report = assemble(
        checkpoints,
        config,
        config_hash,
        provenance,
        dependencies,
        args,
        time.perf_counter() - started,
        state,
    )
    state["status"] = report["status"]
    state["report"] = {"path": str(args.report), "sha256": screen.file_sha256(args.report)}
    screen.atomic_json(state_path, state)
    return report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--selection-only", action="store_true")
    result.add_argument(
        "--index",
        dest="indices",
        type=int,
        action="append",
        help="restrict to a production deletion-manifest class ordinal; repeatable",
    )
    result.add_argument("--device", choices=("mps", "cpu"), default="mps")
    result.add_argument("--total-restarts", type=int, default=8)
    result.add_argument("--restarts-per-wave", type=int, default=8)
    result.add_argument("--steps", type=int, default=800)
    result.add_argument("--learning-rate", type=float, default=0.03)
    result.add_argument("--initialization-scale", type=float, default=0.7)
    result.add_argument("--seedbase", type=int, default=618_181_4099)
    result.add_argument(
        "--no-seed-standard-controls",
        action="store_false",
        dest="seed_standard_controls",
    )
    result.add_argument("--collision-distance", type=float, default=0.12)
    result.add_argument("--collision-weight", type=float, default=0.05)
    result.add_argument("--distinctness-threshold", type=float, default=0.01)
    result.add_argument("--gpu-candidate-rms", type=float, default=3e-4)
    result.add_argument("--cpu-refine-max-nfev", type=int, default=1000)
    result.add_argument("--lm-top-k", type=int, default=2)
    result.add_argument("--lm-candidate-rms", type=float, default=1e-8)
    result.add_argument("--lm-collision-distance", type=float, default=0.05)
    result.add_argument("--lm-collision-weight", type=float, default=100.0)
    result.add_argument("--cpu-workers", type=int, default=1)
    result.add_argument("--class-chunk-size", type=int, default=16)
    result.add_argument(
        "--max-new-chunks",
        type=int,
        help="stop checkpoint-safely after this many newly computed chunks",
    )
    result.add_argument("--checkpoint-root", type=Path, default=DEFAULT_CHECKPOINT_ROOT)
    result.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    result.add_argument("--details", type=Path, default=DEFAULT_DETAILS)
    return result


def validate_args(args: argparse.Namespace, parser_: argparse.ArgumentParser) -> None:
    if args.total_restarts < 1 or args.restarts_per_wave < 1 or args.steps < 1:
        parser_.error("restart and optimizer step counts must be positive")
    if args.cpu_workers < 1 or args.class_chunk_size < 1:
        parser_.error("CPU worker and chunk counts must be positive")
    if args.lm_top_k < 0 or args.cpu_refine_max_nfev < 0:
        parser_.error("LM top-k and evaluation cap must be nonnegative")
    if args.max_new_chunks is not None and args.max_new_chunks < 1:
        parser_.error("max-new-chunks must be positive")
    if args.collision_distance <= 0 or args.lm_collision_distance <= 0:
        parser_.error("collision distances must be positive")


def main() -> None:
    parser_ = parser()
    args = parser_.parse_args()
    validate_args(args, parser_)
    dependencies = dependency_boundary()
    campaign, provenance, _deletions = load_campaign()
    full_metadata = gauge_metadata(campaign)
    full_audit = audit_full_gauge_manifest(full_metadata)
    selected = select_campaign(campaign, args.indices or ())
    selected_metadata = gauge_metadata(selected)
    provenance["full_gauge_manifest"] = full_audit
    provenance["selected_gauge_count"] = len(selected_metadata)
    provenance["selected_gauge_manifest_sha256"] = screen.stable_hash(selected_metadata)
    source_hash = screen.file_sha256(Path(__file__).resolve())
    config = configuration(args, selected, source_hash)
    config_hash = screen.stable_hash(config)
    if args.selection_only:
        print(
            json.dumps(
                {
                    "schema": "d6-18-gpu-mps-cover-v7-selection-v1",
                    "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
                    "config_sha256": config_hash,
                    "classes": len(selected),
                    "gauges": len(selected_metadata),
                    "optimizer_starts": len(selected_metadata) * args.total_restarts,
                    "cpu_lm_attempt_cap": len(selected_metadata) * args.lm_top_k,
                    "cpu_workers": args.cpu_workers,
                    "full_manifest_audit": full_audit,
                    "selected_class_indices": [item["class_index"] for item in selected],
                    "mathematical_rejections": 0,
                    "mathematical_realizability_conclusions": 0,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    report = run(selected, config, config_hash, provenance, dependencies, args)
    if report is not None:
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
