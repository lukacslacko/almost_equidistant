#!/usr/bin/env python3
"""Independent integrity checker for the heuristic d=6 MPS screen.

The checker validates the immutable inputs, producer source, result table,
atomic checkpoints, summary arithmetic, exact standard-18 positive controls,
and every retained numerical witness.  It deliberately makes no rejection or
realizability inference from optimizer failure or floating-point residuals.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import tarfile
import tempfile
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "d6_18_gpu_mps_screen_report.json"
DEFAULT_CORPUS = ROOT / "d6_residue_18_deletions.json"
DEFAULT_OUTPUT = ROOT / "d6_18_gpu_mps_screen_verification.json"
PRODUCER = ROOT / "run_d6_18_gpu_mps_screen.py"
CHECKPOINT_ARCHIVE = (
    ROOT / "d6_18_gpu_mps_checkpoints_3249d4d3a4afd60d.tar.gz"
)
EXPECTED_CHECKPOINT_ARCHIVE_SHA256 = (
    "9836f25e908c175f96203c24921051087c56b210e0812312b3f95b47dbb82df9"
)
CHECKPOINT_ARCHIVE_PREFIX = "3249d4d3a4afd60d/"
EXPECTED_CORPUS_SHA256 = (
    "9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732"
)
INTEGER_COLUMNS = {
    "index",
    "edges",
    "standard_compatible",
    "gauge_size",
    "optimizer_starts",
    "gpu_candidate",
    "cpu_lm_attempts",
    "cpu_lm_successful_terminations",
    "best_distinct_lm_input_rank",
    "first_lm_candidate_input_rank",
    "lm_candidate",
}


def sha256_file(path: Path) -> str:
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
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def parse_results(path: Path, columns: Sequence[str]) -> list[dict]:
    with gzip.open(path, "rt", encoding="ascii") as stream:
        lines = stream.read().splitlines()
    if not lines or lines[0].split("\t") != list(columns):
        raise AssertionError("result columns do not match report")
    rows = []
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) != len(columns):
            raise AssertionError("malformed result row")
        rows.append(
            {
                key: int(value) if key in INTEGER_COLUMNS else float(value)
                for key, value in zip(columns, values)
            }
        )
    return rows


def same_float(left: float, right: float) -> bool:
    return (math.isnan(left) and math.isnan(right)) or left == right


def sign_vectors() -> list[tuple[int, ...]]:
    return [
        tuple(1 if word & (1 << axis) else -1 for axis in range(5))
        for word in range(32)
        if word.bit_count() & 1
    ]


def standard18_scaled() -> list[tuple[int, ...]]:
    # Divide these coordinates by sqrt(8).  Every required unit edge then has
    # squared integer distance 8; the two poles use height +/-sqrt(3), so we
    # check them symbolically below instead of introducing a decimal sqrt.
    return [tuple(vector) + (0,) for vector in sign_vectors()]


def exact_standard_control(record: dict) -> bool:
    permutation = record["standard18_pole_pair_witnesses"][0][
        "embedding_permutation"
    ]
    base = standard18_scaled()
    if sorted(permutation) != list(range(18)):
        return False
    for first in range(18):
        for second in range(first):
            if not (record["adjacency"][first] >> second) & 1:
                continue
            image_first = permutation[first]
            image_second = permutation[second]
            first_pole = image_first >= 16
            second_pole = image_second >= 16
            if not first_pole and not second_pole:
                squared_times_eight = sum(
                    (base[image_first][axis] - base[image_second][axis]) ** 2
                    for axis in range(5)
                )
                if squared_times_eight != 8:
                    return False
            elif first_pole and second_pole:
                # Pole-pole distance squared is 3/2, hence it cannot be a unit
                # edge in a valid embedded control.
                return False
            else:
                # A halfcube vertex has norm squared 5/8 and a pole has norm
                # squared 3/8, with zero cross term: distance squared is one.
                pass
    return True


def witness_metrics(points: Sequence[Sequence[float]], record: dict) -> tuple[float, float, float]:
    errors = []
    minimum_squared = math.inf
    for first in range(18):
        for second in range(first):
            squared = sum(
                (points[first][axis] - points[second][axis]) ** 2
                for axis in range(6)
            )
            minimum_squared = min(minimum_squared, squared)
            if (record["adjacency"][first] >> second) & 1:
                errors.append(squared - 1.0)
    return (
        math.sqrt(sum(value * value for value in errors) / len(errors)),
        max(map(abs, errors)),
        math.sqrt(max(minimum_squared, 0.0)),
    )


def close_metric(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2e-6, abs_tol=2e-8)


def close_mps_metric(left: float, right: float) -> bool:
    """Compare a binary64 replay with a metric originally reduced on MPS.

    Retained MPS endpoints are serialized exactly as their binary32
    coordinates, but their archived metrics were reduced by Torch in
    binary32 and a device-dependent summation order.  Recomputing the same
    expression from those coordinates in Python binary64 therefore need not
    reproduce the final few binary32 ulps.  This envelope is used only for
    the heuristic MPS endpoints; LM witnesses and summaries retain the
    tighter binary64 comparison above.
    """

    return math.isclose(left, right, rel_tol=1e-5, abs_tol=1e-6)


def first_clique(rows: Sequence[int], target: int) -> tuple[int, ...] | None:
    def search(chosen: tuple[int, ...], candidates: int) -> tuple[int, ...] | None:
        if len(chosen) == target:
            return chosen
        if len(chosen) + candidates.bit_count() < target:
            return None
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            answer = search(chosen + (vertex,), candidates & rows[vertex])
            if answer is not None:
                return answer
        return None

    return search((), (1 << 18) - 1)


def deterministic_gauge(record: dict) -> tuple[int, tuple[int, ...]]:
    rows = record["adjacency"]
    seed = first_clique(rows, 7)
    if seed is not None:
        return 7, seed
    seed = first_clique(rows, 6)
    if seed is None:
        raise AssertionError("corpus graph has no K6 gauge")
    return 6, seed


def verify(report_path: Path, corpus_path: Path) -> dict:
    report_hash = sha256_file(report_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    corpus = json.loads(corpus_path.read_text(encoding="utf-8"))
    records = corpus["unique_deletions"]
    checks: dict[str, bool] = {}

    checks["report_schema"] = report.get("schema") == "d6-18-gpu-mps-screen-report-v1"
    checks["heuristic_only_zero_conclusions"] = (
        report.get("proof_status") == "HEURISTIC_ONLY_NO_REJECTIONS"
        and report["summary"].get("mathematical_rejections") == 0
        and report["summary"].get("mathematical_realizability_conclusions") == 0
        and report["config"].get("semantics")
        == "HEURISTIC_ONLY; no numerical rejection or realizability conclusion"
    )
    checks["pinned_corpus"] = (
        sha256_file(corpus_path)
        == report["corpus"]["sha256"]
        == EXPECTED_CORPUS_SHA256
        and len(records) == 12_712
    )
    checks["producer_hash"] = (
        sha256_file(PRODUCER)
        == report["source"]["sha256"]
        == report["config"]["orchestrator_source_sha256"]
    )
    checks["config_hash"] = stable_hash(report["config"]) == report["config_sha256"]

    results_path = ROOT / Path(report["results"]["path"]).name
    checks["results_hash"] = sha256_file(results_path) == report["results"]["sha256"]
    rows = parse_results(results_path, report["results"]["columns"])
    gauges = [deterministic_gauge(record) for record in records]
    checks["complete_ordered_rows"] = (
        len(rows) == len(records) == report["results"]["rows"] == 12_712
        and [row["index"] for row in rows] == list(range(12_712))
    )
    checks["row_graph_metadata"] = all(
        row["edges"] == record["edges"]
        and row["standard_compatible"] == int(record["standard18_compatible"])
        and row["gauge_size"] == gauges[index][0]
        and row["optimizer_starts"] == report["config"]["total_restarts"]
        for index, (row, record) in enumerate(zip(rows, records))
    )

    control_indices = [
        index for index, record in enumerate(records)
        if record["standard18_compatible"]
    ]
    control_set = set(control_indices)
    if report["config"]["seed_standard_controls"]:
        control_metrics_ok = all(
            rows[index]["initial_standard_control_edge_rms"] < 2e-6
            and rows[index]["initial_standard_control_minimum_distance"]
            > report["config"]["distinctness_threshold"]
            for index in control_indices
        )
    else:
        control_metrics_ok = all(
            math.isnan(rows[index]["initial_standard_control_edge_rms"])
            and math.isnan(rows[index]["initial_standard_control_minimum_distance"])
            for index in control_indices
        )
    checks["fourteen_exact_standard_controls"] = (
        len(control_indices) == 14
        and all(exact_standard_control(records[index]) for index in control_indices)
        and control_metrics_ok
        and all(
            math.isnan(row["initial_standard_control_edge_rms"])
            and math.isnan(row["initial_standard_control_minimum_distance"])
            for index, row in enumerate(rows) if index not in control_set
        )
    )

    checkpoint_rows = []
    checkpoint_endpoints = []
    checkpoint_lm_attempts = []
    checkpoint_witnesses = []
    checkpoint_hashes = []
    checkpoint_ok = (
        report["config_sha256"][:16] == CHECKPOINT_ARCHIVE_PREFIX.rstrip("/")
        and sha256_file(CHECKPOINT_ARCHIVE)
        == EXPECTED_CHECKPOINT_ARCHIVE_SHA256
    )
    chunk_size = report["config"]["chunk_size"]
    chunk_names = [
        f"chunk_{start:05d}_{min(12_712, start + chunk_size):05d}.json.gz"
        for start in range(0, 12_712, chunk_size)
    ]
    expected_members = {
        CHECKPOINT_ARCHIVE_PREFIX.rstrip("/"),
        CHECKPOINT_ARCHIVE_PREFIX + "run_state.json",
        *(CHECKPOINT_ARCHIVE_PREFIX + name for name in chunk_names),
    }
    with tarfile.open(CHECKPOINT_ARCHIVE, "r:gz") as archive:
        checkpoint_ok &= set(archive.getnames()) == expected_members
        run_state_member = archive.extractfile(
            CHECKPOINT_ARCHIVE_PREFIX + "run_state.json"
        )
        checkpoint_ok &= run_state_member is not None
        if run_state_member is not None:
            run_state = json.loads(run_state_member.read().decode("utf-8"))
            checkpoint_ok &= (
                run_state.get("schema") == "d6-18-gpu-mps-run-state-v1"
                and run_state.get("proof_status")
                == "HEURISTIC_ONLY_NO_REJECTIONS"
                and run_state.get("config_sha256") == report["config_sha256"]
                and stable_hash(run_state.get("config"))
                == report["config_sha256"]
            )
        for start, name in zip(range(0, 12_712, chunk_size), chunk_names):
            end = min(12_712, start + chunk_size)
            member = archive.extractfile(CHECKPOINT_ARCHIVE_PREFIX + name)
            if member is None:
                checkpoint_ok = False
                continue
            payload = json.loads(gzip.decompress(member.read()).decode("ascii"))
            claimed = payload.pop("checkpoint_sha256", None)
            checkpoint_ok &= (
                claimed is not None
                and stable_hash(payload) == claimed
                and payload.get("schema") == "d6-18-gpu-mps-checkpoint-v1"
                and payload.get("proof_status")
                == "HEURISTIC_ONLY_NO_REJECTIONS"
                and payload.get("config_sha256") == report["config_sha256"]
                and payload.get("range") == [start, end]
                and payload.get("mathematical_rejections") == 0
                and [row.get("index") for row in payload.get("rows", [])]
                == list(range(start, end))
                and [item.get("index") for item in payload.get("lm_attempts", [])]
                == list(range(start, end))
            )
            for row in payload.get("rows", []):
                index = row["index"]
                checkpoint_ok &= (
                    row["gauge_size"] == gauges[index][0]
                    and tuple(row["gauge_seed_old_vertices"])
                    == gauges[index][1]
                )
            checkpoint_hashes.append(claimed)
            checkpoint_rows.extend(payload["rows"])
            checkpoint_endpoints.extend(payload.get("retained_endpoints", []))
            checkpoint_lm_attempts.extend(payload.get("lm_attempts", []))
            checkpoint_witnesses.extend(payload.get("witnesses", []))
    checks["immutable_checkpoint_archive_rehashed"] = (
        checkpoint_ok
        and len(checkpoint_hashes) == report["checkpoints"]["count"]
        and stable_hash(checkpoint_hashes) == report["checkpoints"]["ordered_sha256"]
    )
    checkpoint_table_ok = len(checkpoint_rows) == len(rows)
    for raw, compact in zip(checkpoint_rows, rows):
        for column in report["results"]["columns"]:
            value = raw[column]
            if column in INTEGER_COLUMNS:
                checkpoint_table_ok &= value == compact[column]
            else:
                checkpoint_table_ok &= same_float(float(value), compact[column])
    checks["checkpoints_equal_result_table"] = checkpoint_table_ok

    endpoints_ok = (
        [item["index"] for item in checkpoint_endpoints]
        == sorted(item["index"] for item in checkpoint_endpoints)
        and len({item["index"] for item in checkpoint_endpoints})
        == len(checkpoint_endpoints)
    )
    for item in checkpoint_endpoints:
        index = item["index"]
        metrics = witness_metrics(
            item["coordinates_original_vertex_order"], records[index]
        )
        expected = (
            rows[index]["best_distinct_edge_rms"],
            rows[index]["best_distinct_max_edge_error"],
            rows[index]["best_distinct_minimum_distance"],
        )
        endpoints_ok &= all(
            close_mps_metric(left, right)
            for left, right in zip(metrics, expected)
        )
        endpoints_ok &= all(
            close_mps_metric(left, right)
            for left, right in zip(
                metrics,
                (item["edge_rms"], item["max_edge_error"], item["minimum_distance"]),
            )
        )
    checks["retained_endpoint_metrics_recomputed"] = bool(endpoints_ok)

    lm_attempts_ok = len(checkpoint_lm_attempts) == len(rows)
    for item, row in zip(checkpoint_lm_attempts, rows):
        attempts = item["attempts"]
        lm_attempts_ok &= item["index"] == row["index"]
        lm_attempts_ok &= len(attempts) == row["cpu_lm_attempts"]
        lm_attempts_ok &= (
            sum(attempt["successful_termination"] for attempt in attempts)
            == row["cpu_lm_successful_terminations"]
        )
        best = min(attempts, key=lambda attempt: attempt["edge_rms"], default=None)
        distinct = [
            attempt for attempt in attempts
            if attempt["minimum_distance"] >= report["config"]["distinctness_threshold"]
        ]
        best_distinct = min(
            distinct, key=lambda attempt: attempt["edge_rms"], default=None
        )
        lm_attempts_ok &= close_metric(
            best["edge_rms"] if best else math.inf,
            row["best_lm_edge_rms_including_collisions"],
        )
        lm_attempts_ok &= close_metric(
            best["max_edge_error"] if best else math.inf,
            row["best_lm_max_edge_error_including_collisions"],
        )
        lm_attempts_ok &= close_metric(
            best["minimum_distance"] if best else 0.0,
            row["best_lm_minimum_distance_including_collisions"],
        )
        lm_attempts_ok &= close_metric(
            best_distinct["edge_rms"] if best_distinct else math.inf,
            row["best_distinct_lm_edge_rms"],
        )
        lm_attempts_ok &= close_metric(
            best_distinct["max_edge_error"] if best_distinct else math.inf,
            row["best_distinct_lm_max_edge_error"],
        )
        lm_attempts_ok &= close_metric(
            best_distinct["minimum_distance"] if best_distinct else 0.0,
            row["best_distinct_lm_minimum_distance"],
        )
        lm_attempts_ok &= (
            (best_distinct["input_rank"] if best_distinct else -1)
            == row["best_distinct_lm_input_rank"]
        )
        candidate_ranks = [
            attempt["input_rank"] for attempt in distinct
            if attempt["edge_rms"] <= report["config"]["lm_candidate_rms"]
        ]
        lm_attempts_ok &= min(candidate_ranks, default=-1) == row[
            "first_lm_candidate_input_rank"
        ]
        lm_attempts_ok &= bool(candidate_ranks) == bool(row["lm_candidate"])
    checks["lm_attempt_summaries_recomputed"] = bool(lm_attempts_ok)

    gpu_threshold_counts = {
        format(threshold, ".0e"): sum(
            row["best_distinct_edge_rms"] <= threshold for row in rows
        )
        for threshold in (1e-2, 1e-3, 3e-4, 1e-4, 1e-5)
    }
    lm_threshold_counts = {
        format(threshold, ".0e"): sum(
            row["best_distinct_lm_edge_rms"] <= threshold for row in rows
        )
        for threshold in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
    }
    gpu_candidates = [row["index"] for row in rows if row["gpu_candidate"]]
    lm_candidates = [row["index"] for row in rows if row["lm_candidate"]]
    seeded = 14 if report["config"]["seed_standard_controls"] else 0
    checks["summary_recomputed"] = (
        gpu_threshold_counts == report["summary"]["gpu_distinct_rms_threshold_counts"]
        and lm_threshold_counts == report["summary"]["lm_distinct_rms_threshold_counts"]
        and gpu_candidates == report["summary"]["gpu_candidate_indices"]
        and lm_candidates == report["summary"]["lm_candidate_indices"]
        and report["summary"]["standard_compatible"] == 14
        and report["summary"]["nonstandard_compatible"] == 12_698
        and report["summary"]["gauge_distribution"]
        == {
            "K6": sum(gauge[0] == 6 for gauge in gauges),
            "K7": sum(gauge[0] == 7 for gauge in gauges),
        }
        and report["summary"]["total_optimizer_starts"]
        == 12_712 * report["config"]["total_restarts"]
        and report["summary"]["seeded_exact_control_starts"] == seeded
        and report["summary"]["random_starts"]
        == 12_712 * report["config"]["total_restarts"] - seeded
        and report["summary"]["retained_best_distinct_endpoint_count"]
        == len(checkpoint_endpoints)
        and report["summary"]["cpu_lm_attempts"]
        == sum(row["cpu_lm_attempts"] for row in rows)
    )

    witnesses = report["summary"]["candidate_witnesses"]
    witness_indices = sorted(set(gpu_candidates) | set(lm_candidates))
    witnesses_ok = (
        witnesses == checkpoint_witnesses
        and [item["index"] for item in witnesses] == witness_indices
    )
    for item in witnesses:
        index = item["index"]
        metrics = witness_metrics(
            item["coordinates_original_vertex_order"], records[index]
        )
        witnesses_ok &= all(
            close_metric(left, right)
            for left, right in zip(
                metrics,
                (item["edge_rms"], item["max_edge_error"], item["minimum_distance"]),
            )
        )
        witnesses_ok &= metrics[2] >= report["config"]["distinctness_threshold"]
        if item["source"] == "CPU_float64_analytic_LM":
            witnesses_ok &= rows[index]["lm_candidate"] == 1
            expected = (
                rows[index]["best_distinct_lm_edge_rms"],
                rows[index]["best_distinct_lm_max_edge_error"],
                rows[index]["best_distinct_lm_minimum_distance"],
            )
            witnesses_ok &= item["lm_input_rank"] == rows[index][
                "best_distinct_lm_input_rank"
            ]
        else:
            witnesses_ok &= item["source"] == "MPS_float32"
            witnesses_ok &= item["lm_input_rank"] is None
            expected = (
                rows[index]["best_distinct_edge_rms"],
                rows[index]["best_distinct_max_edge_error"],
                rows[index]["best_distinct_minimum_distance"],
            )
        witnesses_ok &= all(close_metric(left, right) for left, right in zip(metrics, expected))
    checks["candidate_witness_metrics_recomputed"] = bool(witnesses_ok)

    failed = sorted(key for key, passed in checks.items() if not passed)
    return {
        "schema": "d6-18-gpu-mps-screen-verification-v1",
        "status": "PASS" if not failed else "FAIL",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
        "report": report_path.name,
        "report_sha256": report_hash,
        "verifier_source_sha256": sha256_file(Path(__file__).resolve()),
        "checkpoint_archive": {
            "path": CHECKPOINT_ARCHIVE.name,
            "sha256": EXPECTED_CHECKPOINT_ARCHIVE_SHA256,
            "members": len(expected_members),
        },
        "checks": checks,
        "failed_checks": failed,
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = verify(args.report, args.corpus)
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
