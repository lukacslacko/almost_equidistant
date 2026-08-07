#!/usr/bin/env python3
"""Independent structural verifier for the heuristic d=6 18-class screen.

This checker validates provenance, checkpoint coverage, raw-result integrity,
summary counts, and the explicit standard-18 positive controls.  It does not
turn the floating-point search into a proof and intentionally emits zero
mathematical rejections.
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
DEFAULT_REPORT = ROOT / "d6_18_numerical_screen_report.json"
DEFAULT_OUTPUT = ROOT / "d6_18_numerical_screen_verification.json"
DEFAULT_CORPUS = ROOT / "d6_residue_18_deletions.json"
ENGINE_SOURCE = ROOT / "d6_18_lm_screen.c"
INHERITED_SOURCE = ROOT / "lm5.c"
ORCHESTRATOR_SOURCE = ROOT / "run_d6_18_numerical_screen.py"
CHECKPOINT_ARCHIVE = (
    ROOT / "d6_18_numerical_checkpoints_3eaad0fa1ff43461.tar.gz"
)
EXPECTED_CHECKPOINT_ARCHIVE_SHA256 = (
    "3d3bfbdcd71ccbf6d654684934eb472b8f56438f02a5fecb82dd7c12dc240c61"
)
CHECKPOINT_ARCHIVE_PREFIX = (
    "d6_18_numerical_checkpoints/3eaad0fa1ff43461/"
)


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


def sign_vectors() -> list[tuple[int, ...]]:
    return [
        tuple(1 if word & (1 << axis) else -1 for axis in range(5))
        for word in range(32)
        if word.bit_count() % 2 == 1
    ]


def modular_rank(matrix: Sequence[Sequence[int]], prime: int) -> int:
    data = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(len(data[0]) if data else 0):
        pivot = next(
            (row for row in range(rank, len(data)) if data[row][column]), None
        )
        if pivot is None:
            continue
        data[rank], data[pivot] = data[pivot], data[rank]
        inverse = pow(data[rank][column], prime - 2, prime)
        data[rank] = [(entry * inverse) % prime for entry in data[rank]]
        for row in range(rank + 1, len(data)):
            factor = data[row][column]
            if factor:
                data[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(data[row], data[rank])
                ]
        rank += 1
    return rank


def control_rank(record: dict, permutation: Sequence[int]) -> int:
    coordinates = [list(point) + [0] for point in sign_vectors()]
    coordinates += [[0, 0, 0, 0, 0, 4], [0, 0, 0, 0, 0, -4]]
    embedded = [coordinates[permutation[vertex]] for vertex in range(18)]
    matrix = []
    for first in range(18):
        for second in range(first):
            if not (record["adjacency"][first] >> second) & 1:
                continue
            row = [0] * 108
            for axis in range(6):
                difference = embedded[first][axis] - embedded[second][axis]
                row[first * 6 + axis] = difference
                row[second * 6 + axis] = -difference
            matrix.append(row)
    return modular_rank(matrix, 13)


def parse_raw(path: Path, columns: Sequence[str]) -> list[dict]:
    with gzip.open(path, "rt", encoding="ascii") as stream:
        lines = stream.read().splitlines()
    if not lines or lines[0].split("\t") != list(columns):
        raise AssertionError("raw TSV header disagrees with the report")
    integer_columns = {
        "index", "n", "edges", "has_seed", "seed_numerical_rigidity_rank",
        "random_endpoint_numerical_rigidity_rank", "random_near_solution_count",
        "random_distinct_near_solution_count",
    }
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(columns):
            raise AssertionError("malformed raw TSV row")
        row = {
            key: int(value) if key in integer_columns else float(value)
            for key, value in zip(columns, fields)
        }
        rows.append(row)
    return rows


def verify(report_path: Path, corpus_path: Path) -> dict:
    report_hash = sha256_file(report_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    corpus = json.loads(corpus_path.read_text(encoding="utf-8"))
    checks = {}
    checks["report_schema"] = report["schema"] == "d6-18-numerical-screen-report-v1"
    checks["heuristic_status"] = (
        report["proof_status"] == "HEURISTIC_ONLY_NO_REJECTIONS"
        and report["summary"]["mathematical_rejections"] == 0
        and report["summary"]["mathematical_realizability_conclusions"] == 0
    )
    checks["corpus_hash"] = sha256_file(corpus_path) == report["corpus"]["sha256"]
    checks["engine_source_hash"] = (
        sha256_file(ENGINE_SOURCE) == report["engine"]["source_sha256"]
    )
    checks["orchestrator_source_hash"] = (
        sha256_file(ORCHESTRATOR_SOURCE)
        == report["engine"]["orchestrator_sha256"]
        == report["config"]["orchestrator_source_sha256"]
    )
    checks["inherited_source_hash"] = (
        sha256_file(INHERITED_SOURCE)
        == report["engine"]["inherited_source_sha256"]
    )
    checks["config_hash"] = stable_hash(report["config"]) == report["config_sha256"]

    reported_results_path = Path(report["results"]["path"])
    results_path = ROOT / reported_results_path.name
    checks["portable_results_path"] = (
        reported_results_path.name == "d6_18_numerical_screen_results.tsv.gz"
    )
    checks["results_hash"] = sha256_file(results_path) == report["results"]["sha256"]
    rows = parse_raw(results_path, report["results"]["columns"])
    checks["complete_ordered_rows"] = (
        len(rows) == len(corpus["unique_deletions"]) == 12_712
        and [row["index"] for row in rows] == list(range(12_712))
    )
    checks["edge_counts_match_corpus"] = all(
        row["edges"] == record["edges"] and row["n"] == 18
        for row, record in zip(rows, corpus["unique_deletions"])
    )
    expected_controls = [
        index for index, record in enumerate(corpus["unique_deletions"])
        if record["standard18_compatible"]
    ]
    checks["seed_flags_match_standard_controls"] = [
        row["index"] for row in rows if row["has_seed"]
    ] == expected_controls
    controls_by_index = {
        item["index"]: item for item in report["summary"]["control_results"]
    }
    controls_ok = set(controls_by_index) == set(expected_controls)
    for index in expected_controls:
        item = controls_by_index[index]
        record = corpus["unique_deletions"][index]
        permutation = record["standard18_pole_pair_witnesses"][0][
            "embedding_permutation"
        ]
        rank = control_rank(record, permutation)
        controls_ok &= item["mod13_rigidity_rank_lower_bound"] == rank
        controls_ok &= item["seeded_positive_control_pass"]
        controls_ok &= rows[index]["seed_residual"] < 1e-20
        controls_ok &= rows[index]["seed_minimum_distance"] > 1e-3
        controls_ok &= (
            rows[index]["seed_numerical_rigidity_rank"] == rank
        )
    checks["all_standard_controls_pass_and_ranks_recomputed"] = bool(controls_ok)

    near_indices = [
        row["index"] for row in rows
        if row["random_best_distinct_residual"] < 1e-20
        and row["random_best_distinct_minimum_distance"] > 1e-3
    ]
    checks["near_solution_summary_recomputed"] = near_indices == report[
        "summary"
    ]["random_distinct_near_solution_indices"]
    witnesses = report["summary"]["random_distinct_near_solution_witnesses"]
    witnesses_ok = [item["index"] for item in witnesses] == near_indices
    for item in witnesses:
        points = item["coordinates_vertex_major_binary64"]
        record = corpus["unique_deletions"][item["index"]]
        residuals = []
        minimum = math.inf
        for first in range(18):
            for second in range(first):
                squared = sum(
                    (points[first][axis] - points[second][axis]) ** 2
                    for axis in range(6)
                )
                minimum = min(minimum, math.sqrt(squared))
                if record["adjacency"][first] & (1 << second):
                    residuals.append(squared - 1.0)
        witnesses_ok &= abs(
            sum(value * value for value in residuals)
            - item["residual_sum_of_squares_recomputed"]
        ) < 1e-26
        witnesses_ok &= max(map(abs, residuals)) < 1e-12
        witnesses_ok &= abs(
            minimum - item["minimum_pairwise_distance_recomputed"]
        ) < 1e-12
    checks["candidate_coordinates_recomputed"] = bool(witnesses_ok)
    threshold_counts = {
        format(threshold, ".0e"): sum(
            row["random_best_distinct_residual"] < threshold for row in rows
        )
        for threshold in (1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e-2, 1e-1)
    }
    checks["threshold_counts_recomputed"] = threshold_counts == report[
        "summary"
    ]["random_residual_threshold_counts"]

    # The committed gzip tar archive is the immutable checkpoint boundary.
    # This avoids relying on the producer's absolute local checkpoint path and
    # keeps 401 small files out of Git while retaining every chunk, stderr
    # record, graph input, binary, and the atomic run-state manifest.
    checks["checkpoint_archive_hash"] = (
        sha256_file(CHECKPOINT_ARCHIVE)
        == EXPECTED_CHECKPOINT_ARCHIVE_SHA256
    )
    with tarfile.open(CHECKPOINT_ARCHIVE, "r:gz") as archive:
        def archived_bytes(relative: str) -> bytes:
            member = archive.extractfile(CHECKPOINT_ARCHIVE_PREFIX + relative)
            if member is None:
                raise AssertionError(f"missing checkpoint member: {relative}")
            return member.read()

        state_bytes = archived_bytes("run_state.json")
        checks["run_state_hash"] = (
            hashlib.sha256(state_bytes).hexdigest()
            == report["checkpoint"]["run_state_sha256"]
        )
        state = json.loads(state_bytes.decode("utf-8"))
        expected_start = 0
        chunks_ok = state["status"] == "COMPLETE"
        for chunk in state["completed_chunks"]:
            chunks_ok &= chunk["start"] == expected_start
            expected_start = chunk["end"]
            chunks_ok &= (
                hashlib.sha256(archived_bytes(chunk["path"])).hexdigest()
                == chunk["sha256"]
            )
            chunks_ok &= (
                hashlib.sha256(archived_bytes(chunk["stderr_path"])).hexdigest()
                == chunk["stderr_sha256"]
            )
        chunks_ok &= expected_start == 12_712
        checks["checkpoint_chunks_complete_and_hashed"] = bool(chunks_ok)
        checks["checkpoint_binary_hash"] = (
            hashlib.sha256(archived_bytes("d6_18_lm_screen")).hexdigest()
            == report["engine"]["binary_sha256"]
        )
        checks["checkpoint_graph_input_hash"] = (
            hashlib.sha256(archived_bytes("graphs_with_controls.txt")).hexdigest()
            == report["config"]["graph_input_sha256"]
        )

    failed = sorted(key for key, passed in checks.items() if not passed)
    return {
        "schema": "d6-18-numerical-screen-verification-v1",
        "status": "PASS" if not failed else "FAIL",
        "scope": (
            "structural/provenance verification of a heuristic screen; no "
            "non-realizability or realizability conclusion"
        ),
        "report_path": str(report_path.resolve()),
        "report_sha256": report_hash,
        "results_path": str(results_path.resolve()),
        "results_sha256": sha256_file(results_path),
        "graphs_recomputed_structurally": len(rows),
        "positive_controls_recomputed": len(expected_controls),
        "mathematical_rejections": 0,
        "checks": checks,
        "failed_checks": failed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = verify(args.report, args.corpus)
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
