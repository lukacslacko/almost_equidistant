#!/usr/bin/env python3
"""Reproducible 11-worker numerical triage of the d=6 18-deletion corpus.

This is a discovery tool, not a proof engine.  Required unit edges are fitted
by deterministic random-start Levenberg--Marquardt.  Candidate nonedges are
never constrained.  A large residual proves nothing, and even a tiny residual
is only a candidate pending exact/interval reconstruction and distinctness
certification.

Work is divided into small half-open index chunks.  Each successful chunk is
validated and atomically renamed, so an interrupted invocation can be resumed
with the identical configuration.  Known standard-18 compatible classes are
also evaluated at their explicit exact-coordinate embedding as seeded positive
controls; seeded results are kept separate from random-start discoveries.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parent
DEFAULT_CORPUS = ROOT / "d6_residue_18_deletions.json"
DEFAULT_CORPUS_VERIFICATION = ROOT / "d6_residue_18_deletions_verification.json"
ENGINE_SOURCE = ROOT / "d6_18_lm_screen.c"
INHERITED_ENGINE_SOURCE = ROOT / "lm5.c"
DEFAULT_CHECKPOINT_ROOT = ROOT / "d6_18_numerical_checkpoints"
DEFAULT_REPORT = ROOT / "d6_18_numerical_screen_report.json"
DEFAULT_RESULTS = ROOT / "d6_18_numerical_screen_results.tsv.gz"

RESULT_COLUMNS = (
    "index",
    "n",
    "edges",
    "has_seed",
    "seed_residual",
    "seed_minimum_distance",
    "seed_numerical_rigidity_rank",
    "random_best_residual_including_collisions",
    "random_best_distinct_residual",
    "random_best_distinct_minimum_distance",
    "random_endpoint_numerical_rigidity_rank",
    "random_near_solution_count",
    "random_distinct_near_solution_count",
)
NEAR_RESIDUAL = 1e-20
DISTINCT_TOLERANCE = 1e-3
RANK_RELATIVE_TOLERANCE = 1e-10


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


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


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def atomic_gzip_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
                stream.write(text.encode("ascii"))
            raw.flush()
            os.fsync(raw.fileno())
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


def standard_coordinates() -> list[list[float]]:
    scale = 2.0 * math.sqrt(2.0)
    points = [[coordinate / scale for coordinate in point] + [0.0]
              for point in sign_vectors()]
    height = math.sqrt(3.0) / scale
    return points + [[0.0] * 5 + [height], [0.0] * 5 + [-height]]


def modular_rank(matrix: Sequence[Sequence[int]], prime: int) -> int:
    data = [[entry % prime for entry in row] for row in matrix]
    if not data:
        return 0
    rank = 0
    for column in range(len(data[0])):
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
        if rank == len(data):
            break
    return rank


def control_mod13_rank(record: dict, permutation: Sequence[int]) -> int:
    # Standard coordinates scaled into Q(sqrt(3)), reduced modulo 13 with
    # sqrt(3)=4.  A nonzero modular minor is an exact characteristic-zero
    # lower bound.  It is not labeled an exact rank unless a matching upper
    # bound has separately been supplied.
    coordinates = [list(point) + [0] for point in sign_vectors()]
    coordinates += [[0, 0, 0, 0, 0, 4], [0, 0, 0, 0, 0, -4]]
    embedded = [coordinates[permutation[vertex]] for vertex in range(18)]
    rows = []
    for first in range(18):
        for second in range(first):
            if not (record["adjacency"][first] >> second) & 1:
                continue
            row = [0] * 108
            for axis in range(6):
                difference = embedded[first][axis] - embedded[second][axis]
                row[first * 6 + axis] = difference
                row[second * 6 + axis] = -difference
            rows.append(row)
    return modular_rank(rows, 13)


def graph_input(corpus: dict) -> tuple[str, dict[int, dict]]:
    coordinates = standard_coordinates()
    lines = []
    controls = {}
    for index, record in enumerate(corpus["unique_deletions"]):
        fields = ["18", *(str(mask) for mask in record["adjacency"])]
        if record["standard18_compatible"]:
            witness = record["standard18_pole_pair_witnesses"][0]
            permutation = witness["embedding_permutation"]
            seed = [coordinate
                    for vertex in range(18)
                    for coordinate in coordinates[permutation[vertex]]]
            fields += ["1", *(format(value, ".17g") for value in seed)]
            controls[index] = {
                "class_id": record["class_id"],
                "embedding_permutation": permutation,
                "pole_pair": witness["pole_pair"],
                "mod13_rigidity_rank_lower_bound": control_mod13_rank(
                    record, permutation
                ),
            }
        else:
            fields.append("0")
        lines.append(" ".join(fields))
    return "\n".join(lines) + "\n", controls


def parse_result_line(line: str) -> dict:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != len(RESULT_COLUMNS):
        raise ValueError(f"expected {len(RESULT_COLUMNS)} fields, got {len(fields)}")
    integers = {0, 1, 2, 3, 6, 10, 11, 12}
    values = [int(value) if position in integers else float(value)
              for position, value in enumerate(fields)]
    return dict(zip(RESULT_COLUMNS, values))


def validate_chunk(path: Path, start: int, end: int) -> list[dict]:
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) != end - start:
        raise ValueError(
            f"{path}: expected {end - start} rows, found {len(lines)}"
        )
    parsed = [parse_result_line(line) for line in lines]
    if [row["index"] for row in parsed] != list(range(start, end)):
        raise ValueError(f"{path}: noncontiguous or wrong indices")
    if any(row["n"] != 18 for row in parsed):
        raise ValueError(f"{path}: unexpected graph order")
    return parsed


def compile_engine(
    compiler: str, flags: Sequence[str], binary: Path
) -> tuple[list[str], str]:
    binary.parent.mkdir(parents=True, exist_ok=True)
    temporary = binary.with_name(binary.name + f".tmp.{os.getpid()}")
    command = [compiler, *flags, str(ENGINE_SOURCE), "-lm", "-o", str(temporary)]
    completed = subprocess.run(
        command, cwd=ROOT, text=True, capture_output=True, check=False
    )
    if completed.returncode:
        raise RuntimeError(
            f"compiler exited {completed.returncode}:\n{completed.stdout}{completed.stderr}"
        )
    os.replace(temporary, binary)
    return command, sha256_file(binary)


def run_chunk(
    binary: Path,
    graph_path: Path,
    checkpoint: Path,
    stderr_path: Path,
    start: int,
    end: int,
    restarts: int,
    seedbase: int,
) -> dict:
    if checkpoint.exists():
        validate_chunk(checkpoint, start, end)
        return {
            "start": start,
            "end": end,
            "path": checkpoint.name,
            "sha256": sha256_file(checkpoint),
            "stderr_path": stderr_path.name if stderr_path.exists() else None,
            "stderr_sha256": sha256_file(stderr_path) if stderr_path.exists() else None,
            "resumed": True,
        }
    temporary = checkpoint.with_name(checkpoint.name + f".tmp.{os.getpid()}")
    temporary_stderr = stderr_path.with_name(
        stderr_path.name + f".tmp.{os.getpid()}"
    )
    command = [
        str(binary), str(graph_path), str(start), str(end), str(restarts),
        str(seedbase),
    ]
    with temporary.open("w", encoding="ascii") as output, \
            temporary_stderr.open("w", encoding="ascii") as errors:
        completed = subprocess.run(
            command, cwd=ROOT, stdout=output, stderr=errors, check=False
        )
        output.flush()
        os.fsync(output.fileno())
        errors.flush()
        os.fsync(errors.fileno())
    if completed.returncode:
        raise RuntimeError(
            f"chunk [{start},{end}) exited {completed.returncode}; "
            f"diagnostics retained at {temporary_stderr}"
        )
    validate_chunk(temporary, start, end)
    os.replace(temporary, checkpoint)
    os.replace(temporary_stderr, stderr_path)
    return {
        "start": start,
        "end": end,
        "path": checkpoint.name,
        "sha256": sha256_file(checkpoint),
        "stderr_path": stderr_path.name,
        "stderr_sha256": sha256_file(stderr_path),
        "resumed": False,
    }


def quantile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def numerical_affine_rank(points: Sequence[Sequence[float]], tolerance: float = 1e-9) -> int:
    matrix = [
        [point[axis] - points[0][axis] for axis in range(6)]
        for point in points[1:]
    ]
    rank = 0
    scale = max(abs(value) for row in matrix for value in row)
    while rank < len(matrix) and rank < 6:
        pivot_row, pivot_column, pivot_size = -1, -1, 0.0
        for row in range(rank, len(matrix)):
            for column in range(rank, 6):
                if abs(matrix[row][column]) > pivot_size:
                    pivot_row, pivot_column = row, column
                    pivot_size = abs(matrix[row][column])
        if pivot_size <= tolerance * scale:
            break
        matrix[rank], matrix[pivot_row] = matrix[pivot_row], matrix[rank]
        for row in matrix:
            row[rank], row[pivot_column] = row[pivot_column], row[rank]
        pivot = matrix[rank][rank]
        for row in range(rank + 1, len(matrix)):
            factor = matrix[row][rank] / pivot
            for column in range(rank, 6):
                matrix[row][column] -= factor * matrix[rank][column]
        rank += 1
    return rank


def collect_candidate_witnesses(
    run_directory: Path,
    ranges: Sequence[tuple[int, int]],
    corpus: dict,
) -> list[dict]:
    raw = {}
    for start, end in ranges:
        path = run_directory / f"chunk_{start:05d}_{end:05d}.stderr"
        for line in path.read_text(encoding="ascii").splitlines():
            if not line:
                continue
            fields = line.split("\t")
            if fields[0] != "HEURISTIC-WITNESS" or len(fields) != 5 + 108:
                raise RuntimeError(f"unrecognized engine diagnostic in {path}: {line[:80]}")
            index = int(fields[1])
            if index in raw:
                raise RuntimeError(f"duplicate numerical witness for index {index}")
            raw[index] = {
                "residual": float(fields[2]),
                "minimum_distance": float(fields[3]),
                "numerical_rigidity_rank": int(fields[4]),
                "flat_coordinates": [float(value) for value in fields[5:]],
            }
    standard = standard_coordinates()
    answer = []
    for index in sorted(raw):
        record = corpus["unique_deletions"][index]
        flat = raw[index]["flat_coordinates"]
        points = [flat[6 * vertex : 6 * vertex + 6] for vertex in range(18)]
        edge_residuals = []
        squared_distances = []
        unit_pairs = []
        extra_unit_pairs = []
        minimum = math.inf
        for first in range(18):
            for second in range(first):
                squared = sum(
                    (points[first][axis] - points[second][axis]) ** 2
                    for axis in range(6)
                )
                squared_distances.append(squared)
                minimum = min(minimum, math.sqrt(squared))
                required = bool(record["adjacency"][first] & (1 << second))
                if required:
                    edge_residuals.append(squared - 1.0)
                if abs(squared - 1.0) < 1e-8:
                    unit_pairs.append([second, first])
                    if not required:
                        extra_unit_pairs.append([second, first])
        distance_histogram = Counter(round(value, 8) for value in squared_distances)
        standard_distance_discrepancy = None
        if record["standard18_compatible"]:
            permutation = record["standard18_pole_pair_witnesses"][0][
                "embedding_permutation"
            ]
            embedded = [standard[permutation[vertex]] for vertex in range(18)]
            discrepancies = []
            for first in range(18):
                for second in range(first):
                    observed = sum(
                        (points[first][axis] - points[second][axis]) ** 2
                        for axis in range(6)
                    )
                    expected = sum(
                        (embedded[first][axis] - embedded[second][axis]) ** 2
                        for axis in range(6)
                    )
                    discrepancies.append(abs(observed - expected))
            standard_distance_discrepancy = max(discrepancies)
        answer.append({
            "index": index,
            "class_id": record["class_id"],
            "standard18_compatible": record["standard18_compatible"],
            "heuristic_only": True,
            "residual_sum_of_squares_recomputed": sum(
                value * value for value in edge_residuals
            ),
            "maximum_required_edge_squared_distance_error": max(
                map(abs, edge_residuals)
            ),
            "minimum_pairwise_distance_recomputed": minimum,
            "numerical_affine_rank": numerical_affine_rank(points),
            "numerical_rigidity_rank": raw[index]["numerical_rigidity_rank"],
            "unit_pairs_at_squared_distance_tolerance_1e-8": len(unit_pairs),
            "extra_unit_pairs_at_squared_distance_tolerance_1e-8": extra_unit_pairs,
            "squared_distance_histogram_rounded_8_decimals": {
                format(value, ".8g"): count
                for value, count in sorted(distance_histogram.items())
            },
            "maximum_squared_distance_discrepancy_from_recorded_standard_embedding": (
                standard_distance_discrepancy
            ),
            "coordinates_vertex_major_binary64": points,
        })
    return answer


def summarize(
    rows: Sequence[dict], corpus: dict, controls: dict[int, dict],
    candidate_witnesses: Sequence[dict],
) -> dict:
    thresholds = [1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e-2, 1e-1]
    finite_distinct = [
        row["random_best_distinct_residual"] for row in rows
        if math.isfinite(row["random_best_distinct_residual"])
    ]
    control_results = []
    for index in sorted(controls):
        row = rows[index]
        metadata = controls[index]
        control_results.append({
            "index": index,
            **metadata,
            "seed_residual": row["seed_residual"],
            "seed_minimum_distance": row["seed_minimum_distance"],
            "seed_numerical_rigidity_rank": row[
                "seed_numerical_rigidity_rank"
            ],
            "random_best_distinct_residual": row[
                "random_best_distinct_residual"
            ],
            "random_best_distinct_minimum_distance": row[
                "random_best_distinct_minimum_distance"
            ],
            "random_endpoint_numerical_rigidity_rank": row[
                "random_endpoint_numerical_rigidity_rank"
            ],
            "random_found_near_solution": row[
                "random_best_distinct_residual"
            ] < NEAR_RESIDUAL,
            "seeded_positive_control_pass": (
                row["seed_residual"] < NEAR_RESIDUAL
                and row["seed_minimum_distance"] > DISTINCT_TOLERANCE
            ),
        })
    near = [
        row for row in rows
        if row["random_best_distinct_residual"] < NEAR_RESIDUAL
        and row["random_best_distinct_minimum_distance"] > DISTINCT_TOLERANCE
    ]
    lowest = sorted(
        rows, key=lambda row: row["random_best_distinct_residual"]
    )[:100]
    lowest_records = []
    for row in lowest:
        record = corpus["unique_deletions"][row["index"]]
        lowest_records.append({
            "index": row["index"],
            "class_id": record["class_id"],
            "edges": row["edges"],
            "contains_K7": record["contains_K7"],
            "standard18_compatible": record["standard18_compatible"],
            "random_best_distinct_residual": row[
                "random_best_distinct_residual"
            ],
            "random_best_distinct_minimum_distance": row[
                "random_best_distinct_minimum_distance"
            ],
            "random_endpoint_numerical_rigidity_rank": row[
                "random_endpoint_numerical_rigidity_rank"
            ],
        })
    by_stratum: dict[str, list[float]] = defaultdict(list)
    for row, record in zip(rows, corpus["unique_deletions"]):
        stratum = "contains_K7" if record["contains_K7"] else "K6_without_K7"
        by_stratum[stratum].append(row["random_best_distinct_residual"])
    return {
        "interpretation": (
            "HEURISTIC ONLY. No graph is rejected by this screen. A small "
            "residual is a reconstruction target, while a large residual is "
            "only failure of these starts and this optimizer. Numerical "
            "Jacobian rank is triage, not an exact rigidity certificate."
        ),
        "graphs_screened": len(rows),
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
        "graphs_still_unresolved_by_proof_after_this_screen": len(rows),
        "seeded_positive_controls": len(control_results),
        "seeded_positive_controls_passed": sum(
            item["seeded_positive_control_pass"] for item in control_results
        ),
        "control_results": control_results,
        "random_distinct_near_solution_indices": [row["index"] for row in near],
        "random_distinct_near_solution_count": len(near),
        "random_distinct_near_solution_noncontrol_indices": [
            row["index"] for row in near if row["index"] not in controls
        ],
        "random_distinct_near_solution_witnesses": list(candidate_witnesses),
        "random_residual_threshold_counts": {
            format(threshold, ".0e"): sum(
                value < threshold for value in finite_distinct
            ) for threshold in thresholds
        },
        "random_best_distinct_residual_quantiles": {
            "minimum": min(finite_distinct),
            "q01": quantile(finite_distinct, 0.01),
            "q10": quantile(finite_distinct, 0.10),
            "median": quantile(finite_distinct, 0.50),
            "q90": quantile(finite_distinct, 0.90),
            "q99": quantile(finite_distinct, 0.99),
            "maximum": max(finite_distinct),
        },
        "random_endpoint_rank_histogram_all_endpoints": {
            str(rank): count for rank, count in sorted(Counter(
                row["random_endpoint_numerical_rigidity_rank"] for row in rows
            ).items())
        },
        "strata": {
            stratum: {
                "count": len(values),
                "minimum_random_best_distinct_residual": min(values),
                "median_random_best_distinct_residual": quantile(values, 0.5),
                "below_1e-8": sum(value < 1e-8 for value in values),
                "below_1e-20": sum(value < 1e-20 for value in values),
            }
            for stratum, values in sorted(by_stratum.items())
        },
        "lowest_100_random_residuals": lowest_records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument(
        "--corpus-verification", type=Path, default=DEFAULT_CORPUS_VERIFICATION
    )
    parser.add_argument("--checkpoint-root", type=Path, default=DEFAULT_CHECKPOINT_ROOT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--restarts", type=int, default=6)
    parser.add_argument("--seedbase", type=int, default=600180001)
    args = parser.parse_args()
    if args.workers < 1 or args.chunk_size < 1 or args.restarts < 0:
        parser.error("workers and chunk-size must be positive; restarts nonnegative")

    corpus_hash = sha256_file(args.corpus)
    verification_hash = sha256_file(args.corpus_verification)
    corpus = json.loads(args.corpus.read_text(encoding="utf-8"))
    verification = json.loads(args.corpus_verification.read_text(encoding="utf-8"))
    if corpus.get("schema") != "d6-residue-18-deletion-manifest-v1":
        raise RuntimeError("unexpected corpus schema")
    if verification.get("status") != "PASS":
        raise RuntimeError("the deletion corpus verification is not PASS")
    if len(corpus["unique_deletions"]) != 12_712:
        raise RuntimeError("unexpected deletion corpus size")

    compiler = shutil.which("clang")
    if compiler is None:
        raise RuntimeError("clang not found")
    compiler_version = subprocess.run(
        [compiler, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    flags = ["-O3", "-std=c11", "-Wall", "-Wextra", "-Werror"]
    graph_text, controls = graph_input(corpus)
    graph_hash = hashlib.sha256(graph_text.encode("ascii")).hexdigest()
    config = {
        "schema": "d6-18-numerical-screen-config-v1",
        "corpus_sha256": corpus_hash,
        "corpus_verification_sha256": verification_hash,
        "orchestrator_source_sha256": sha256_file(Path(__file__).resolve()),
        "engine_source_sha256": sha256_file(ENGINE_SOURCE),
        "inherited_lm_source_sha256": sha256_file(INHERITED_ENGINE_SOURCE),
        "compiler_path": compiler,
        "compiler_version": compiler_version,
        "compiler_flags": flags,
        "dimension": 6,
        "graph_input_sha256": graph_hash,
        "graph_count": len(corpus["unique_deletions"]),
        "standard18_seeded_control_indices": sorted(controls),
        "workers": args.workers,
        "chunk_size": args.chunk_size,
        "random_restarts_per_graph": args.restarts,
        "seedbase": args.seedbase,
        "seed_formula": (
            "seedbase*1000003 + corpus_index*7919 + restart_index + 1, "
            "computed modulo uint64"
        ),
        "lm_iterations_per_restart": 400,
        "near_residual_threshold": NEAR_RESIDUAL,
        "distinctness_screen_threshold": DISTINCT_TOLERANCE,
        "numerical_rank_relative_pivot_threshold": RANK_RELATIVE_TOLERANCE,
        "semantics": "heuristic numerical triage only; zero mathematical rejections",
    }
    config_hash = stable_hash(config)
    run_directory = args.checkpoint_root / config_hash[:16]
    run_directory.mkdir(parents=True, exist_ok=True)
    graph_path = run_directory / "graphs_with_controls.txt"
    if graph_path.exists():
        if sha256_file(graph_path) != graph_hash:
            raise RuntimeError(f"resume graph input hash mismatch: {graph_path}")
    else:
        atomic_text(graph_path, graph_text)
    binary = run_directory / "d6_18_lm_screen"
    compile_command, binary_hash = compile_engine(compiler, flags, binary)

    state_path = run_directory / "run_state.json"
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state["config_sha256"] != config_hash:
            raise RuntimeError("checkpoint run-state configuration mismatch")
    else:
        state = {
            "schema": "d6-18-numerical-screen-run-state-v1",
            "status": "RUNNING",
            "started_utc": utc_now(),
            "config": config,
            "config_sha256": config_hash,
            "binary_sha256": binary_hash,
            "compile_command": compile_command,
            "completed_chunks": [],
        }
        atomic_json(state_path, state)

    count = len(corpus["unique_deletions"])
    ranges = [
        (start, min(start + args.chunk_size, count))
        for start in range(0, count, args.chunk_size)
    ]
    print(
        f"screening {count} graphs: workers={args.workers} "
        f"restarts={args.restarts} chunks={len(ranges)} run={run_directory.name}",
        flush=True,
    )
    started_monotonic = time.monotonic()
    chunks_by_start = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {}
        for start, end in ranges:
            checkpoint = run_directory / f"chunk_{start:05d}_{end:05d}.tsv"
            errors = run_directory / f"chunk_{start:05d}_{end:05d}.stderr"
            future = pool.submit(
                run_chunk, binary, graph_path, checkpoint, errors,
                start, end, args.restarts, args.seedbase,
            )
            futures[future] = (start, end)
        finished = 0
        for future in concurrent.futures.as_completed(futures):
            record = future.result()
            chunks_by_start[record["start"]] = record
            finished += 1
            state["completed_chunks"] = [
                chunks_by_start[start] for start in sorted(chunks_by_start)
            ]
            state["last_checkpoint_utc"] = utc_now()
            atomic_json(state_path, state)
            if finished == 1 or finished % 20 == 0 or finished == len(ranges):
                completed_graphs = sum(
                    item["end"] - item["start"]
                    for item in chunks_by_start.values()
                )
                elapsed = time.monotonic() - started_monotonic
                print(
                    f"checkpoint {finished}/{len(ranges)}: "
                    f"{completed_graphs}/{count} graphs; {elapsed:.1f}s",
                    flush=True,
                )

    rows = []
    raw_lines = ["\t".join(RESULT_COLUMNS)]
    for start, end in ranges:
        path = run_directory / f"chunk_{start:05d}_{end:05d}.tsv"
        chunk_rows = validate_chunk(path, start, end)
        rows.extend(chunk_rows)
        raw_lines.extend(path.read_text(encoding="ascii").splitlines())
    if [row["index"] for row in rows] != list(range(count)):
        raise AssertionError("merged result indices are incomplete")
    raw_text = "\n".join(raw_lines) + "\n"
    atomic_gzip_text(args.results, raw_text)
    candidate_witnesses = collect_candidate_witnesses(
        run_directory, ranges, corpus
    )
    near_indices = [
        row["index"] for row in rows
        if row["random_best_distinct_residual"] < NEAR_RESIDUAL
        and row["random_best_distinct_minimum_distance"] > DISTINCT_TOLERANCE
    ]
    if [item["index"] for item in candidate_witnesses] != near_indices:
        raise RuntimeError("candidate coordinate diagnostics do not match result rows")

    elapsed = time.monotonic() - started_monotonic
    state.update({
        "status": "COMPLETE",
        "completed_utc": utc_now(),
        "elapsed_this_invocation_seconds": elapsed,
        "completed_chunks": [chunks_by_start[start] for start in sorted(chunks_by_start)],
        "results_path": str(args.results.resolve()),
        "results_sha256": sha256_file(args.results),
    })
    atomic_json(state_path, state)
    report = {
        "schema": "d6-18-numerical-screen-report-v1",
        "status": "COMPLETE",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
        "created_utc": utc_now(),
        "config": config,
        "config_sha256": config_hash,
        "corpus": {
            "path": str(args.corpus.resolve()),
            "sha256": corpus_hash,
            "verification_path": str(args.corpus_verification.resolve()),
            "verification_sha256": verification_hash,
            "verification_status": verification["status"],
        },
        "engine": {
            "orchestrator_path": str(Path(__file__).resolve()),
            "orchestrator_sha256": config["orchestrator_source_sha256"],
            "source_path": str(ENGINE_SOURCE.resolve()),
            "source_sha256": config["engine_source_sha256"],
            "inherited_source_path": str(INHERITED_ENGINE_SOURCE.resolve()),
            "inherited_source_sha256": config["inherited_lm_source_sha256"],
            "binary_sha256": binary_hash,
            "compile_command": compile_command,
            "compiler_version": compiler_version,
            "floating_point_assumptions": (
                "ordinary IEEE-754 binary64 arithmetic and platform libm; "
                "results are not certificates and have no ulp trust claim"
            ),
        },
        "machine": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "logical_cpu_count_seen": os.cpu_count(),
            "workers": args.workers,
            "gpu_used": False,
            "gpu_note": (
                "Each task is a small branch-heavy 108-variable LM solve; "
                "11 independent CPU processes give the practical parallelism."
            ),
        },
        "checkpoint": {
            "directory": str(run_directory.resolve()),
            "run_state_path": str(state_path.resolve()),
            "run_state_sha256": sha256_file(state_path),
            "chunks": len(ranges),
            "chunk_size": args.chunk_size,
            "resume_command": " ".join(
                [
                    sys.executable, str(Path(__file__).resolve()),
                    "--workers", str(args.workers),
                    "--chunk-size", str(args.chunk_size),
                    "--restarts", str(args.restarts),
                    "--seedbase", str(args.seedbase),
                ]
            ),
        },
        "results": {
            "path": str(args.results.resolve()),
            "sha256": sha256_file(args.results),
            "format": "deterministic gzip TSV",
            "columns": list(RESULT_COLUMNS),
        },
        "timing": {
            "elapsed_this_invocation_seconds": elapsed,
            "graphs_per_second_this_invocation": count / elapsed,
        },
        "summary": summarize(rows, corpus, controls, candidate_witnesses),
    }
    atomic_json(args.report, report)
    print(
        f"COMPLETE report={args.report} sha256={sha256_file(args.report)} "
        f"results_sha256={report['results']['sha256']}",
        flush=True,
    )
    print(
        json.dumps({
            key: report["summary"][key]
            for key in (
                "graphs_screened", "seeded_positive_controls_passed",
                "random_distinct_near_solution_count",
                "random_distinct_near_solution_noncontrol_indices",
            )
        }, sort_keys=True),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
