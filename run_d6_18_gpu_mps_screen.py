#!/usr/bin/env python3
"""Batched Apple-MPS realization triage for all 12,712 d=6 deletions.

This is a HEURISTIC-ONLY discovery screen.  Each graph is gauged by fixing a
deterministically chosen unit K6 to a centered regular simplex.  The other 12
vertices are optimized in simultaneous random starts.  Required edges are
fitted to squared distance one; candidate nonedges remain unconstrained.  A
soft short-distance repulsion and separate distinct/collapsed scores prevent
coincident numerical roots from being mistaken for realizations.

Large residual, optimizer failure, or exhaustion of every start proves
nothing.  A small residual is only a candidate for exact/interval follow-up.
Chunks are atomic and resumable, and every random stream is a deterministic
function of the corpus index, seed base, and restart wave.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import gzip
import hashlib
import json
import math
import multiprocessing
import os
import platform
import shlex
import sys
import tempfile
import time
from pathlib import Path
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
CORPUS_NAME = "d6_residue_18_deletions.json"
CORPUS_VERIFICATION_NAME = "d6_residue_18_deletions_verification.json"
EXPECTED_CORPUS_SHA256 = (
    "9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732"
)
EXPECTED_CORPUS_VERIFICATION_SHA256 = (
    "50ab5d88441b2869e60a05f3e1485db0b25008c3a816960a20a3bb21c1c4defd"
)
GRAPH_COUNT = 12_712
STANDARD_COMPATIBLE_COUNT = 14
NONSTANDARD_COUNT = 12_698
N = 18
D = 6
SEED_SIZE = 6
OUTSIDE = 12
PAIR_COUNT = N * (N - 1) // 2

DEFAULT_CHECKPOINT_ROOT = ROOT / "d6_18_gpu_mps_checkpoints"
DEFAULT_REPORT = ROOT / "d6_18_gpu_mps_screen_report.json"
DEFAULT_RESULTS = ROOT / "d6_18_gpu_mps_screen_results.tsv.gz"
DEFAULT_BENCHMARK = ROOT / "d6_18_gpu_mps_benchmark.json"

RESULT_COLUMNS = (
    "index",
    "edges",
    "standard_compatible",
    "gauge_size",
    "optimizer_starts",
    "best_edge_rms_including_collisions",
    "best_max_edge_error_including_collisions",
    "best_minimum_distance_including_collisions",
    "best_distinct_edge_rms",
    "best_distinct_max_edge_error",
    "best_distinct_minimum_distance",
    "initial_standard_control_edge_rms",
    "initial_standard_control_minimum_distance",
    "gpu_candidate",
    "cpu_lm_attempts",
    "cpu_lm_successful_terminations",
    "best_lm_edge_rms_including_collisions",
    "best_lm_max_edge_error_including_collisions",
    "best_lm_minimum_distance_including_collisions",
    "best_distinct_lm_edge_rms",
    "best_distinct_lm_max_edge_error",
    "best_distinct_lm_minimum_distance",
    "best_distinct_lm_input_rank",
    "first_lm_candidate_input_rank",
    "lm_candidate",
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def file_sha256(path: Path) -> str:
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
    descriptor, name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    except BaseException:
        try:
            os.unlink(name)
        except FileNotFoundError:
            pass
        raise


def atomic_gzip_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
                stream.write(
                    (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode(
                        "ascii"
                    )
                )
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(name, path)
    except BaseException:
        try:
            os.unlink(name)
        except FileNotFoundError:
            pass
        raise


def atomic_gzip_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
                stream.write(value.encode("ascii"))
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(name, path)
    except BaseException:
        try:
            os.unlink(name)
        except FileNotFoundError:
            pass
        raise


def validate_rows(value: object, label: str) -> tuple[int, ...]:
    if (
        not isinstance(value, list)
        or len(value) != N
        or any(type(row) is not int for row in value)
    ):
        raise ValueError(f"{label}: malformed adjacency")
    rows = tuple(value)
    for vertex, row in enumerate(rows):
        if row < 0 or row >= 1 << N or row & (1 << vertex):
            raise ValueError(f"{label}: invalid row {vertex}")
        for other in range(N):
            if bool(row & (1 << other)) != bool(rows[other] & (1 << vertex)):
                raise ValueError(f"{label}: asymmetric adjacency")
    return rows


def find_first_clique(rows: Sequence[int], target: int = SEED_SIZE) -> tuple[int, ...]:
    def search(chosen: tuple[int, ...], candidates: int) -> tuple[int, ...] | None:
        if len(chosen) == target:
            return chosen
        if len(chosen) + candidates.bit_count() < target:
            return None
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            result = search(chosen + (vertex,), candidates & rows[vertex])
            if result is not None:
                return result
        return None

    answer = search((), (1 << N) - 1)
    if answer is None:
        raise ValueError("18-deletion graph has no K6 gauge seed")
    return answer


def simplex_coordinates() -> np.ndarray:
    # First six vertices of one fixed regular K7 simplex in R^6.  Using the
    # same frame permits K6 and K7 gauges in a single padded MPS batch.
    return full_simplex_coordinates()[:SEED_SIZE]


def full_simplex_coordinates() -> np.ndarray:
    # Helmert basis for 1^perp in R^7, divided by sqrt(2), gives seven rows
    # whose pairwise distances are one.
    basis = np.zeros((7, 6), dtype=np.float64)
    for column in range(6):
        denominator = math.sqrt((column + 1) * (column + 2))
        basis[: column + 1, column] = 1.0 / denominator
        basis[column + 1, column] = -(column + 1) / denominator
    return basis / math.sqrt(2.0)


def standard_coordinates() -> np.ndarray:
    base = []
    for word in range(32):
        if word.bit_count() & 1:
            base.append(
                [
                    (1.0 if word & (1 << axis) else -1.0) / math.sqrt(8.0)
                    for axis in range(5)
                ]
                + [0.0]
            )
    height = math.sqrt(3.0 / 8.0)
    return np.asarray(base + [[0.0] * 5 + [height], [0.0] * 5 + [-height]])


def common_sphere_parameters(adjacency: np.ndarray, gauge_size: int) -> list[dict]:
    seed_points = full_simplex_coordinates()[:gauge_size]
    parameters = []
    for slot in range(OUTSIDE):
        if gauge_size == 7 and slot == 0:
            parameters.append({"fixed": seed_points[6]})
            continue
        vertex = SEED_SIZE + slot
        neighbours = np.flatnonzero(adjacency[vertex, :gauge_size] > 0.5)
        k = len(neighbours)
        if k == 0:
            parameters.append({"unconstrained": True})
            continue
        if k == 7:
            raise ValueError("K8 cannot be initialized in R^6")
        centre = seed_points[neighbours].mean(axis=0)
        if k == 1:
            complement = np.eye(D)
        else:
            differences = (
                seed_points[neighbours[1:]] - seed_points[neighbours[0]]
            )
            _, _, right_t = np.linalg.svd(differences, full_matrices=True)
            complement = right_t[k - 1 :].T
        if complement.shape[1] != 7 - k:
            raise AssertionError("regular-simplex common-sphere dimension drift")
        parameters.append(
            {
                "centre": centre,
                "complement": complement,
                "radius": math.sqrt((k + 1.0) / (2.0 * k)),
            }
        )
    return parameters


def align_standard_seed(
    record: dict, new_to_old: Sequence[int], simplex: np.ndarray
) -> np.ndarray | None:
    if not record["standard18_compatible"]:
        return None
    permutation = record["standard18_pole_pair_witnesses"][0][
        "embedding_permutation"
    ]
    source = standard_coordinates()[permutation]
    ordered = source[list(new_to_old)]
    seed = ordered[:SEED_SIZE]
    centre = seed.mean(axis=0)
    target_centre = simplex.mean(axis=0)
    covariance = (seed - centre).T @ (simplex - target_centre)
    left, _, right_t = np.linalg.svd(covariance, full_matrices=True)
    rotation = left @ right_t
    aligned = (ordered - centre) @ rotation + target_centre
    if np.max(np.abs(aligned[:SEED_SIZE] - simplex)) > 1e-10:
        raise AssertionError("standard control K6 alignment failed")
    return aligned


def enumerate_cliques(
    rows: Sequence[int], target: int, limit: int | None = None
) -> list[tuple[int, ...]]:
    answer: list[tuple[int, ...]] = []

    def search(chosen: tuple[int, ...], candidates: int) -> None:
        if limit is not None and len(answer) >= limit:
            return
        if len(chosen) == target:
            answer.append(chosen)
            return
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            remaining = candidates & rows[vertex]
            if len(chosen) + 1 + remaining.bit_count() >= target:
                search(chosen + (vertex,), remaining)

    search((), (1 << N) - 1)
    return answer


def reorder_record_with_seed(
    index: int, record: dict, seed: Sequence[int]
) -> dict:
    rows = validate_rows(record["adjacency"], record["class_id"])
    seed = tuple(seed)
    if len(seed) not in (6, 7) or tuple(sorted(seed)) != seed:
        raise ValueError("gauge seed must be a sorted K6 or K7")
    for first in seed:
        for second in seed:
            if first > second and not rows[first] & (1 << second):
                raise ValueError("gauge seed is not a clique")
    gauge_size = len(seed)
    remaining = tuple(vertex for vertex in range(N) if vertex not in seed)
    new_to_old = seed + remaining
    old_to_new = {old: new for new, old in enumerate(new_to_old)}
    adjacency = np.zeros((N, N), dtype=np.float32)
    edges = 0
    for old_first in range(N):
        for old_second in range(old_first):
            if rows[old_first] & (1 << old_second):
                first = old_to_new[old_first]
                second = old_to_new[old_second]
                adjacency[first, second] = adjacency[second, first] = 1.0
                edges += 1
    if edges != record["edges"]:
        raise AssertionError("edge count drift")
    simplex = simplex_coordinates()
    aligned = align_standard_seed(record, new_to_old, simplex)
    return {
        "index": index,
        "class_id": record["class_id"],
        "edges": edges,
        "standard_compatible": bool(record["standard18_compatible"]),
        "gauge_size": gauge_size,
        "seed_old_vertices": list(seed),
        "new_to_old": list(new_to_old),
        "adjacency": adjacency,
        "standard_aligned": aligned,
        "sphere_parameters": common_sphere_parameters(adjacency, gauge_size),
    }


def reorder_record(index: int, record: dict) -> dict:
    rows = validate_rows(record["adjacency"], record["class_id"])
    try:
        seed = find_first_clique(rows, 7)
    except ValueError:
        seed = find_first_clique(rows, 6)
    return reorder_record_with_seed(index, record, seed)


def load_corpus(path: Path) -> tuple[dict, list[dict]]:
    if file_sha256(path) != EXPECTED_CORPUS_SHA256:
        raise ValueError("pinned deletion corpus hash mismatch")
    verification_path = path.with_name(CORPUS_VERIFICATION_NAME)
    if file_sha256(verification_path) != EXPECTED_CORPUS_VERIFICATION_SHA256:
        raise ValueError("pinned deletion verification hash mismatch")
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    if verification.get("status") != "PASS":
        raise ValueError("deletion verification is not PASS")
    corpus = json.loads(path.read_text(encoding="utf-8"))
    if corpus.get("schema") != "d6-residue-18-deletion-manifest-v1":
        raise ValueError("unexpected deletion corpus schema")
    if len(corpus.get("unique_deletions", [])) != GRAPH_COUNT:
        raise ValueError("deletion corpus count drift")
    records = [
        reorder_record(index, record)
        for index, record in enumerate(corpus["unique_deletions"])
    ]
    if sum(record["standard_compatible"] for record in records) != STANDARD_COMPATIBLE_COUNT:
        raise ValueError("standard-compatible count drift")
    return corpus, records


def graph_seed(seedbase: int, corpus_index: int, wave: int) -> int:
    # PCG64 accepts arbitrary nonnegative Python integers.  Keeping the formula
    # in the report makes every graph/wave stream independent of batching.
    return (seedbase * 1_000_003 + corpus_index * 7_919 + wave * 104_729) % (1 << 64)


def seed_sphere_initial(record: dict, generator: np.random.Generator, scale: float) -> np.ndarray:
    """Place every outside vertex exactly on its seed-neighbour sphere."""
    values = np.empty((OUTSIDE, D), dtype=np.float64)
    for slot, parameters in enumerate(record["sphere_parameters"]):
        if "fixed" in parameters:
            values[slot] = parameters["fixed"]
            continue
        if parameters.get("unconstrained"):
            values[slot] = generator.normal(0.0, scale, size=D)
            continue
        complement = parameters["complement"]
        coefficients = generator.normal(size=complement.shape[1])
        norm = np.linalg.norm(coefficients)
        if norm == 0.0:
            coefficients[0] = 1.0
            norm = 1.0
        direction = complement @ (coefficients / norm)
        values[slot] = parameters["centre"] + parameters["radius"] * direction
    return values.astype(np.float32)


def initial_outside(
    records: Sequence[dict],
    restarts: int,
    wave: int,
    seedbase: int,
    scale: float,
    seed_standard_controls: bool,
    initialization: str,
) -> np.ndarray:
    answer = np.empty((len(records), restarts, OUTSIDE, D), dtype=np.float32)
    for batch_index, record in enumerate(records):
        generator = np.random.Generator(
            np.random.PCG64(graph_seed(seedbase, record["index"], wave))
        )
        if initialization == "seed_sphere":
            values = np.stack(
                [seed_sphere_initial(record, generator, scale) for _ in range(restarts)]
            )
        elif initialization == "gaussian":
            values = generator.normal(
                0.0, scale, size=(restarts, OUTSIDE, D)
            ).astype(np.float32)
        else:
            raise ValueError(f"unknown initialization {initialization}")
        if (
            seed_standard_controls
            and wave == 0
            and record["standard_aligned"] is not None
        ):
            values[0] = record["standard_aligned"][SEED_SIZE:].astype(np.float32)
        answer[batch_index] = values
    return answer


def numpy_metrics(points: np.ndarray, adjacency: np.ndarray) -> tuple[float, float, float]:
    differences = points[:, None, :] - points[None, :, :]
    distance_squared = np.sum(differences * differences, axis=2)
    edge_values = distance_squared[np.triu(adjacency.astype(bool), 1)] - 1.0
    pair_values = distance_squared[np.triu(np.ones((N, N), dtype=bool), 1)]
    return (
        float(np.sqrt(np.mean(edge_values * edge_values))),
        float(np.max(np.abs(edge_values))),
        float(np.sqrt(np.min(np.maximum(pair_values, 0.0)))),
    )


def optimize_wave(
    records: Sequence[dict],
    *,
    restarts: int,
    wave: int,
    steps: int,
    learning_rate: float,
    seedbase: int,
    initialization_scale: float,
    collision_distance: float,
    collision_weight: float,
    device_name: str,
    seed_standard_controls: bool,
    initialization: str,
) -> dict:
    import torch

    if device_name == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("Apple MPS unavailable; run outside the sandbox")
    torch.set_num_threads(1)
    device = torch.device(device_name)
    dtype = torch.float32
    count = len(records)
    adjacency_cpu = np.stack([record["adjacency"] for record in records])
    edge_counts_cpu = adjacency_cpu.sum(axis=(1, 2)) / 2.0
    adjacency = torch.as_tensor(adjacency_cpu, dtype=dtype, device=device)
    edge_counts = torch.as_tensor(edge_counts_cpu, dtype=dtype, device=device)
    upper = torch.triu(torch.ones((N, N), dtype=dtype, device=device), diagonal=1)
    simplex = torch.as_tensor(simplex_coordinates(), dtype=dtype, device=device)
    fixed = simplex[None, None, :, :].expand(count, restarts, -1, -1)
    initial = initial_outside(
        records, restarts, wave, seedbase, initialization_scale,
        seed_standard_controls, initialization,
    )
    outside = torch.nn.Parameter(torch.as_tensor(initial, dtype=dtype, device=device))
    free_mask_cpu = np.ones((count, 1, OUTSIDE, 1), dtype=np.float32)
    anchors_cpu = np.zeros((count, 1, OUTSIDE, D), dtype=np.float32)
    seventh = full_simplex_coordinates()[6].astype(np.float32)
    for local, record in enumerate(records):
        if record["gauge_size"] == 7:
            free_mask_cpu[local, 0, 0, 0] = 0.0
            anchors_cpu[local, 0, 0] = seventh
    free_mask = torch.as_tensor(free_mask_cpu, dtype=dtype, device=device)
    anchors = torch.as_tensor(anchors_cpu, dtype=dtype, device=device)
    optimizer = torch.optim.Adam([outside], lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=max(steps, 1), eta_min=learning_rate * 0.02
    )
    collision_squared = collision_distance * collision_distance

    def losses() -> tuple[object, object, object, object]:
        effective_outside = outside * free_mask + anchors * (1.0 - free_mask)
        points = torch.cat((fixed, effective_outside), dim=2)
        gram = points @ points.transpose(2, 3)
        norm = torch.diagonal(gram, dim1=2, dim2=3)
        distances = (norm[:, :, :, None] + norm[:, :, None, :] - 2.0 * gram).clamp_min(0.0)
        error = distances - 1.0
        edge_loss = (
            error.square() * adjacency[:, None, :, :]
        ).sum(dim=(2, 3)) / (2.0 * edge_counts[:, None])
        # This is a discovery bias, not a valid geometric constraint.  It is
        # smooth away from exact coincidence and is reported explicitly.
        repulsion = torch.exp(-distances / collision_squared) * upper
        collision_loss = repulsion.sum(dim=(2, 3)) / PAIR_COUNT
        objective = edge_loss + collision_weight * collision_loss
        return points, distances, edge_loss, objective

    initial_controls = {}
    with torch.no_grad():
        points0, distances0, edge0, _ = losses()
        for local, record in enumerate(records):
            if (
                seed_standard_controls
                and record["standard_compatible"]
                and wave == 0
            ):
                pair = distances0[local, 0][upper.bool()]
                initial_controls[record["index"]] = {
                    "edge_rms": float(torch.sqrt(edge0[local, 0]).cpu()),
                    "minimum_distance": float(torch.sqrt(pair.min()).cpu()),
                }

    started = time.perf_counter()
    for _ in range(steps):
        optimizer.zero_grad(set_to_none=True)
        _, _, _, objective = losses()
        objective.mean().backward()
        torch.nn.utils.clip_grad_norm_([outside], max_norm=10.0)
        optimizer.step()
        scheduler.step()
    if device_name == "mps":
        torch.mps.synchronize()
    elapsed = time.perf_counter() - started

    with torch.no_grad():
        points, distances, edge_loss, _ = losses()
        error = (distances - 1.0).abs() * adjacency[:, None, :, :]
        max_edge = error.amax(dim=(2, 3))
        pair_distances = torch.sqrt(
            torch.where(
                upper.bool(),
                distances,
                torch.full_like(distances, float("inf")),
            ).amin(dim=(2, 3))
        )
        edge_rms = torch.sqrt(edge_loss)
        points_cpu = points.cpu().numpy()
        edge_rms_cpu = edge_rms.cpu().numpy()
        max_edge_cpu = max_edge.cpu().numpy()
        minimum_cpu = pair_distances.cpu().numpy()

    return {
        "wave": wave,
        "restarts": restarts,
        "steps": steps,
        "wall_seconds": elapsed,
        "graph_start_steps": count * restarts * steps,
        "graph_start_steps_per_second": count * restarts * steps / max(elapsed, 1e-12),
        "edge_rms": edge_rms_cpu,
        "max_edge_error": max_edge_cpu,
        "minimum_distance": minimum_cpu,
        "coordinates": points_cpu,
        "initial_controls": initial_controls,
        "device": str(device),
        "dtype": str(dtype),
    }


def cpu_lm_refine(
    points: np.ndarray,
    adjacency: np.ndarray,
    gauge_size: int,
    max_nfev: int,
    collision_distance: float,
    collision_weight: float,
) -> dict:
    from scipy.optimize import least_squares

    edges = np.asarray([
        (first, second)
        for first in range(N)
        for second in range(first)
        if adjacency[first, second]
    ], dtype=np.int64)
    pairs = np.asarray(
        [(first, second) for first in range(N) for second in range(first)],
        dtype=np.int64,
    )
    fixed = full_simplex_coordinates()[:gauge_size]
    free_vertices = N - gauge_size
    collision_scale = math.sqrt(collision_weight)

    def residual(vector: np.ndarray) -> np.ndarray:
        current = np.concatenate((fixed, vector.reshape(free_vertices, D)), axis=0)
        differences = current[edges[:, 0]] - current[edges[:, 1]]
        edge_residual = np.einsum("ij,ij->i", differences, differences) - 1.0
        pair_differences = current[pairs[:, 0]] - current[pairs[:, 1]]
        pair_squared = np.einsum("ij,ij->i", pair_differences, pair_differences)
        collision_residual = collision_scale * np.maximum(
            collision_distance - np.sqrt(np.maximum(pair_squared, 0.0)), 0.0
        )
        return np.concatenate((edge_residual, collision_residual))

    def jacobian(vector: np.ndarray) -> np.ndarray:
        current = np.concatenate((fixed, vector.reshape(free_vertices, D)), axis=0)
        differences = current[edges[:, 0]] - current[edges[:, 1]]
        answer = np.zeros(
            (len(edges) + len(pairs), free_vertices * D), dtype=np.float64
        )
        for row, (first, second) in enumerate(edges):
            if first >= gauge_size:
                start = (first - gauge_size) * D
                answer[row, start : start + D] = 2.0 * differences[row]
            if second >= gauge_size:
                start = (second - gauge_size) * D
                answer[row, start : start + D] = -2.0 * differences[row]
        pair_differences = current[pairs[:, 0]] - current[pairs[:, 1]]
        pair_squared = np.einsum("ij,ij->i", pair_differences, pair_differences)
        for offset, (first, second) in enumerate(pairs):
            distance = math.sqrt(max(pair_squared[offset], 0.0))
            if distance >= collision_distance:
                continue
            row = len(edges) + offset
            direction = pair_differences[offset] / max(distance, 1e-12)
            if first >= gauge_size:
                start = (first - gauge_size) * D
                answer[row, start : start + D] = (
                    -collision_scale * direction
                )
            if second >= gauge_size:
                start = (second - gauge_size) * D
                answer[row, start : start + D] = (
                    collision_scale * direction
                )
        return answer

    solution = least_squares(
        residual,
        points[gauge_size:].astype(np.float64).ravel(),
        jac=jacobian,
        method="lm",
        x_scale="jac",
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
        max_nfev=max_nfev,
    )
    refined = np.concatenate((fixed, solution.x.reshape(free_vertices, D)), axis=0)
    metrics = numpy_metrics(refined, adjacency)
    return {
        "coordinates": refined,
        "edge_rms": metrics[0],
        "max_edge_error": metrics[1],
        "minimum_distance": metrics[2],
        "successful_termination": bool(solution.success),
        "status": int(solution.status),
        "nfev": int(solution.nfev),
        "njev": int(solution.njev) if solution.njev is not None else -1,
        "optimality": float(solution.optimality),
    }


def cpu_lm_refine_graph(task: tuple) -> tuple[int, list[dict]]:
    """Process-pool entry: refine one graph's ranked MPS endpoints."""
    (
        local, adjacency, gauge_size, candidates, max_nfev,
        collision_distance, collision_weight,
    ) = task
    results = []
    for input_rank, candidate in enumerate(candidates, start=1):
        refined = cpu_lm_refine(
            candidate["coordinates"], adjacency, gauge_size, max_nfev,
            collision_distance, collision_weight,
        )
        refined["input_rank"] = input_rank
        refined["input_edge_rms"] = candidate["edge_rms"]
        refined["input_minimum_distance"] = candidate["minimum_distance"]
        results.append(refined)
    return local, results


def restore_original_order(points: np.ndarray, new_to_old: Sequence[int]) -> list[list[float]]:
    answer = np.empty_like(points)
    for new, old in enumerate(new_to_old):
        answer[old] = points[new]
    return [[float(value) for value in row] for row in answer]


def process_batch(
    records: Sequence[dict],
    config: dict,
    executor: concurrent.futures.Executor | None = None,
) -> dict:
    best_all: list[dict | None] = [None] * len(records)
    top_distinct: list[list[dict]] = [[] for _ in records]
    initial_controls: dict[int, dict] = {}
    wave_profiles = []
    total_restarts = config["total_restarts"]
    per_wave = config["restarts_per_wave"]
    waves = math.ceil(total_restarts / per_wave)
    retain_count = max(1, config["lm_top_k"])

    for wave in range(waves):
        restarts = min(per_wave, total_restarts - wave * per_wave)
        result = optimize_wave(
            records,
            restarts=restarts,
            wave=wave,
            steps=config["optimizer_steps"],
            learning_rate=config["learning_rate"],
            seedbase=config["seedbase"],
            initialization_scale=config["initialization_scale"],
            collision_distance=config["collision_distance"],
            collision_weight=config["collision_weight"],
            device_name=config["device"],
            seed_standard_controls=config["seed_standard_controls"],
            initialization=config["initialization_method"],
        )
        initial_controls.update(result["initial_controls"])
        wave_profiles.append(
            {
                key: result[key]
                for key in (
                    "wave",
                    "restarts",
                    "steps",
                    "wall_seconds",
                    "graph_start_steps",
                    "graph_start_steps_per_second",
                    "device",
                    "dtype",
                )
            }
        )
        for local in range(len(records)):
            for restart in range(restarts):
                candidate = {
                    "edge_rms": float(result["edge_rms"][local, restart]),
                    "max_edge_error": float(result["max_edge_error"][local, restart]),
                    "minimum_distance": float(result["minimum_distance"][local, restart]),
                    "wave": wave,
                    "restart_in_wave": restart,
                    "coordinates": result["coordinates"][local, restart].copy(),
                }
                if best_all[local] is None or candidate["edge_rms"] < best_all[local]["edge_rms"]:
                    best_all[local] = candidate
                if candidate["minimum_distance"] >= config["distinctness_threshold"]:
                    top_distinct[local].append(candidate)
                    top_distinct[local].sort(key=lambda item: item["edge_rms"])
                    del top_distinct[local][retain_count:]

    lm_started = time.perf_counter()
    lm_by_local: list[list[dict]] = [[] for _ in records]
    if config["lm_top_k"] > 0 and config["cpu_refine_max_nfev"] > 0:
        tasks = [
            (
                local,
                record["adjacency"],
                record["gauge_size"],
                top_distinct[local][: config["lm_top_k"]],
                config["cpu_refine_max_nfev"],
                config["lm_collision_distance"],
                config["lm_collision_weight"],
            )
            for local, record in enumerate(records)
        ]
        mapped = executor.map(cpu_lm_refine_graph, tasks) if executor else map(
            cpu_lm_refine_graph, tasks
        )
        for local, refinements in mapped:
            lm_by_local[local] = refinements
    lm_profile = {
        "wall_seconds": time.perf_counter() - lm_started,
        "workers": config["cpu_workers"],
        "top_k": config["lm_top_k"],
        "attempts": sum(map(len, lm_by_local)),
        "method": "SciPy float64 least_squares(method=lm), analytic Jacobian, and short-distance collision hinge",
    }

    rows = []
    witnesses = []
    retained_endpoints = []
    lm_attempts = []
    for local, record in enumerate(records):
        all_candidate = best_all[local]
        distinct = top_distinct[local][0] if top_distinct[local] else None
        if all_candidate is None:
            raise AssertionError("empty GPU restart collection")
        control = initial_controls.get(record["index"])
        gpu_candidate = bool(
            distinct is not None
            and distinct["edge_rms"] <= config["gpu_candidate_rms"]
        )
        refinements = lm_by_local[local]
        best_lm = min(refinements, key=lambda item: item["edge_rms"], default=None)
        distinct_lm = [
            item for item in refinements
            if item["minimum_distance"] >= config["distinctness_threshold"]
        ]
        best_distinct_lm = min(
            distinct_lm, key=lambda item: item["edge_rms"], default=None
        )
        lm_candidate = bool(
            best_distinct_lm is not None
            and best_distinct_lm["edge_rms"] <= config["lm_candidate_rms"]
        )
        candidate_input_ranks = [
            item["input_rank"] for item in distinct_lm
            if item["edge_rms"] <= config["lm_candidate_rms"]
        ]
        row = {
            "index": record["index"],
            "class_id": record["class_id"],
            "edges": record["edges"],
            "standard_compatible": int(record["standard_compatible"]),
            "gauge_size": record["gauge_size"],
            "gauge_seed_old_vertices": record["seed_old_vertices"],
            "optimizer_starts": config["total_restarts"],
            "best_edge_rms_including_collisions": all_candidate["edge_rms"],
            "best_max_edge_error_including_collisions": all_candidate["max_edge_error"],
            "best_minimum_distance_including_collisions": all_candidate["minimum_distance"],
            "best_distinct_edge_rms": distinct["edge_rms"] if distinct else math.inf,
            "best_distinct_max_edge_error": distinct["max_edge_error"] if distinct else math.inf,
            "best_distinct_minimum_distance": distinct["minimum_distance"] if distinct else 0.0,
            "initial_standard_control_edge_rms": control["edge_rms"] if control else math.nan,
            "initial_standard_control_minimum_distance": control["minimum_distance"] if control else math.nan,
            "gpu_candidate": int(gpu_candidate),
            "cpu_lm_attempts": len(refinements),
            "cpu_lm_successful_terminations": sum(
                item["successful_termination"] for item in refinements
            ),
            "best_lm_edge_rms_including_collisions": best_lm["edge_rms"] if best_lm else math.inf,
            "best_lm_max_edge_error_including_collisions": best_lm["max_edge_error"] if best_lm else math.inf,
            "best_lm_minimum_distance_including_collisions": best_lm["minimum_distance"] if best_lm else 0.0,
            "best_distinct_lm_edge_rms": best_distinct_lm["edge_rms"] if best_distinct_lm else math.inf,
            "best_distinct_lm_max_edge_error": best_distinct_lm["max_edge_error"] if best_distinct_lm else math.inf,
            "best_distinct_lm_minimum_distance": best_distinct_lm["minimum_distance"] if best_distinct_lm else 0.0,
            "best_distinct_lm_input_rank": best_distinct_lm["input_rank"] if best_distinct_lm else -1,
            "first_lm_candidate_input_rank": min(candidate_input_ranks, default=-1),
            "lm_candidate": int(lm_candidate),
        }
        rows.append(row)
        if distinct is not None:
            retained_endpoints.append(
                {
                    "index": record["index"],
                    "edge_rms": distinct["edge_rms"],
                    "max_edge_error": distinct["max_edge_error"],
                    "minimum_distance": distinct["minimum_distance"],
                    "wave": distinct["wave"],
                    "restart_in_wave": distinct["restart_in_wave"],
                    "coordinates_original_vertex_order": restore_original_order(
                        distinct["coordinates"], record["new_to_old"]
                    ),
                }
            )
        lm_attempts.append(
            {
                "index": record["index"],
                "attempts": [
                    {
                        key: item[key]
                        for key in (
                            "input_rank", "input_edge_rms", "input_minimum_distance",
                            "edge_rms", "max_edge_error", "minimum_distance",
                            "successful_termination", "status", "nfev", "njev",
                            "optimality",
                        )
                    }
                    for item in refinements
                ],
            }
        )
        if lm_candidate or gpu_candidate:
            selected = best_distinct_lm if lm_candidate else distinct
            selected_points = selected["coordinates"]
            selected_metrics = (
                selected["edge_rms"], selected["max_edge_error"],
                selected["minimum_distance"],
            ) if lm_candidate else (
                distinct["edge_rms"],
                distinct["max_edge_error"],
                distinct["minimum_distance"],
            )
            witnesses.append(
                {
                    "index": record["index"],
                    "class_id": record["class_id"],
                    "source": "CPU_float64_analytic_LM" if lm_candidate else "MPS_float32",
                    "edge_rms": selected_metrics[0],
                    "max_edge_error": selected_metrics[1],
                    "minimum_distance": selected_metrics[2],
                    "lm_input_rank": selected["input_rank"] if lm_candidate else None,
                    "coordinates_original_vertex_order": restore_original_order(
                        selected_points, record["new_to_old"]
                    ),
                }
            )
    return {
        "rows": rows,
        "witnesses": witnesses,
        "retained_endpoints": retained_endpoints,
        "lm_attempts": lm_attempts,
        "wave_profiles": wave_profiles,
        "lm_profile": lm_profile,
    }


def configuration(args: argparse.Namespace, source_hash: str) -> dict:
    return {
        "schema": "d6-18-gpu-mps-screen-config-v1",
        "corpus_sha256": EXPECTED_CORPUS_SHA256,
        "corpus_verification_sha256": EXPECTED_CORPUS_VERIFICATION_SHA256,
        "orchestrator_source_sha256": source_hash,
        "graph_count": GRAPH_COUNT,
        "dimension": D,
        "gauge": "lexicographically first K7 when present, otherwise first K6; K6 frame is the first six vertices of the fixed regular K7 simplex",
        "total_restarts": args.total_restarts,
        "restarts_per_wave": args.restarts_per_wave,
        "optimizer": "PyTorch Adam with cosine learning-rate decay",
        "optimizer_steps": args.steps,
        "learning_rate": args.learning_rate,
        "initialization_method": args.initialization,
        "initialization": (
            "each outside vertex sampled exactly on the common unit sphere of its fixed-simplex neighbours"
            if args.initialization == "seed_sphere" else
            "independent NumPy PCG64 Gaussian outside coordinates"
        ),
        "seed_standard_controls": args.seed_standard_controls,
        "control_seeding": (
            "for each of 14 standard-compatible graphs only, optimizer start 0 of wave 0 is replaced by an aligned exact standard-18 realization"
            if args.seed_standard_controls else
            "disabled for unseeded recovery calibration"
        ),
        "initialization_scale": args.initialization_scale,
        "seedbase": args.seedbase,
        "graph_seed_formula": "(seedbase*1000003 + corpus_index*7919 + wave*104729) mod 2^64",
        "collision_handling": "soft exp(-distance_squared/collision_distance^2) repulsion plus separate distinct score",
        "collision_distance": args.collision_distance,
        "collision_weight": args.collision_weight,
        "distinctness_threshold": args.distinctness_threshold,
        "gpu_candidate_rms": args.gpu_candidate_rms,
        "cpu_refine_max_nfev": args.cpu_refine_max_nfev,
        "cpu_refinement": "top-k distinct MPS endpoints per graph; SciPy float64 LM with exact analytic Jacobian and a zero-outside-radius distance collision hinge",
        "lm_top_k": args.lm_top_k,
        "lm_candidate_rms": args.lm_candidate_rms,
        "lm_collision_distance": args.lm_collision_distance,
        "lm_collision_weight": args.lm_collision_weight,
        "cpu_workers": args.cpu_workers,
        "chunk_size": args.chunk_size,
        "device": args.device,
        "dtype": "torch.float32 on accelerator; candidate refinement uses NumPy/SciPy float64",
        "semantics": "HEURISTIC_ONLY; no numerical rejection or realizability conclusion",
    }


def checkpoint_path(run_directory: Path, start: int, end: int) -> Path:
    return run_directory / f"chunk_{start:05d}_{end:05d}.json.gz"


def load_checkpoint(path: Path, config_hash: str, start: int, end: int) -> dict:
    with gzip.open(path, "rt", encoding="ascii") as stream:
        payload = json.load(stream)
    if payload.get("schema") != "d6-18-gpu-mps-checkpoint-v1":
        raise ValueError(f"{path}: checkpoint schema")
    if payload.get("config_sha256") != config_hash:
        raise ValueError(f"{path}: checkpoint configuration mismatch")
    if payload.get("range") != [start, end]:
        raise ValueError(f"{path}: checkpoint range mismatch")
    rows = payload.get("rows")
    if not isinstance(rows, list) or [row.get("index") for row in rows] != list(range(start, end)):
        raise ValueError(f"{path}: checkpoint row coverage")
    return payload


def benchmark_indices(records: Sequence[dict], count: int, seedbase: int) -> list[int]:
    controls = [record["index"] for record in records if record["standard_compatible"]]
    nonstandard = [record["index"] for record in records if not record["standard_compatible"]]
    wanted = max(0, count - len(controls))
    ranked = sorted(
        nonstandard,
        key=lambda index: hashlib.sha256(f"{seedbase}:{index}".encode("ascii")).digest(),
    )[:wanted]
    return sorted(controls + ranked)


def make_executor(config: dict) -> concurrent.futures.ProcessPoolExecutor | None:
    if config["cpu_workers"] <= 1 or config["lm_top_k"] <= 0:
        return None
    return concurrent.futures.ProcessPoolExecutor(
        max_workers=config["cpu_workers"],
        mp_context=multiprocessing.get_context("spawn"),
    )


def run_benchmark(
    records: Sequence[dict], config: dict, config_hash: str, args: argparse.Namespace
) -> dict:
    indices = benchmark_indices(records, args.benchmark_count, args.seedbase)
    selected = [records[index] for index in indices]
    started = time.perf_counter()
    executor = make_executor(config)
    try:
        result = process_batch(selected, config, executor)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    elapsed = time.perf_counter() - started
    profile_steps = sum(
        wave["graph_start_steps"] for wave in result["wave_profiles"]
    )
    gpu_candidates = [
        row["index"] for row in result["rows"] if row["gpu_candidate"]
    ]
    lm_candidates = [
        row["index"] for row in result["rows"] if row["lm_candidate"]
    ]
    control_recovery = [
        {
            key: row[key]
            for key in (
                "index", "best_distinct_edge_rms", "cpu_lm_attempts",
                "best_distinct_lm_edge_rms", "best_distinct_lm_minimum_distance",
                "best_distinct_lm_input_rank", "first_lm_candidate_input_rank",
                "lm_candidate",
            )
        }
        for row in result["rows"] if row["standard_compatible"]
    ]
    prefix_values = sorted(
        set(
            value for value in (1, 2, 4, 8, 16, 32, 64, 128, config["lm_top_k"])
            if 0 < value <= config["lm_top_k"]
        )
    )
    control_prefix_recovery = {
        str(prefix): sum(
            0 < row["first_lm_candidate_input_rank"] <= prefix
            for row in result["rows"] if row["standard_compatible"]
        )
        for prefix in prefix_values
    }
    report = {
        "schema": "d6-18-gpu-mps-benchmark-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
        "created_utc": utc_now(),
        "config": config,
        "config_sha256": config_hash,
        "indices": indices,
        "graphs": len(indices),
        "gauge_distribution": {
            "K6": sum(row["gauge_size"] == 6 for row in result["rows"]),
            "K7": sum(row["gauge_size"] == 7 for row in result["rows"]),
        },
        "wall_seconds": elapsed,
        "graph_starts": len(indices) * config["total_restarts"],
        "graph_starts_per_second": (
            len(indices) * config["total_restarts"] / max(elapsed, 1e-12)
        ),
        "graph_start_steps_per_second": profile_steps / max(elapsed, 1e-12),
        "projected_full_wall_seconds": elapsed * GRAPH_COUNT / max(len(indices), 1),
        "gpu_candidate_indices": gpu_candidates,
        "lm_candidate_indices": lm_candidates,
        "standard_control_initial_pass": (
            all(
                row["initial_standard_control_edge_rms"] < 2e-6
                and row["initial_standard_control_minimum_distance"]
                > config["distinctness_threshold"]
                for row in result["rows"] if row["standard_compatible"]
            ) if config["seed_standard_controls"] else None
        ),
        "standard_control_lm_recovery": control_recovery,
        "standard_control_lm_recovered": sum(
            row["lm_candidate"] for row in result["rows"]
            if row["standard_compatible"]
        ),
        "standard_control_lm_prefix_recovery": control_prefix_recovery,
        "lm_candidate_count_by_gauge": {
            "K6": sum(
                row["lm_candidate"] for row in result["rows"]
                if row["gauge_size"] == 6
            ),
            "K7": sum(
                row["lm_candidate"] for row in result["rows"]
                if row["gauge_size"] == 7
            ),
        },
        "wave_profiles": result["wave_profiles"],
        "lm_profile": result["lm_profile"],
        "machine": {
            "platform": platform.platform(),
            "python": sys.version,
        },
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }
    atomic_json(args.benchmark_output, report)
    return report


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return format(value, ".17g")


def assemble(
    checkpoints: Sequence[dict], config: dict, config_hash: str, args: argparse.Namespace
) -> dict:
    rows = [row for checkpoint in checkpoints for row in checkpoint["rows"]]
    witnesses = [
        witness for checkpoint in checkpoints for witness in checkpoint["witnesses"]
    ]
    if [row["index"] for row in rows] != list(range(GRAPH_COUNT)):
        raise ValueError("cannot assemble incomplete GPU sweep")
    header = "\t".join(RESULT_COLUMNS)
    lines = [header]
    integer_fields = {
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
    for row in rows:
        lines.append(
            "\t".join(
                str(row[column])
                if column in integer_fields
                else format_float(float(row[column]))
                for column in RESULT_COLUMNS
            )
        )
    atomic_gzip_text(args.results, "\n".join(lines) + "\n")
    threshold_counts = {
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
    lm_candidates = [row["index"] for row in rows if row["lm_candidate"]]
    report = {
        "schema": "d6-18-gpu-mps-screen-report-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
        "created_utc": utc_now(),
        "config": config,
        "config_sha256": config_hash,
        "corpus": {
            "path": CORPUS_NAME,
            "sha256": EXPECTED_CORPUS_SHA256,
            "verification_path": CORPUS_VERIFICATION_NAME,
            "verification_sha256": EXPECTED_CORPUS_VERIFICATION_SHA256,
        },
        "source": {
            "path": Path(__file__).name,
            "sha256": config["orchestrator_source_sha256"],
        },
        "results": {
            "path": args.results.name,
            "sha256": file_sha256(args.results),
            "columns": list(RESULT_COLUMNS),
            "rows": len(rows),
        },
        "checkpoints": {
            "directory": str(args.checkpoint_root / config_hash[:16]),
            "count": len(checkpoints),
            "ordered_sha256": stable_hash(
                [checkpoint["checkpoint_sha256"] for checkpoint in checkpoints]
            ),
        },
        "summary": {
            "graphs": len(rows),
            "gauge_distribution": {
                "K6": sum(row["gauge_size"] == 6 for row in rows),
                "K7": sum(row["gauge_size"] == 7 for row in rows),
            },
            "standard_compatible": sum(row["standard_compatible"] for row in rows),
            "nonstandard_compatible": NONSTANDARD_COUNT,
            "total_optimizer_starts": len(rows) * config["total_restarts"],
            "seeded_exact_control_starts": (
                STANDARD_COMPATIBLE_COUNT if config["seed_standard_controls"] else 0
            ),
            "random_starts": (
                len(rows) * config["total_restarts"]
                - (STANDARD_COMPATIBLE_COUNT if config["seed_standard_controls"] else 0)
            ),
            "gpu_distinct_rms_threshold_counts": threshold_counts,
            "lm_distinct_rms_threshold_counts": lm_threshold_counts,
            "gpu_candidate_indices": [row["index"] for row in rows if row["gpu_candidate"]],
            "lm_candidate_indices": lm_candidates,
            "retained_best_distinct_endpoint_count": sum(
                len(checkpoint["retained_endpoints"]) for checkpoint in checkpoints
            ),
            "cpu_lm_attempts": sum(row["cpu_lm_attempts"] for row in rows),
            "candidate_witnesses": witnesses,
            "mathematical_rejections": 0,
            "mathematical_realizability_conclusions": 0,
        },
        "trust": {
            "floating_point": "PyTorch/MPS float32 and optional SciPy float64; no rigorous rounding claim",
            "optimizer": "failure to find a solution has no mathematical meaning",
            "collision_bias": "soft repulsion may miss genuine near-collision realizations",
            "nonedges": "unconstrained and may also be unit",
        },
    }
    atomic_json(args.report, report)
    return report


def run_full(records: Sequence[dict], config: dict, config_hash: str, args: argparse.Namespace) -> dict | None:
    run_directory = args.checkpoint_root / config_hash[:16]
    run_directory.mkdir(parents=True, exist_ok=True)
    state = {
        "schema": "d6-18-gpu-mps-run-state-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
        "config": config,
        "config_sha256": config_hash,
        "pid": os.getpid(),
        "command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "source_sha256": config["orchestrator_source_sha256"],
        "started_utc": utc_now(),
        "checkpoint_directory": str(run_directory),
        "resume_command": shlex.join(
            [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]
        ),
    }
    atomic_json(run_directory / "run_state.json", state)
    checkpoints = []
    executor = make_executor(config)
    for start in range(0, GRAPH_COUNT, config["chunk_size"]):
        end = min(GRAPH_COUNT, start + config["chunk_size"])
        path = checkpoint_path(run_directory, start, end)
        if path.exists():
            payload = load_checkpoint(path, config_hash, start, end)
            checkpoints.append(payload)
            print(f"GPU resume [{start},{end})", flush=True)
            continue
        print(f"GPU start [{start},{end})", flush=True)
        started = time.perf_counter()
        try:
            result = process_batch(records[start:end], config, executor)
        except BaseException:
            if executor is not None:
                executor.shutdown(wait=False, cancel_futures=True)
            raise
        payload = {
            "schema": "d6-18-gpu-mps-checkpoint-v1",
            "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
            "config_sha256": config_hash,
            "range": [start, end],
            "created_utc": utc_now(),
            "wall_seconds": time.perf_counter() - started,
            "rows": result["rows"],
            "witnesses": result["witnesses"],
            "retained_endpoints": result["retained_endpoints"],
            "lm_attempts": result["lm_attempts"],
            "wave_profiles": result["wave_profiles"],
            "lm_profile": result["lm_profile"],
            "mathematical_rejections": 0,
        }
        # Hash the mathematical payload, excluding the self-referential field.
        payload["checkpoint_sha256"] = stable_hash(payload)
        atomic_gzip_json(path, payload)
        checkpoints.append(load_checkpoint(path, config_hash, start, end))
        print(
            f"GPU done [{start},{end}) {payload['wall_seconds']:.2f}s "
            f"gpu_candidates={sum(row['gpu_candidate'] for row in payload['rows'])} "
            f"lm_candidates={sum(row['lm_candidate'] for row in payload['rows'])}",
            flush=True,
        )
    if executor is not None:
        executor.shutdown(wait=True)
    if len(checkpoints) == math.ceil(GRAPH_COUNT / config["chunk_size"]):
        report = assemble(checkpoints, config, config_hash, args)
        print(json.dumps({"status": "COMPLETE", "summary": report["summary"]}, indent=2))
        return report
    return None


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    mode = result.add_mutually_exclusive_group(required=True)
    mode.add_argument("--benchmark", action="store_true")
    mode.add_argument("--full", action="store_true")
    result.add_argument("--corpus", type=Path, default=ROOT / CORPUS_NAME)
    result.add_argument("--checkpoint-root", type=Path, default=DEFAULT_CHECKPOINT_ROOT)
    result.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    result.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    result.add_argument("--benchmark-output", type=Path, default=DEFAULT_BENCHMARK)
    result.add_argument("--benchmark-count", type=int, default=64)
    result.add_argument("--device", choices=("mps", "cpu"), default="mps")
    result.add_argument("--total-restarts", type=int, default=128)
    result.add_argument("--restarts-per-wave", type=int, default=32)
    result.add_argument("--steps", type=int, default=1000)
    result.add_argument("--learning-rate", type=float, default=0.03)
    result.add_argument("--initialization-scale", type=float, default=0.7)
    result.add_argument(
        "--initialization", choices=("seed_sphere", "gaussian"),
        default="seed_sphere",
    )
    result.add_argument("--seedbase", type=int, default=600_180_777)
    result.add_argument(
        "--no-seed-standard-controls",
        action="store_false",
        dest="seed_standard_controls",
        help="disable the exact seeded start to measure unseeded recovery",
    )
    result.add_argument("--collision-distance", type=float, default=0.12)
    result.add_argument("--collision-weight", type=float, default=0.05)
    result.add_argument("--distinctness-threshold", type=float, default=0.01)
    result.add_argument("--gpu-candidate-rms", type=float, default=3e-4)
    result.add_argument("--cpu-refine-max-nfev", type=int, default=2000)
    result.add_argument("--lm-top-k", type=int, default=8)
    result.add_argument("--lm-candidate-rms", type=float, default=1e-8)
    result.add_argument("--lm-collision-distance", type=float, default=0.05)
    result.add_argument("--lm-collision-weight", type=float, default=100.0)
    result.add_argument("--cpu-workers", type=int, default=2)
    result.add_argument("--chunk-size", type=int, default=64)
    return result


def main() -> None:
    args = parser().parse_args()
    if args.total_restarts <= 0 or args.restarts_per_wave <= 0:
        raise ValueError("restart counts must be positive")
    if args.steps <= 0 or args.chunk_size <= 0:
        raise ValueError("steps and chunk size must be positive")
    if args.lm_top_k < 0 or args.cpu_workers <= 0:
        raise ValueError("LM top-k must be nonnegative and CPU workers positive")
    if args.lm_collision_distance <= 0 or args.lm_collision_weight < 0:
        raise ValueError("LM collision distance must be positive and weight nonnegative")
    source_hash = file_sha256(Path(__file__))
    config = configuration(args, source_hash)
    config_hash = stable_hash(config)
    _, records = load_corpus(args.corpus)
    if args.benchmark:
        report = run_benchmark(records, config, config_hash, args)
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        run_full(records, config, config_hash, args)


if __name__ == "__main__":
    main()
