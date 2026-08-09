#!/usr/bin/env python3
"""Restartable Apple-MPS triage of the published 260-case d=6 n=19 residue.

This is a numerical discovery campaign, not a proof campaign.  It fits only
the required unit edges.  Candidate nonedges are absent from the edge
residual and may also have distance one.  Failure to find a low-residual,
distinct endpoint has no mathematical meaning, while any endpoint found here
still needs exact or rigorous interval certification.

The frozen input physically retains the verified 261-case v8 residue.  The
verified exact rejection 3936435 is tagged as a negative calibration and is
never counted among the 260 active unresolved cases.  Every K7 seed is used
on the twelve physical K7 records; every K6 seed is used on the 249 K6-only
records.  The resulting 7,721 variants have chunk-independent random streams
and atomic checkpoints.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import math
import multiprocessing
import os
import platform
import shlex
import subprocess
import sys
import time
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np

import run_d6_18_gpu_mps_screen as kernel


ROOT = Path(__file__).resolve().parent
N = 19
D = 6
SEED_SIZE = 6
OUTSIDE = N - SEED_SIZE
PAIR_COUNT = N * (N - 1) // 2
VARIANT_STRIDE = 4_096
CONTROL_VARIANT_BASE = 1 << 62

V8_MANIFEST = ROOT / "d6_current_residue_manifest_v8.json"
V8_VERIFICATION = ROOT / "d6_current_residue_manifest_v8_verification.json"
K7_INCREMENT_REPORT = ROOT / "d6_k7_triple_petals_3936435_increment_report.json"
K7_INCREMENT_VERIFICATION = (
    ROOT / "d6_k7_triple_petals_3936435_increment_verification.json"
)

EXPECTED_KERNEL_SOURCE_SHA256 = (
    "cd57a6193fba1a9b05b88be566598413bd7ff79a13cdc7b5530df90a14293627"
)
EXPECTED_INPUT_SHA256 = {
    "v8_manifest": "9ea10a7794f033e66152c477a130f11c6c2862e87b4bdd705c14021521c18285",
    "v8_verification": "117ac6833a24eb69cec7a514445bb4c2913c3f3ac7c413aa81fad19fe305602e",
    "k7_increment_report": "b13240a303f3e304dd2d095cc6dd28f3caa1f6571e0ba210637cb7e0ab0f168a",
    "k7_increment_verification": "d07916b38ec2e640d06eaf44fe351e92d32b32f522455d5b5ff008d233ad7c0e",
}
EXPECTED_PHYSICAL_GRAPH_COUNT = 261
EXPECTED_ACTIVE_GRAPH_COUNT = 260
EXPECTED_PHYSICAL_CLASS_COUNTS = {"K7": 12, "K6_only": 249}
EXPECTED_ACTIVE_CLASS_COUNTS = {"K7": 11, "K6_only": 249}
EXPECTED_CAMPAIGN_INDICES_SHA256 = (
    "cdec394afee80a54e049e19fcadcfd0eb382b192872fae2d0e24ce42b1e8cec1"
)
EXPECTED_PHYSICAL_INDICES_SHA256 = (
    "2391a93a3629517363106603bdad00be9b6f960966d9089ef37d93be2988c213"
)
EXPECTED_PHYSICAL_TAGGED_GRAPHS_SHA256 = (
    "fbe20738cec808da2f2f9b4fe2f4f12399c109d85836d6bf6e539115e20700ca"
)
EXPECTED_GAUGE_COUNT = 7_721
EXPECTED_GAUGE_POPULATIONS = {
    "K7_active": 34,
    "K7_exact_negative_calibration": 3,
    "K6_only_active": 7_684,
}
EXPECTED_GAUGE_MANIFEST_SHA256 = (
    "a545686fc30b230473bcfcc091a826bd3d781fe0079163d4f1a8a31d04c10ca3"
)
EXPECTED_VARIANT_IDS_SHA256 = (
    "7aa068c6725716ed4dcd034b949609e4d2c6464c06a8ad9f34ca14668a34777f"
)
EXPECTED_GAUGE_CENSUS_SHA256 = (
    "c8af715a608e5b3b1e37515b90228ae8f41839a0c90a010d5d84c6fe1e8baf50"
)
EXPECTED_INTRINSIC_DOF_PROFILE_SHA256 = (
    "07be54ffed9f0c21f8cd326eb18043b10aef953c537670581b3a561ddbff13ca"
)
EXACTLY_REMOVED_INDEX = 3_936_435

DEFAULT_CHECKPOINT_ROOT = ROOT / ".runs/d6_19_gpu_mps_residue_v1"
DEFAULT_REPORT = ROOT / ".runs/d6_19_gpu_mps_residue_v1_report.json"
DEFAULT_DETAILS = ROOT / ".runs/d6_19_gpu_mps_residue_v1_details.json.gz"
LM_THREAD_ENVIRONMENT = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}

_BASE_CPU_LM_REFINE_GRAPH = kernel.cpu_lm_refine_graph
_BASE_OPTIMIZE_WAVE = kernel.optimize_wave


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


@contextmanager
def n19_kernel_context() -> Iterator[None]:
    """Temporarily specialize the pinned generic numerical kernel to n=19.

    The imported n=18 implementation computes all array sizes from module
    constants at call time.  Keeping the specialization local prevents this
    runner from contaminating accelerator-free tests of the n=18 campaign.
    """

    saved = {
        "N": kernel.N,
        "OUTSIDE": kernel.OUTSIDE,
        "PAIR_COUNT": kernel.PAIR_COUNT,
        "cpu_lm_refine_graph": kernel.cpu_lm_refine_graph,
        "optimize_wave": kernel.optimize_wave,
    }
    kernel.N = N
    kernel.OUTSIDE = OUTSIDE
    kernel.PAIR_COUNT = PAIR_COUNT
    kernel.cpu_lm_refine_graph = cpu_lm_refine_graph_n19
    kernel.optimize_wave = optimize_wave_manifold_n19
    try:
        yield
    finally:
        kernel.N = saved["N"]
        kernel.OUTSIDE = saved["OUTSIDE"]
        kernel.PAIR_COUNT = saved["PAIR_COUNT"]
        kernel.cpu_lm_refine_graph = saved["cpu_lm_refine_graph"]
        kernel.optimize_wave = saved["optimize_wave"]


def cpu_lm_refine_graph_n19(task: tuple) -> tuple[int, list[dict]]:
    """Spawn-safe entry that activates n=19 constants in every LM worker."""

    with n19_kernel_context():
        return _BASE_CPU_LM_REFINE_GRAPH(task)


def manifold_parameters(adjacency: np.ndarray, gauge_size: int) -> tuple[list[dict], dict]:
    """Build intrinsic parameters and resolve zero-dimensional K7 links.

    For k fixed-simplex neighbours, the unit-sphere intersection has an
    orthonormal coefficient space of dimension 7-k.  Its sphere therefore has
    6-k free dimensions.  In a K7 gauge, k=6 gives two points: the missing
    simplex vertex and its reflection through the K6 face.  Distinctness
    excludes the former, so the latter is fixed rather than randomized.
    """

    with n19_kernel_context():
        raw = kernel.common_sphere_parameters(adjacency, gauge_size)
    seed_points = kernel.full_simplex_coordinates()[:gauge_size]
    answer = []
    s0_reflections = 0
    for slot, original in enumerate(raw):
        if "fixed" in original:
            answer.append(
                {
                    "mode": "fixed",
                    "point": np.asarray(original["fixed"], dtype=np.float64),
                    "reason": "seventh_gauge_vertex",
                }
            )
            continue
        if original.get("unconstrained"):
            answer.append({"mode": "cartesian", "intrinsic_dimension": D})
            continue
        centre = np.asarray(original["centre"], dtype=np.float64)
        basis = np.asarray(original["complement"], dtype=np.float64)
        radius = float(original["radius"])
        coefficient_dimension = basis.shape[1]
        require(1 <= coefficient_dimension <= D, "coefficient-space dimension")
        if coefficient_dimension == 1:
            options = (centre + radius * basis[:, 0], centre - radius * basis[:, 0])
            if gauge_size == 7:
                vertex = SEED_SIZE + slot
                neighbours = set(np.flatnonzero(adjacency[vertex, :gauge_size] > 0.5))
                missing = sorted(set(range(gauge_size)) - neighbours)
                require(len(missing) == 1, "K7 S0 link does not miss exactly one seed")
                missing_point = seed_points[missing[0]]
                distances = [float(np.linalg.norm(option - missing_point)) for option in options]
                duplicate = int(np.argmin(distances))
                reflected = 1 - duplicate
                require(distances[duplicate] < 2e-10, "K7 S0 duplicate root drift")
                require(distances[reflected] > 0.5, "K7 S0 reflected root drift")
                answer.append(
                    {
                        "mode": "fixed",
                        "point": options[reflected],
                        "reason": "distinct_reflection_of_missing_K7_seed",
                        "missing_seed_slot": missing[0],
                    }
                )
                s0_reflections += 1
            else:
                # The two roots are congruent under the global reflection that
                # fixes a K6 gauge pointwise.  No production K6-only gauge in
                # the pinned census reaches this branch, but fixing one root
                # is complete modulo that residual gauge symmetry.
                answer.append(
                    {
                        "mode": "fixed",
                        "point": options[0],
                        "reason": "K6_gauge_reflection_symmetry_representative",
                    }
                )
            continue
        answer.append(
            {
                "mode": "sphere",
                "centre": centre,
                "basis": basis,
                "radius": radius,
                "coefficient_dimension": coefficient_dimension,
                "intrinsic_dimension": coefficient_dimension - 1,
            }
        )

    fixed_points = [
        (slot, parameters["point"])
        for slot, parameters in enumerate(answer)
        if parameters["mode"] == "fixed"
    ]
    for slot, point in fixed_points:
        for seed_slot, seed_point in enumerate(kernel.simplex_coordinates()):
            require(
                float(np.linalg.norm(point - seed_point)) > 1e-9,
                f"forced seed collision at outside slot {slot}, seed slot {seed_slot}",
            )
    for position, (first_slot, first) in enumerate(fixed_points):
        for second_slot, second in fixed_points[:position]:
            require(
                float(np.linalg.norm(first - second)) > 1e-9,
                f"forced-twin collision at outside slots {second_slot},{first_slot}",
            )
    intrinsic_free_dof = sum(
        parameters.get("intrinsic_dimension", 0) for parameters in answer
    )
    return answer, {
        "intrinsic_free_dof": intrinsic_free_dof,
        "s0_reflection_count": s0_reflections,
        "cartesian_outside_vertices": sum(
            parameters["mode"] == "cartesian" for parameters in answer
        ),
        "sphere_outside_vertices": sum(
            parameters["mode"] == "sphere" for parameters in answer
        ),
        "fixed_outside_vertices": sum(
            parameters["mode"] == "fixed" for parameters in answer
        ),
    }


def initial_manifold_parameters(
    records: Sequence[dict],
    restarts: int,
    wave: int,
    seedbase: int,
    scale: float,
    seed_standard_controls: bool,
) -> np.ndarray:
    values = np.zeros((len(records), restarts, OUTSIDE, D), dtype=np.float32)
    for local, record in enumerate(records):
        generator = np.random.Generator(
            np.random.PCG64(kernel.graph_seed(seedbase, record["index"], wave))
        )
        for restart in range(restarts):
            for slot, parameters in enumerate(record["manifold_parameters"]):
                if parameters["mode"] == "cartesian":
                    values[local, restart, slot] = generator.normal(0.0, scale, size=D)
                elif parameters["mode"] == "sphere":
                    dimension = parameters["coefficient_dimension"]
                    coefficients = generator.normal(size=dimension)
                    norm = float(np.linalg.norm(coefficients))
                    if norm == 0.0:
                        coefficients[0] = 1.0
                    values[local, restart, slot, :dimension] = coefficients
        if (
            seed_standard_controls
            and wave == 0
            and record["standard_aligned"] is not None
        ):
            outside = record["standard_aligned"][SEED_SIZE:]
            for slot, parameters in enumerate(record["manifold_parameters"]):
                point = outside[slot]
                if parameters["mode"] == "cartesian":
                    values[local, 0, slot] = point
                elif parameters["mode"] == "sphere":
                    coefficients = (
                        parameters["basis"].T
                        @ (point - parameters["centre"])
                        / parameters["radius"]
                    )
                    dimension = parameters["coefficient_dimension"]
                    require(
                        abs(float(np.linalg.norm(coefficients)) - 1.0) < 2e-8,
                        "control is off the seed-link sphere",
                    )
                    values[local, 0, slot, :dimension] = coefficients
                else:
                    require(
                        float(np.linalg.norm(point - parameters["point"])) < 2e-8,
                        "control disagrees with a forced manifold point",
                    )
    return values


def manifold_points_numpy(record: dict, raw: np.ndarray) -> np.ndarray:
    """Accelerator-free evaluator used by controls and focused tests."""

    require(raw.shape[-2:] == (OUTSIDE, D), "raw manifold parameter shape")
    outside = np.empty_like(raw, dtype=np.float64)
    for slot, parameters in enumerate(record["manifold_parameters"]):
        if parameters["mode"] == "cartesian":
            outside[..., slot, :] = raw[..., slot, :]
        elif parameters["mode"] == "fixed":
            outside[..., slot, :] = parameters["point"]
        else:
            dimension = parameters["coefficient_dimension"]
            coefficients = raw[..., slot, :dimension].astype(np.float64)
            norms = np.linalg.norm(coefficients, axis=-1, keepdims=True)
            require(bool(np.all(norms > 1e-14)), "zero sphere coefficient vector")
            unit = coefficients / norms
            direction = np.einsum("ij,...j->...i", parameters["basis"], unit)
            outside[..., slot, :] = (
                parameters["centre"] + parameters["radius"] * direction
            )
    fixed = np.broadcast_to(
        kernel.simplex_coordinates(), raw.shape[:-2] + (SEED_SIZE, D)
    )
    return np.concatenate((fixed, outside), axis=-2)


def optimize_wave_manifold_n19(
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
    """Adam on intrinsic seed-link manifolds, followed by Cartesian LM later."""

    import torch

    require(initialization == "seed_sphere", "manifold kernel requires seed_sphere")
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
    fixed = torch.as_tensor(kernel.simplex_coordinates(), dtype=dtype, device=device)
    fixed = fixed[None, None, :, :].expand(count, restarts, -1, -1)

    mode = np.zeros((count, OUTSIDE), dtype=np.int64)
    centres = np.zeros((count, OUTSIDE, D), dtype=np.float32)
    bases = np.zeros((count, OUTSIDE, D, D), dtype=np.float32)
    coefficient_masks = np.zeros((count, OUTSIDE, D), dtype=np.float32)
    radii = np.zeros((count, OUTSIDE, 1), dtype=np.float32)
    anchors = np.zeros((count, OUTSIDE, D), dtype=np.float32)
    for local, record in enumerate(records):
        for slot, parameters in enumerate(record["manifold_parameters"]):
            if parameters["mode"] == "cartesian":
                mode[local, slot] = 0
            elif parameters["mode"] == "sphere":
                mode[local, slot] = 1
                dimension = parameters["coefficient_dimension"]
                centres[local, slot] = parameters["centre"]
                bases[local, slot, :, :dimension] = parameters["basis"]
                coefficient_masks[local, slot, :dimension] = 1.0
                radii[local, slot, 0] = parameters["radius"]
            else:
                mode[local, slot] = 2
                anchors[local, slot] = parameters["point"]

    initial = initial_manifold_parameters(
        records,
        restarts,
        wave,
        seedbase,
        initialization_scale,
        seed_standard_controls,
    )
    raw = torch.nn.Parameter(torch.as_tensor(initial, dtype=dtype, device=device))
    mode_t = torch.as_tensor(mode, device=device)
    centre_t = torch.as_tensor(centres, dtype=dtype, device=device)[:, None]
    basis_t = torch.as_tensor(bases, dtype=dtype, device=device)
    coefficient_mask_t = torch.as_tensor(
        coefficient_masks, dtype=dtype, device=device
    )[:, None]
    radius_t = torch.as_tensor(radii, dtype=dtype, device=device)[:, None]
    anchor_t = torch.as_tensor(anchors, dtype=dtype, device=device)[:, None]
    cartesian_mask = (mode_t == 0).to(dtype=dtype)[:, None, :, None]
    sphere_mask = (mode_t == 1).to(dtype=dtype)[:, None, :, None]
    fixed_mask = (mode_t == 2).to(dtype=dtype)[:, None, :, None]
    optimizer = torch.optim.Adam([raw], lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=max(steps, 1), eta_min=learning_rate * 0.02
    )
    collision_squared = collision_distance * collision_distance

    def losses() -> tuple[object, object, object, object]:
        coefficients = raw * coefficient_mask_t
        norms = torch.linalg.vector_norm(coefficients, dim=-1, keepdim=True).clamp_min(1e-12)
        unit = coefficients / norms
        directions = torch.einsum("boij,broj->broi", basis_t, unit)
        sphere_points = centre_t + radius_t * directions
        outside = (
            cartesian_mask * raw
            + sphere_mask * sphere_points
            + fixed_mask * anchor_t
        )
        points = torch.cat((fixed, outside), dim=2)
        gram = points @ points.transpose(2, 3)
        norm = torch.diagonal(gram, dim1=2, dim2=3)
        distances = (
            norm[:, :, :, None] + norm[:, :, None, :] - 2.0 * gram
        ).clamp_min(0.0)
        error = distances - 1.0
        edge_loss = (error.square() * adjacency[:, None]).sum(dim=(2, 3)) / (
            2.0 * edge_counts[:, None]
        )
        repulsion = torch.exp(-distances / collision_squared) * upper
        collision_loss = repulsion.sum(dim=(2, 3)) / PAIR_COUNT
        objective = edge_loss + collision_weight * collision_loss
        return points, distances, edge_loss, objective

    initial_controls = {}
    with torch.no_grad():
        _, distances0, edge0, _ = losses()
        for local, record in enumerate(records):
            if seed_standard_controls and record["standard_compatible"] and wave == 0:
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
        torch.nn.utils.clip_grad_norm_([raw], max_norm=10.0)
        optimizer.step()
        # The forward map is radially invariant on each coefficient sphere.
        # Project after every Adam step so irrelevant radial momentum cannot
        # grow without bound; padded coordinates and fixed slots stay zero.
        with torch.no_grad():
            coefficients = raw * coefficient_mask_t
            coefficient_norms = torch.linalg.vector_norm(
                coefficients, dim=-1, keepdim=True
            )
            projected = coefficients / coefficient_norms.clamp_min(1e-12)
            fallback = torch.zeros_like(projected)
            fallback[..., 0] = 1.0
            projected = torch.where(
                (coefficient_norms < 1e-12) & sphere_mask.bool(),
                fallback,
                projected,
            )
            raw.copy_(cartesian_mask * raw + sphere_mask * projected)
        scheduler.step()
    if device_name == "mps":
        torch.mps.synchronize()
    elapsed = time.perf_counter() - started

    with torch.no_grad():
        points, distances, edge_loss, _ = losses()
        error = (distances - 1.0).abs() * adjacency[:, None]
        max_edge = error.amax(dim=(2, 3))
        pair_distances = torch.sqrt(
            torch.where(
                upper.bool(), distances, torch.full_like(distances, float("inf"))
            ).amin(dim=(2, 3))
        )
        points_cpu = points.cpu().numpy()
        edge_rms_cpu = torch.sqrt(edge_loss).cpu().numpy()
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
        "parameterization": "intrinsic_common_sphere_with_forced_S0_reflections",
    }


def dependency_boundary() -> dict:
    paths = {
        "numerical_kernel": Path(kernel.__file__).resolve(),
        "v8_manifest": V8_MANIFEST,
        "v8_verification": V8_VERIFICATION,
        "k7_increment_report": K7_INCREMENT_REPORT,
        "k7_increment_verification": K7_INCREMENT_VERIFICATION,
    }
    hashes = {name: kernel.file_sha256(path) for name, path in paths.items()}
    require(
        hashes["numerical_kernel"] == EXPECTED_KERNEL_SOURCE_SHA256,
        "imported numerical-kernel source hash drift",
    )
    for name, expected in EXPECTED_INPUT_SHA256.items():
        require(hashes[name] == expected, f"{name} hash drift")
    return {
        name: {"path": str(paths[name]), "sha256": hashes[name]}
        for name in sorted(paths)
    }


def validate_rows(value: object, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == N
        and all(type(row) is int for row in value),
        f"{label}: malformed adjacency",
    )
    rows = tuple(value)
    for vertex, row in enumerate(rows):
        require(0 <= row < 1 << N, f"{label}: row outside n=19 mask")
        require(not row & (1 << vertex), f"{label}: loop at {vertex}")
        for other in range(N):
            require(
                bool(row & (1 << other)) == bool(rows[other] & (1 << vertex)),
                f"{label}: asymmetric adjacency",
            )
    return rows


def enumerate_cliques(rows: Sequence[int], target: int) -> list[tuple[int, ...]]:
    answer: list[tuple[int, ...]] = []

    def search(chosen: tuple[int, ...], candidates: int) -> None:
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


def edge_count(rows: Sequence[int]) -> int:
    return sum(row.bit_count() for row in rows) // 2


def intrinsic_profile(rows: Sequence[int], seed: Sequence[int]) -> dict:
    seed_set = set(seed)
    free_dof = 0
    s0 = 0
    for vertex in range(N):
        if vertex in seed_set:
            continue
        neighbours = sum(bool(rows[vertex] & (1 << fixed)) for fixed in seed)
        require(neighbours < 7, "K8 entered n=19 manifold gauge")
        free_dof += D if neighbours == 0 else D - neighbours
        s0 += neighbours == 6
    return {"intrinsic_free_dof": free_dof, "s0_reflection_count": s0}


def load_campaign() -> tuple[list[dict], dict]:
    """Load all 261 physical records and tag the exact-negative calibration."""

    manifest = json.loads(V8_MANIFEST.read_text(encoding="utf-8"))
    verification = json.loads(V8_VERIFICATION.read_text(encoding="utf-8"))
    increment = json.loads(K7_INCREMENT_REPORT.read_text(encoding="utf-8"))
    increment_verification = json.loads(
        K7_INCREMENT_VERIFICATION.read_text(encoding="utf-8")
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v8",
        "v8 manifest schema drift",
    )
    require(verification.get("status") == "PASS", "v8 verification is not PASS")
    require(
        verification.get("counts", {}).get("combined") == 261,
        "v8 verified count drift",
    )
    require(
        increment.get("status") == "COMPLETE_EXACT_REJECTION",
        "K7 increment is not a complete exact rejection",
    )
    require(
        increment_verification.get("status") == "PASS"
        and increment_verification.get("conclusion", {}).get("rejected_indices")
        == [EXACTLY_REMOVED_INDEX],
        "K7 increment verification boundary drift",
    )
    require(
        increment.get("summary", {}).get("rejected_indices")
        == [EXACTLY_REMOVED_INDEX],
        "K7 increment rejected-index drift",
    )

    campaign = []
    for source_class in ("K7", "K6_only"):
        source = manifest["classes"][source_class]
        require(len(source["graphs"]) == source["count"], f"{source_class} count")
        for raw in source["graphs"]:
            candidate_index = int(raw["index"])
            rows = validate_rows(raw["adjacency"], str(candidate_index))
            campaign_role = (
                "exact_negative_calibration"
                if candidate_index == EXACTLY_REMOVED_INDEX
                else "active_unresolved"
            )
            if campaign_role == "exact_negative_calibration":
                require(source_class == "K7", "exact removal moved class")
            campaign.append(
                {
                    "candidate_index": candidate_index,
                    "source_class": source_class,
                    "campaign_role": campaign_role,
                    "class_id": f"{source_class}:{candidate_index}",
                    "adjacency": list(rows),
                    "edges": edge_count(rows),
                    "control_coordinates": None,
                }
            )

    physical_class_counts = dict(Counter(item["source_class"] for item in campaign))
    active = [item for item in campaign if item["campaign_role"] == "active_unresolved"]
    active_class_counts = dict(Counter(item["source_class"] for item in active))
    require(
        len(campaign) == EXPECTED_PHYSICAL_GRAPH_COUNT,
        "physical campaign graph-count drift",
    )
    require(len(active) == EXPECTED_ACTIVE_GRAPH_COUNT, "active graph-count drift")
    require(
        physical_class_counts == EXPECTED_PHYSICAL_CLASS_COUNTS,
        "physical class-count drift",
    )
    require(
        active_class_counts == EXPECTED_ACTIVE_CLASS_COUNTS,
        "active class-count drift",
    )
    require(
        kernel.stable_hash([item["candidate_index"] for item in active])
        == EXPECTED_CAMPAIGN_INDICES_SHA256,
        "active campaign index root drift",
    )
    require(
        kernel.stable_hash([item["candidate_index"] for item in campaign])
        == EXPECTED_PHYSICAL_INDICES_SHA256,
        "physical campaign index root drift",
    )
    tagged_graphs = [
        {
            "source_class": item["source_class"],
            "campaign_role": item["campaign_role"],
            "index": item["candidate_index"],
            "adjacency": item["adjacency"],
        }
        for item in campaign
    ]
    require(
        kernel.stable_hash(tagged_graphs) == EXPECTED_PHYSICAL_TAGGED_GRAPHS_SHA256,
        "physical tagged-graph root drift",
    )
    require(
        [item["candidate_index"] for item in campaign if item["campaign_role"] != "active_unresolved"]
        == [EXACTLY_REMOVED_INDEX],
        "exact-negative calibration tagging drift",
    )
    provenance = {
        "input_construction": (
            "all verified v8 records retained physically; verified exact K7 "
            "rejection 3936435 tagged as a non-active negative calibration"
        ),
        "physical_graph_count": len(campaign),
        "active_unresolved_graph_count": len(active),
        "physical_class_counts": physical_class_counts,
        "active_class_counts": active_class_counts,
        "ordered_candidate_indices_sha256": EXPECTED_CAMPAIGN_INDICES_SHA256,
        "ordered_physical_indices_sha256": EXPECTED_PHYSICAL_INDICES_SHA256,
        "physical_tagged_graphs_sha256": EXPECTED_PHYSICAL_TAGGED_GRAPHS_SHA256,
        "exact_negative_calibration_indices": [EXACTLY_REMOVED_INDEX],
        "prospective_source_only_rejections_excluded": (
            "none; all source-only prospective rejections remain in this campaign"
        ),
        "candidate_nonedges": "unconstrained and may also have distance one",
    }
    return campaign, provenance


def gauge_metadata(campaign: Sequence[dict]) -> list[dict]:
    """Enumerate all K7 gauges in K7 and all K6 gauges in K6-only cases."""

    answer = []
    seen_variant_ids: set[int] = set()
    for item in campaign:
        target = 7 if item["source_class"] == "K7" else 6
        rows = validate_rows(item["adjacency"], item["class_id"])
        seeds = enumerate_cliques(rows, target)
        require(bool(seeds), f"{item['class_id']}: missing K{target} gauge")
        require(len(seeds) < VARIANT_STRIDE, "variant stride exhausted")
        if item["source_class"] == "K6_only":
            require(
                not enumerate_cliques(rows, 7),
                f"{item['class_id']}: K6-only graph contains a K7",
            )
        for gauge_ordinal, seed in enumerate(seeds):
            variant_id = item["candidate_index"] * VARIANT_STRIDE + gauge_ordinal
            require(variant_id not in seen_variant_ids, "variant-id collision")
            seen_variant_ids.add(variant_id)
            profile = intrinsic_profile(rows, seed)
            answer.append(
                {
                    "variant_id": variant_id,
                    "candidate_index": item["candidate_index"],
                    "source_class": item["source_class"],
                    "campaign_role": item["campaign_role"],
                    "gauge_ordinal": gauge_ordinal,
                    "gauge_size": target,
                    "seed": list(seed),
                    **profile,
                }
            )
    return answer


def audit_full_gauge_manifest(metadata: Sequence[dict]) -> dict:
    populations = dict(
        Counter(
            (
                "K7_exact_negative_calibration"
                if item["campaign_role"] == "exact_negative_calibration"
                else f"{item['source_class']}_active"
            )
            for item in metadata
        )
    )
    require(len(metadata) == EXPECTED_GAUGE_COUNT, "gauge-count drift")
    require(populations == EXPECTED_GAUGE_POPULATIONS, "gauge populations drift")
    require(
        kernel.stable_hash(list(metadata)) == EXPECTED_GAUGE_MANIFEST_SHA256,
        "gauge-manifest root drift",
    )
    require(
        kernel.stable_hash([item["variant_id"] for item in metadata])
        == EXPECTED_VARIANT_IDS_SHA256,
        "variant-id root drift",
    )
    census = []
    for candidate_index in dict.fromkeys(item["candidate_index"] for item in metadata):
        selected = [item for item in metadata if item["candidate_index"] == candidate_index]
        census.append(
            {
                "candidate_index": candidate_index,
                "source_class": selected[0]["source_class"],
                "campaign_role": selected[0]["campaign_role"],
                "gauge_size": selected[0]["gauge_size"],
                "gauge_count": len(selected),
            }
        )
    require(
        kernel.stable_hash(census) == EXPECTED_GAUGE_CENSUS_SHA256,
        "per-graph gauge census drift",
    )
    dof_profile = [
        {
            "candidate_index": item["candidate_index"],
            "gauge_ordinal": item["gauge_ordinal"],
            "intrinsic_free_dof": item["intrinsic_free_dof"],
            "s0_vertices": item["s0_reflection_count"],
        }
        for item in metadata
    ]
    require(
        kernel.stable_hash(dof_profile) == EXPECTED_INTRINSIC_DOF_PROFILE_SHA256,
        "intrinsic-DOF profile drift",
    )
    active_k7_dof = [
        item["intrinsic_free_dof"]
        for item in metadata
        if item["source_class"] == "K7" and item["campaign_role"] == "active_unresolved"
    ]
    active_k6_dof = [
        item["intrinsic_free_dof"]
        for item in metadata
        if item["source_class"] == "K6_only"
    ]
    require((min(active_k7_dof), max(active_k7_dof)) == (25, 31), "K7 DOF range")
    require((min(active_k6_dof), max(active_k6_dof)) == (27, 36), "K6 DOF range")
    require(
        sum(bool(item["s0_reflection_count"]) for item in metadata) == 6,
        "S0 gauge census drift",
    )
    return {
        "gauge_count": len(metadata),
        "population_counts": populations,
        "gauge_manifest_sha256": EXPECTED_GAUGE_MANIFEST_SHA256,
        "variant_ids_sha256": EXPECTED_VARIANT_IDS_SHA256,
        "per_graph_census_sha256": EXPECTED_GAUGE_CENSUS_SHA256,
        "intrinsic_dof_profile_sha256": EXPECTED_INTRINSIC_DOF_PROFILE_SHA256,
        "active_intrinsic_dof_ranges": {"K7": [25, 31], "K6_only": [27, 36]},
        "gauges_with_forced_S0_reflections": 6,
    }


def align_control(
    coordinates: np.ndarray, new_to_old: Sequence[int], gauge_size: int
) -> np.ndarray:
    ordered = coordinates[list(new_to_old)]
    source = ordered[:gauge_size]
    target = kernel.full_simplex_coordinates()[:gauge_size]
    source_centre = source.mean(axis=0)
    target_centre = target.mean(axis=0)
    covariance = (source - source_centre).T @ (target - target_centre)
    left, _, right_t = np.linalg.svd(covariance, full_matrices=True)
    rotation = left @ right_t
    aligned = (ordered - source_centre) @ rotation + target_centre
    require(
        np.max(np.abs(aligned[:gauge_size] - target)) < 2e-10,
        "positive-control gauge alignment failed",
    )
    return aligned


def reorder_record(item: dict, info: dict) -> dict:
    rows = validate_rows(item["adjacency"], item["class_id"])
    seed = tuple(info["seed"])
    require(len(seed) == info["gauge_size"], "gauge-size mismatch")
    for position, first in enumerate(seed):
        for second in seed[:position]:
            require(bool(rows[first] & (1 << second)), "gauge seed is not a clique")
    remaining = tuple(vertex for vertex in range(N) if vertex not in seed)
    new_to_old = seed + remaining
    old_to_new = {old: new for new, old in enumerate(new_to_old)}
    adjacency = np.zeros((N, N), dtype=np.float32)
    for old_first in range(N):
        for old_second in range(old_first):
            if rows[old_first] & (1 << old_second):
                first = old_to_new[old_first]
                second = old_to_new[old_second]
                adjacency[first, second] = adjacency[second, first] = 1.0
    aligned = None
    if item.get("control_coordinates") is not None:
        aligned = align_control(
            np.asarray(item["control_coordinates"], dtype=np.float64),
            new_to_old,
            info["gauge_size"],
        )
    parameters, parameter_audit = manifold_parameters(adjacency, info["gauge_size"])
    if "intrinsic_free_dof" in info:
        require(
            info["intrinsic_free_dof"] == parameter_audit["intrinsic_free_dof"]
            and info["s0_reflection_count"] == parameter_audit["s0_reflection_count"],
            "combinatorial/manifold intrinsic profile mismatch",
        )
    return {
        "index": info["variant_id"],
        "class_id": item["class_id"],
        "campaign_role": item.get("campaign_role", "positive_control"),
        "edges": item["edges"],
        "standard_compatible": aligned is not None,
        "gauge_size": info["gauge_size"],
        "seed_old_vertices": list(seed),
        "new_to_old": list(new_to_old),
        "adjacency": adjacency,
        "standard_aligned": aligned,
        # Kept for the inherited Cartesian LM/result layer; Adam uses only the
        # intrinsic structure below.
        "sphere_parameters": parameters,
        "manifold_parameters": parameters,
        **parameter_audit,
    }


def make_variants(items: Sequence[dict]) -> tuple[list[dict], list[dict]]:
    metadata = gauge_metadata(items)
    by_index = {item["candidate_index"]: item for item in items}
    variants = [
        reorder_record(by_index[info["candidate_index"]], info) for info in metadata
    ]
    return variants, metadata


def rows_from_coordinates(coordinates: np.ndarray) -> list[int]:
    require(coordinates.shape == (N, D), "positive-control coordinate shape")
    rows = [0] * N
    for first in range(N):
        for second in range(first):
            squared = float(np.sum((coordinates[first] - coordinates[second]) ** 2))
            if abs(squared - 1.0) <= 1e-10:
                rows[first] |= 1 << second
                rows[second] |= 1 << first
    return rows


def positive_control_variants() -> tuple[list[dict], list[dict]]:
    """Build exact K6 and K7 seeded controls on the final n=19 kernel."""

    standard18 = kernel.standard_coordinates()
    standard18_rows = rows_from_coordinates(
        np.vstack((standard18, np.asarray([[0, 0, 0, 0, 0, 3.0]])))
    )
    # Remove the temporary isolated row from the exact standard-18 unit graph.
    standard18_rows = [row & ((1 << 18) - 1) for row in standard18_rows[:18]] + [0]
    standard_seeds = enumerate_cliques(standard18_rows, 6)
    require(len(standard_seeds) == 32, "standard18 K6 census drift")

    tether = standard18[0].copy()
    tether[5] += 1.0
    standard = np.vstack((standard18, tether))
    standard_rows = list(standard18_rows)
    standard_rows[0] |= 1 << 18
    standard_rows[18] |= 1 << 0
    require(
        abs(float(np.sum((standard[18] - standard[0]) ** 2)) - 1.0) < 1e-12,
        "light tether control edge drift",
    )
    standard_item = {
        "candidate_index": -1,
        "source_class": "positive_control",
        "campaign_role": "positive_control",
        "class_id": "positive:standard18_plus_light_tether",
        "adjacency": standard_rows,
        "edges": edge_count(standard_rows),
        "control_coordinates": standard.tolist(),
    }

    k6_seed = standard_seeds[0]
    seed_points = standard18[list(k6_seed)]
    centre = seed_points.mean(axis=0)
    differences = seed_points[1:] - seed_points[0]
    _, _, right_t = np.linalg.svd(differences, full_matrices=True)
    normal = right_t[-1]
    radius = math.sqrt(7.0 / 12.0)
    apex_options = (centre + radius * normal, centre - radius * normal)
    apex = max(
        apex_options,
        key=lambda point: min(float(np.linalg.norm(point - old)) for old in standard18),
    )
    require(
        min(float(np.linalg.norm(apex - old)) for old in standard18) > 1e-6,
        "reflected-apex control is not distinct",
    )
    k7_coordinates = np.vstack((standard18, apex))
    k7_rows = list(standard18_rows)
    for vertex in k6_seed:
        require(
            abs(float(np.sum((apex - standard18[vertex]) ** 2)) - 1.0) < 2e-12,
            "reflected-apex seed edge drift",
        )
        k7_rows[vertex] |= 1 << 18
        k7_rows[18] |= 1 << vertex
    k7_item = {
        "candidate_index": -2,
        "source_class": "positive_control",
        "campaign_role": "positive_control",
        "class_id": "positive:standard18_plus_reflected_K7_apex",
        "adjacency": k7_rows,
        "edges": edge_count(k7_rows),
        "control_coordinates": k7_coordinates.tolist(),
    }

    metadata = [
        {
            "variant_id": CONTROL_VARIANT_BASE,
            "candidate_index": -1,
            "source_class": "positive_control",
            "campaign_role": "positive_control",
            "gauge_ordinal": 0,
            "gauge_size": 6,
            "seed": list(standard_seeds[0]),
            "control_name": "standard18_plus_light_tether_first_K6",
        },
        {
            "variant_id": CONTROL_VARIANT_BASE + 1,
            "candidate_index": -1,
            "source_class": "positive_control",
            "campaign_role": "positive_control",
            "gauge_ordinal": 31,
            "gauge_size": 6,
            "seed": list(standard_seeds[-1]),
            "control_name": "standard18_plus_light_tether_last_K6",
        },
        {
            "variant_id": CONTROL_VARIANT_BASE + 2,
            "candidate_index": -2,
            "source_class": "positive_control",
            "campaign_role": "positive_control",
            "gauge_ordinal": 0,
            "gauge_size": 7,
            "seed": list(k6_seed) + [18],
            "control_name": "standard18_plus_distinct_reflected_K7_apex",
        },
    ]
    by_index = {-1: standard_item, -2: k7_item}
    variants = [reorder_record(by_index[info["candidate_index"]], info) for info in metadata]
    return variants, metadata


def process_batch(
    records: Sequence[dict],
    config: dict,
    executor: concurrent.futures.Executor | None,
) -> dict:
    with n19_kernel_context():
        return kernel.process_batch(records, config, executor)


def decorate_result(result: dict, metadata: Sequence[dict]) -> dict:
    lookup = {item["variant_id"]: item for item in metadata}
    require(len(lookup) == len(metadata), "duplicate metadata variant id")

    def decorate(items: Sequence[dict]) -> list[dict]:
        answer = []
        for original in items:
            item = dict(original)
            variant_id = int(item.pop("index"))
            require(variant_id in lookup, "kernel returned unknown variant")
            info = lookup[variant_id]
            item.update(
                {
                    "variant_id": variant_id,
                    "candidate_index": info["candidate_index"],
                    "source_class": info["source_class"],
                    "campaign_role": info["campaign_role"],
                    "gauge_ordinal": info["gauge_ordinal"],
                    "gauge_size": info["gauge_size"],
                    "gauge_seed_old_vertices": info["seed"],
                }
            )
            if "control_name" in info:
                item["control_name"] = info["control_name"]
            answer.append(item)
        return answer

    decorated = {
        key: decorate(result[key])
        for key in ("rows", "witnesses", "retained_endpoints", "lm_attempts")
    }
    decorated["wave_profiles"] = result["wave_profiles"]
    decorated["lm_profile"] = result["lm_profile"]
    return decorated


def select_campaign(campaign: Sequence[dict], indices: Sequence[int]) -> list[dict]:
    if not indices:
        return list(campaign)
    requested = set(map(int, indices))
    available = {item["candidate_index"] for item in campaign}
    require(not requested - available, f"unknown candidate indices: {sorted(requested - available)}")
    return [item for item in campaign if item["candidate_index"] in requested]


def configuration(args: argparse.Namespace, selected: Sequence[dict], source_hash: str) -> dict:
    return {
        "schema": "d6-19-gpu-mps-residue-v1-config-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS_NO_REALIZABILITY_CONCLUSIONS",
        "orchestrator_source_sha256": source_hash,
        "numerical_kernel_source_sha256": EXPECTED_KERNEL_SOURCE_SHA256,
        "input_sha256": EXPECTED_INPUT_SHA256,
        "full_active_indices_sha256": EXPECTED_CAMPAIGN_INDICES_SHA256,
        "full_physical_indices_sha256": EXPECTED_PHYSICAL_INDICES_SHA256,
        "full_gauge_manifest_sha256": EXPECTED_GAUGE_MANIFEST_SHA256,
        "selected_candidate_indices": [item["candidate_index"] for item in selected],
        "selected_candidate_indices_sha256": kernel.stable_hash(
            [item["candidate_index"] for item in selected]
        ),
        "dimension": D,
        "vertices": N,
        "gauge_policy": (
            "all unit K7 seeds for K7-class graphs; all unit K6 seeds for "
            "K6-only graphs; lexicographic within candidate"
        ),
        "variant_id_formula": "4096 * published_candidate_index + gauge_ordinal",
        "total_restarts": args.total_restarts,
        "restarts_per_wave": args.restarts_per_wave,
        "optimizer": "PyTorch Adam with cosine learning-rate decay",
        "primary_parameterization": (
            "intrinsic common-unit-sphere coordinates; all required links to "
            "the fixed K6/K7 seed are exact by construction; K7 S0 uses the "
            "unique distinct reflected root"
        ),
        "optimizer_steps": args.steps,
        "learning_rate": args.learning_rate,
        "initialization_method": "seed_sphere",
        "initialization_scale": args.initialization_scale,
        "seed_standard_controls": True,
        "positive_controls": (
            "two gauges of exact standard18+light-tether and one exact "
            "standard18+distinct-reflected-K7-apex gauge, once per configuration"
        ),
        "seedbase": args.seedbase,
        "graph_seed_formula": (
            "(seedbase*1000003 + variant_id*7919 + wave*104729) mod 2^64"
        ),
        "collision_distance": args.collision_distance,
        "collision_weight": args.collision_weight,
        "collision_semantics": (
            "all-pair short-distance repulsion is a heuristic distinctness "
            "bias, not a nonedge distance constraint"
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
        "lm_worker_thread_environment": LM_THREAD_ENVIRONMENT,
        "device": args.device,
        "dtype": "MPS/CPU float32 Adam and SciPy float64 LM refinement",
        "candidate_nonedges": "absent from edge residual; unconstrained and may be unit",
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }


def checkpoint_path(run_directory: Path, start: int, end: int) -> Path:
    return run_directory / f"classes_{start:04d}_{end:04d}.json.gz"


def load_checkpoint(
    path: Path,
    config_hash: str,
    start: int,
    end: int,
    candidate_indices: Sequence[int],
) -> dict:
    with gzip.open(path, "rt", encoding="ascii") as stream:
        payload = json.load(stream)
    claimed = payload.pop("checkpoint_sha256", None)
    require(claimed == kernel.stable_hash(payload), f"{path}: payload hash")
    payload["checkpoint_sha256"] = claimed
    require(
        payload.get("schema") == "d6-19-gpu-mps-residue-v1-checkpoint-v1",
        f"{path}: schema",
    )
    require(payload.get("config_sha256") == config_hash, f"{path}: config")
    require(payload.get("class_range") == [start, end], f"{path}: range")
    require(payload.get("candidate_indices") == list(candidate_indices), f"{path}: indices")
    rows = payload.get("rows", [])
    metadata = payload.get("gauge_metadata", [])
    require(len(rows) == len(metadata) == payload.get("gauge_count"), f"{path}: gauges")
    require(
        [row["variant_id"] for row in rows]
        == [item["variant_id"] for item in metadata],
        f"{path}: variant order",
    )
    return payload


def control_checkpoint_path(run_directory: Path) -> Path:
    return run_directory / "positive_controls.json.gz"


def guard_output_target(path: Path, config_hash: str, *, compressed: bool) -> None:
    """Refuse to overwrite a completed artifact from another configuration."""

    if not path.exists():
        return
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="ascii") as stream:
        payload = json.load(stream)
    require(
        payload.get("config_sha256") == config_hash,
        f"refusing to overwrite output from another configuration: {path}",
    )


def audit_positive_controls(payload: dict, config: dict) -> dict:
    rows = payload["rows"]
    require(len(rows) == 3, "positive-control result count")
    exact = [
        row
        for row in rows
        if math.isfinite(row["initial_standard_control_edge_rms"])
        and row["initial_standard_control_edge_rms"] <= 2e-6
        and row["initial_standard_control_minimum_distance"]
        >= config["distinctness_threshold"]
    ]
    final_recoveries = sum(bool(row["lm_candidate"]) for row in rows)
    final_recovery_required = bool(
        config["lm_top_k"] > 0 and config["cpu_refine_max_nfev"] > 0
    )
    if len(exact) != len(rows):
        status = "FAIL_EXACT_SEEDED_INPUT"
    elif final_recovery_required and final_recoveries != len(rows):
        status = "FAIL_FINAL_KERNEL_RECOVERY"
    elif final_recovery_required:
        status = "PASS"
    else:
        status = "PASS_SEEDED_INPUT_ONLY"
    return {
        "status": status,
        "variants": len(rows),
        "exact_seeded_input_passes": len(exact),
        "final_recovery_required": final_recovery_required,
        "lm_candidate_recoveries": final_recoveries,
        "best_initial_edge_rms": min(
            row["initial_standard_control_edge_rms"] for row in rows
        ),
        "worst_initial_edge_rms": max(
            row["initial_standard_control_edge_rms"] for row in rows
        ),
        "minimum_initial_distance": min(
            row["initial_standard_control_minimum_distance"] for row in rows
        ),
        "mathematical_conclusion": False,
    }


def ensure_positive_controls(
    run_directory: Path,
    config: dict,
    config_hash: str,
    executor: concurrent.futures.Executor | None,
) -> tuple[dict, dict]:
    path = control_checkpoint_path(run_directory)
    if path.exists():
        with gzip.open(path, "rt", encoding="ascii") as stream:
            payload = json.load(stream)
        claimed = payload.pop("checkpoint_sha256", None)
        require(claimed == kernel.stable_hash(payload), "positive-control checkpoint hash")
        payload["checkpoint_sha256"] = claimed
        require(
            payload.get("schema") == "d6-19-gpu-mps-residue-v1-controls-v1"
            and payload.get("config_sha256") == config_hash,
            "positive-control checkpoint boundary mismatch",
        )
    else:
        variants, metadata = positive_control_variants()
        started = time.perf_counter()
        decorated = decorate_result(process_batch(variants, config, executor), metadata)
        payload = {
            "schema": "d6-19-gpu-mps-residue-v1-controls-v1",
            "proof_status": "NUMERICAL_POSITIVE_CONTROLS_NO_THEOREM_CLAIM",
            "config_sha256": config_hash,
            "gauge_metadata": metadata,
            "wall_seconds": time.perf_counter() - started,
            **decorated,
            "mathematical_rejections": 0,
            "mathematical_realizability_conclusions": 0,
        }
        payload["checkpoint_sha256"] = kernel.stable_hash(payload)
        kernel.atomic_gzip_json(path, payload)
    audit = audit_positive_controls(payload, config)
    require(
        audit["status"] in ("PASS", "PASS_SEEDED_INPUT_ONLY"),
        "positive controls failed; campaign blocked",
    )
    return payload, audit


def set_lm_thread_environment() -> None:
    """Set limits early enough for fresh spawn workers to inherit them."""

    # Spawned workers inherit this environment before importing NumPy/SciPy,
    # unlike limits set only in the initializer after module import.
    for name, value in LM_THREAD_ENVIRONMENT.items():
        os.environ[name] = value


def make_executor(config: dict) -> concurrent.futures.ProcessPoolExecutor | None:
    if config["cpu_workers"] <= 1 or config["lm_top_k"] <= 0:
        return None
    set_lm_thread_environment()
    return concurrent.futures.ProcessPoolExecutor(
        max_workers=config["cpu_workers"],
        mp_context=multiprocessing.get_context("spawn"),
        initializer=limit_lm_worker_threads,
    )


def limit_lm_worker_threads() -> None:
    """Keep each of the ten LM processes to one BLAS/OpenMP thread."""

    for name, value in LM_THREAD_ENVIRONMENT.items():
        os.environ[name] = value
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=1)
    except ImportError:
        # The environment variables still constrain libraries initialized by
        # SciPy after this spawn initializer.
        pass


def summarize_candidates(rows: Sequence[dict]) -> list[dict]:
    grouped: dict[int, list[dict]] = {}
    for row in rows:
        grouped.setdefault(int(row["candidate_index"]), []).append(row)
    answer = []
    for candidate_index in sorted(grouped):
        selected = sorted(grouped[candidate_index], key=lambda row: row["gauge_ordinal"])
        best = min(
            selected,
            key=lambda row: (
                float(row["best_distinct_lm_edge_rms"]),
                float(row["best_distinct_edge_rms"]),
                row["gauge_ordinal"],
            ),
        )
        candidates = [row for row in selected if row["lm_candidate"] or row["gpu_candidate"]]
        answer.append(
            {
                "candidate_index": candidate_index,
                "source_class": selected[0]["source_class"],
                "campaign_role": selected[0]["campaign_role"],
                "edges": selected[0]["edges"],
                "gauge_size": selected[0]["gauge_size"],
                "gauge_count": len(selected),
                "optimizer_starts": sum(row["optimizer_starts"] for row in selected),
                "cpu_lm_attempts": sum(row["cpu_lm_attempts"] for row in selected),
                "candidate_gauges": len(candidates),
                "heuristic_status": (
                    "NUMERICAL_CANDIDATE_REQUIRES_CERTIFICATION"
                    if candidates
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


def assemble(
    checkpoints: Sequence[dict],
    control_payload: dict,
    control_audit: dict,
    config: dict,
    config_hash: str,
    provenance: dict,
    dependencies: dict,
    args: argparse.Namespace,
    state: dict,
    invocation_wall_seconds: float,
) -> dict:
    guard_output_target(args.details, config_hash, compressed=True)
    guard_output_target(args.report, config_hash, compressed=False)
    keys = ("rows", "witnesses", "retained_endpoints", "lm_attempts")
    merged = {
        key: [item for checkpoint in checkpoints for item in checkpoint[key]]
        for key in keys
    }
    merged["rows"].sort(key=lambda row: (row["candidate_index"], row["gauge_ordinal"]))
    merged["witnesses"].sort(key=lambda row: (row["candidate_index"], row["gauge_ordinal"]))
    merged["retained_endpoints"].sort(
        key=lambda row: (row["candidate_index"], row["gauge_ordinal"])
    )
    merged["lm_attempts"].sort(key=lambda row: (row["candidate_index"], row["gauge_ordinal"]))
    summaries = summarize_candidates(merged["rows"])
    active_summaries = [
        summary for summary in summaries if summary["campaign_role"] == "active_unresolved"
    ]
    negative = [
        summary
        for summary in summaries
        if summary["campaign_role"] == "exact_negative_calibration"
    ]
    require(len(negative) in (0, 1), "exact-negative calibration summary count")
    negative_audit = {
        "included": bool(negative),
        "candidate_index": EXACTLY_REMOVED_INDEX,
        "status": (
            "NUMERICAL_CANDIDATE_ALERT_EXACT_REJECTION_STILL_CONTROLS"
            if negative and negative[0]["candidate_gauges"]
            else "PASS_NO_NUMERICAL_CANDIDATE"
            if negative
            else "NOT_IN_SELECTED_SUBSET"
        ),
        "candidate_gauges": negative[0]["candidate_gauges"] if negative else 0,
        "mathematical_conclusion": False,
    }
    details = {
        "schema": "d6-19-gpu-mps-residue-v1-details-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS_NO_REALIZABILITY_CONCLUSIONS",
        "config_sha256": config_hash,
        **merged,
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }
    kernel.atomic_gzip_json(args.details, details)
    report = {
        "schema": "d6-19-gpu-mps-residue-v1-report-v1",
        "status": (
            "COMPLETE_HEURISTIC_NEGATIVE_CALIBRATION_ALERT"
            if negative_audit["status"].startswith("NUMERICAL_CANDIDATE_ALERT")
            else "COMPLETE_HEURISTIC_ONLY"
        ),
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS_NO_REALIZABILITY_CONCLUSIONS",
        "warning": (
            "Failure in every gauge proves nothing. Any numerical candidate "
            "requires exact or rigorous interval certification."
        ),
        "config": config,
        "config_sha256": config_hash,
        "producer_source": {
            "path": str(Path(__file__).resolve()),
            "sha256": kernel.file_sha256(Path(__file__).resolve()),
        },
        "dependency_boundary": dependencies,
        "input_provenance": provenance,
        "git_head_at_launch": state["git_head_at_launch"],
        "launch_command": state["launch_command"],
        "assembly_command": shlex.join(sys.argv),
        "environment": state["environment"],
        "invocation_wall_seconds": invocation_wall_seconds,
        "checkpoint_compute_wall_seconds": sum(
            checkpoint["wall_seconds"] for checkpoint in checkpoints
        ),
        "physical_candidate_count": len(summaries),
        "active_unresolved_candidate_count": len(active_summaries),
        "gauge_count": len(merged["rows"]),
        "optimizer_starts": sum(row["optimizer_starts"] for row in merged["rows"]),
        "cpu_lm_attempts": sum(row["cpu_lm_attempts"] for row in merged["rows"]),
        "numerical_candidate_count": sum(
            summary["heuristic_status"].startswith("NUMERICAL_CANDIDATE")
            for summary in active_summaries
        ),
        "candidate_summaries": summaries,
        "positive_control_audit": control_audit,
        "exact_negative_calibration_audit": negative_audit,
        "positive_control_checkpoint_sha256": control_payload["checkpoint_sha256"],
        "optional_nonedge_semantics_audit": {
            "status": "PASS_BY_PINNED_OBJECTIVE_CONSTRUCTION_AND_FOCUSED_TEST",
            "unit_fit_mask": "adjacency entries equal to one only",
            "candidate_nonedges": "not fitted; unconstrained and may also be unit",
            "all_pair_term": "soft short-distance collision bias only",
        },
        "details": {
            "path": str(args.details),
            "sha256": kernel.file_sha256(args.details),
            "rows": len(merged["rows"]),
            "witnesses_with_coordinates": len(merged["witnesses"]),
            "retained_endpoints_with_coordinates": len(merged["retained_endpoints"]),
        },
        "checkpoints": {
            "count": len(checkpoints),
            "ordered_sha256": kernel.stable_hash(
                [checkpoint["checkpoint_sha256"] for checkpoint in checkpoints]
            ),
        },
        "trust": {
            "floating_point": "PyTorch/MPS float32 and SciPy float64; no rigorous rounding",
            "optimizer": "non-recovery has no mathematical meaning",
            "collision_bias": "discovery aid that can alter numerical basins",
            "nonedges": "unconstrained and may also be unit",
        },
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }
    kernel.atomic_json(args.report, report)
    return report


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


def run(
    selected: Sequence[dict],
    config: dict,
    config_hash: str,
    provenance: dict,
    dependencies: dict,
    args: argparse.Namespace,
) -> dict | None:
    if config["cpu_workers"] > 1 and config["lm_top_k"] > 0:
        set_lm_thread_environment()
    run_directory = args.checkpoint_root / config_hash[:16]
    run_directory.mkdir(parents=True, exist_ok=True)
    state_path = run_directory / "run_state.json"
    expected_indices = [item["candidate_index"] for item in selected]
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        require(
            state.get("schema") == "d6-19-gpu-mps-residue-v1-run-state-v1"
            and state.get("config_sha256") == config_hash
            and state.get("candidate_indices") == expected_indices,
            "existing run-state boundary mismatch",
        )
    else:
        state = {
            "schema": "d6-19-gpu-mps-residue-v1-run-state-v1",
            "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
            "config_sha256": config_hash,
            "physical_candidate_count": len(selected),
            "active_unresolved_candidate_count": sum(
                item["campaign_role"] == "active_unresolved" for item in selected
            ),
            "candidate_indices": expected_indices,
            "checkpoint_directory": str(run_directory),
            "started_at_utc": kernel.utc_now(),
            "pid_at_launch": os.getpid(),
            "git_head_at_launch": git_head(),
            "launch_command": shlex.join(sys.argv),
            "source_sha256": config["orchestrator_source_sha256"],
            "environment": {
                "python": sys.version,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "cpu_count": os.cpu_count(),
                "lm_worker_thread_environment": {
                    name: os.environ.get(name) for name in LM_THREAD_ENVIRONMENT
                },
            },
            "completed_class_chunks": 0,
            "invocations": [],
            "mathematical_rejections": 0,
            "mathematical_realizability_conclusions": 0,
        }
    state.setdefault("invocations", []).append(
        {"started_at_utc": kernel.utc_now(), "git_head": git_head(), "command": shlex.join(sys.argv)}
    )
    kernel.atomic_json(state_path, state)

    checkpoints = []
    launched = 0
    started = time.perf_counter()
    executor = make_executor(config)
    try:
        control_payload, control_audit = ensure_positive_controls(
            run_directory, config, config_hash, executor
        )
        print(
            f"positive controls {control_audit['status']} "
            f"exact={control_audit['exact_seeded_input_passes']}/3",
            flush=True,
        )
        if args.control_only:
            return {
                "status": "CONTROL_ONLY_COMPLETE",
                "config_sha256": config_hash,
                "positive_control_audit": control_audit,
                "checkpoint": str(control_checkpoint_path(run_directory)),
                "mathematical_rejections": 0,
                "mathematical_realizability_conclusions": 0,
            }
        for start in range(0, len(selected), config["class_chunk_size"]):
            end = min(start + config["class_chunk_size"], len(selected))
            candidate_indices = [item["candidate_index"] for item in selected[start:end]]
            path = checkpoint_path(run_directory, start, end)
            if path.exists():
                payload = load_checkpoint(path, config_hash, start, end, candidate_indices)
            else:
                if args.max_new_chunks is not None and launched >= args.max_new_chunks:
                    break
                variants, metadata = make_variants(selected[start:end])
                chunk_started = time.perf_counter()
                decorated = decorate_result(process_batch(variants, config, executor), metadata)
                payload = {
                    "schema": "d6-19-gpu-mps-residue-v1-checkpoint-v1",
                    "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
                    "config_sha256": config_hash,
                    "class_range": [start, end],
                    "candidate_indices": candidate_indices,
                    "gauge_count": len(metadata),
                    "gauge_metadata": metadata,
                    "completed_at_utc": kernel.utc_now(),
                    "wall_seconds": time.perf_counter() - chunk_started,
                    **decorated,
                    "mathematical_rejections": 0,
                    "mathematical_realizability_conclusions": 0,
                }
                payload["checkpoint_sha256"] = kernel.stable_hash(payload)
                kernel.atomic_gzip_json(path, payload)
                payload = load_checkpoint(path, config_hash, start, end, candidate_indices)
                launched += 1
                print(
                    f"classes {start}:{end} gauges={payload['gauge_count']} "
                    f"lm_candidates={sum(row['lm_candidate'] for row in payload['rows'])}",
                    flush=True,
                )
            checkpoints.append(payload)
            state["completed_class_chunks"] = len(checkpoints)
            state["completed_candidates"] = sum(
                len(checkpoint["candidate_indices"]) for checkpoint in checkpoints
            )
            state["completed_gauges"] = sum(
                checkpoint["gauge_count"] for checkpoint in checkpoints
            )
            kernel.atomic_json(state_path, state)
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
        control_payload,
        control_audit,
        config,
        config_hash,
        provenance,
        dependencies,
        args,
        state,
        time.perf_counter() - started,
    )
    state["status"] = report["status"]
    state["report"] = {"path": str(args.report), "sha256": kernel.file_sha256(args.report)}
    kernel.atomic_json(state_path, state)
    return report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--selection-only", action="store_true")
    result.add_argument("--control-only", action="store_true")
    result.add_argument(
        "--index",
        dest="indices",
        type=int,
        action="append",
        help="restrict to a published candidate index; repeatable",
    )
    result.add_argument("--device", choices=("mps", "cpu"), default="mps")
    result.add_argument("--total-restarts", type=int, default=16)
    result.add_argument("--restarts-per-wave", type=int, default=16)
    result.add_argument("--steps", type=int, default=1000)
    result.add_argument("--learning-rate", type=float, default=0.03)
    result.add_argument("--initialization-scale", type=float, default=0.7)
    result.add_argument("--seedbase", type=int, default=619_260_7721)
    result.add_argument("--collision-distance", type=float, default=0.12)
    result.add_argument("--collision-weight", type=float, default=0.05)
    result.add_argument("--distinctness-threshold", type=float, default=0.01)
    result.add_argument("--gpu-candidate-rms", type=float, default=3e-4)
    result.add_argument("--cpu-refine-max-nfev", type=int, default=1500)
    result.add_argument("--lm-top-k", type=int, default=2)
    result.add_argument("--lm-candidate-rms", type=float, default=1e-8)
    result.add_argument("--lm-collision-distance", type=float, default=0.05)
    result.add_argument("--lm-collision-weight", type=float, default=100.0)
    result.add_argument("--cpu-workers", type=int, default=10)
    result.add_argument("--class-chunk-size", type=int, default=8)
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
    if args.selection_only and args.control_only:
        parser_.error("selection-only and control-only are mutually exclusive")
    if args.total_restarts < 1 or args.restarts_per_wave < 1 or args.steps < 1:
        parser_.error("restart and optimizer step counts must be positive")
    if args.cpu_workers < 1 or args.class_chunk_size < 1:
        parser_.error("CPU worker and class chunk counts must be positive")
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
    campaign, provenance = load_campaign()
    full_metadata = gauge_metadata(campaign)
    full_audit = audit_full_gauge_manifest(full_metadata)
    selected = select_campaign(campaign, args.indices or ())
    selected_metadata = gauge_metadata(selected)
    provenance["full_gauge_manifest"] = full_audit
    provenance["selected_gauge_count"] = len(selected_metadata)
    provenance["selected_gauge_manifest_sha256"] = kernel.stable_hash(selected_metadata)
    source_hash = kernel.file_sha256(Path(__file__).resolve())
    config = configuration(args, selected, source_hash)
    config_hash = kernel.stable_hash(config)
    if args.selection_only:
        print(
            json.dumps(
                {
                    "schema": "d6-19-gpu-mps-residue-v1-selection-v1",
                    "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
                    "config_sha256": config_hash,
                    "physical_candidates": len(selected),
                    "active_unresolved_candidates": sum(
                        item["campaign_role"] == "active_unresolved" for item in selected
                    ),
                    "exact_negative_calibrations": sum(
                        item["campaign_role"] == "exact_negative_calibration"
                        for item in selected
                    ),
                    "gauges": len(selected_metadata),
                    "optimizer_starts": len(selected_metadata) * args.total_restarts,
                    "cpu_lm_attempt_cap": len(selected_metadata) * args.lm_top_k,
                    "cpu_workers": args.cpu_workers,
                    "full_manifest_audit": full_audit,
                    "selected_candidate_indices": [item["candidate_index"] for item in selected],
                    "positive_control_variants": 3,
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
