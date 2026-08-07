#!/usr/bin/env python3
"""Batched Apple-MPS triage for the unresolved K7 affine-H systems.

This file is intentionally heuristic.  It first reproduces the exact prior
cover filters and exact strict-affine-H certificates on CPU.  Only the still
unresolved saturating-clique systems are sent to one padded MPS batch with
multiple simultaneous restarts.  The GPU minimizes the special-inverse
identities

    (1-S) H_ij + s_i s_j = 0                         (i != j)
    (1-S) H_ii - s_i(1-s_i) = 0,
    s=H1, S=1^T H 1,

together with the necessary strict affine inequalities.  Small residual is
not a realization and large residual is not a rejection.  The output is only
a ranking of targets for later exact algebra or interval work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as prior
import d6_k7_special_h_reference as href


Q = Fraction
MAX_U = 7


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _strict_inequalities(
    u: int,
    remainder: Sequence[int],
    columns: dict[int, tuple[int, ...]],
) -> list[tuple[tuple[Q, ...], Q, dict]]:
    output: list[tuple[tuple[Q, ...], Q, dict]] = []
    for i in range(u):
        output.append((
            href.matrix_entry_coefficients(u, i, i), Q(0),
            {"kind": "H_diagonal_positive", "basis_indices": [i]},
        ))
    for i, j in combinations(range(u), 2):
        output.append((
            tuple(-value for value in href.matrix_entry_coefficients(u, i, j)),
            Q(0),
            {"kind": "H_offdiagonal_negative", "basis_indices": [i, j]},
        ))
    for i in range(u):
        output.append((
            href.row_sum_coefficients(u, i), Q(0),
            {"kind": "H_row_sum_positive", "basis_indices": [i]},
        ))
    total = href.total_sum_coefficients(u)
    output.append((total, Q(0), {"kind": "H_total_sum_positive"}))
    output.append((
        tuple(-value for value in total), Q(-1),
        {"kind": "H_total_sum_below_one"},
    ))
    for vertex in remainder:
        output.append((
            href.quadratic_form_coefficients(columns[vertex]), Q(1),
            {"kind": "outside_diagonal_above_one", "vertices": [vertex]},
        ))
    for first, second in combinations(remainder, 2):
        difference = tuple(
            a - b for a, b in zip(columns[first], columns[second])
        )
        if any(difference):
            output.append((
                href.quadratic_form_coefficients(difference), Q(0),
                {
                    "kind": "column_difference_form_positive",
                    "vertices": [first, second],
                },
            ))
    return output


def build_problem(
    graph_n: Sequence[int],
    clique_mask: int,
    metadata: dict,
) -> dict:
    clique = list(prior.bits(clique_mask))
    u = len(clique)
    remainder = [
        vertex for vertex in range(len(graph_n))
        if not (clique_mask & (1 << vertex))
    ]
    columns = {
        vertex: tuple(
            int(bool(graph_n[vertex] & (1 << basis))) for basis in clique
        )
        for vertex in remainder
    }
    equation_rows = []
    equation_rhs = []
    for first, second in combinations(remainder, 2):
        equation_rows.append(
            href.bilinear_form_coefficients(columns[first], columns[second])
        )
        equation_rhs.append(int(bool(graph_n[first] & (1 << second))))
    system = href.exact_rref(
        equation_rows, equation_rhs, u * (u + 1) // 2
    )
    if not system.consistent:
        raise AssertionError("unresolved H problem unexpectedly inconsistent")
    origin, basis = system.parameterization()
    inequalities = _strict_inequalities(u, remainder, columns)
    affine_a = []
    affine_c = []
    for functional, lower_bound, _ in inequalities:
        constant = sum(a * b for a, b in zip(functional, origin))
        affine_a.append([
            float(sum(a * b for a, b in zip(functional, direction)))
            for direction in basis
        ])
        affine_c.append(float(lower_bound - constant))

    origin_matrix = [0.0] * (MAX_U * MAX_U)
    basis_matrix = [
        [0.0] * len(basis) for _ in range(MAX_U * MAX_U)
    ]
    for coordinate, (i, j) in enumerate(href.symmetric_coordinates(u)):
        origin_matrix[i * MAX_U + j] = float(origin[coordinate])
        origin_matrix[j * MAX_U + i] = float(origin[coordinate])
        for parameter, direction in enumerate(basis):
            basis_matrix[i * MAX_U + j][parameter] = float(
                direction[coordinate]
            )
            basis_matrix[j * MAX_U + i][parameter] = float(
                direction[coordinate]
            )
    return {
        **metadata,
        "u": u,
        "clique": clique,
        "remainder": remainder,
        "equation_rank": system.rank,
        "h_variables": system.variables,
        "free_parameters": len(basis),
        "inequalities": len(inequalities),
        "affine_A": affine_a,
        "affine_c": affine_c,
        "origin_matrix": origin_matrix,
        "basis_matrix": basis_matrix,
    }


def graph_problems(graph: dict) -> list[dict]:
    adj = tuple(graph["adjacency"])
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    output = []
    for seed_mask in prior.clique_masks(adj, 7):
        seed, outside, defects, ladj, eligible = prior.seed_instance(adj, seed_mask)
        total_term_rank = prior.matching_size(defects)
        for zmask in prior.eligible_covers(ladj, eligible):
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            nvertices = [
                outside[i] for i in range(len(outside))
                if not (zmask & (1 << i))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            for clique_mask in prior.clique_masks(
                graph_n, analysis.k_rank_upper
            ):
                exact = href.assess_saturating_clique(graph_n, clique_mask)
                if exact.failed:
                    continue
                output.append(build_problem(
                    graph_n,
                    clique_mask,
                    {
                        "graph_index": graph.get("index"),
                        "seed": seed,
                        "zmask": zmask,
                        "nvertices": nvertices,
                    },
                ))
    return output


def collect_problems(graphs: Sequence[dict], workers: int) -> list[dict]:
    if workers == 1:
        nested = [graph_problems(graph) for graph in graphs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            nested = list(executor.map(graph_problems, graphs, chunksize=1))
    return [problem for group in nested for problem in group]


def mps_optimize(
    problems: Sequence[dict],
    *,
    restarts: int,
    steps: int,
    learning_rate: float,
    seed: int,
) -> dict:
    import torch

    if not torch.backends.mps.is_available():
        raise RuntimeError("Apple MPS is unavailable; run outside the sandbox")
    torch.manual_seed(seed)
    count = len(problems)
    if not count:
        return {"systems": 0, "targets": []}
    max_parameters = max(problem["free_parameters"] for problem in problems)
    max_inequalities = max(problem["inequalities"] for problem in problems)
    device = torch.device("mps")
    dtype = torch.float32

    affine_a = torch.zeros(
        (count, max_inequalities, max_parameters), dtype=dtype
    )
    affine_c = torch.zeros((count, max_inequalities), dtype=dtype)
    inequality_mask = torch.zeros((count, max_inequalities), dtype=dtype)
    origin = torch.zeros((count, MAX_U * MAX_U), dtype=dtype)
    basis = torch.zeros(
        (count, MAX_U * MAX_U, max_parameters), dtype=dtype
    )
    matrix_mask = torch.zeros((count, MAX_U, MAX_U), dtype=dtype)
    parameter_mask = torch.zeros((count, max_parameters), dtype=dtype)
    for index, problem in enumerate(problems):
        d = problem["free_parameters"]
        m = problem["inequalities"]
        u = problem["u"]
        if d:
            affine_a[index, :m, :d] = torch.tensor(
                problem["affine_A"], dtype=dtype
            )
            for flat, row in enumerate(problem["basis_matrix"]):
                basis[index, flat, :d] = torch.tensor(row, dtype=dtype)
            parameter_mask[index, :d] = 1
        affine_c[index, :m] = torch.tensor(problem["affine_c"], dtype=dtype)
        inequality_mask[index, :m] = 1
        origin[index] = torch.tensor(problem["origin_matrix"], dtype=dtype)
        matrix_mask[index, :u, :u] = 1

    # Positive row scaling preserves every strict inequality and improves the
    # common batched optimizer's conditioning.
    scale = torch.maximum(
        torch.linalg.vector_norm(affine_a, dim=2), affine_c.abs()
    ).clamp_min(1.0)
    affine_a = affine_a / scale[:, :, None]
    affine_c = affine_c / scale
    affine_a = affine_a.to(device)
    affine_c = affine_c.to(device)
    inequality_mask = inequality_mask.to(device)
    origin = origin.to(device)
    basis = basis.to(device)
    matrix_mask = matrix_mask.to(device)
    parameter_mask = parameter_mask.to(device)

    parameters = torch.nn.Parameter(
        0.25 * torch.randn(
            (count, restarts, max_parameters), device=device, dtype=dtype
        )
    )
    with torch.no_grad():
        parameters[:, 0, :] = 0
        parameters.mul_(parameter_mask[:, None, :])
    optimizer = torch.optim.Adam([parameters], lr=learning_rate)
    identity_mask = matrix_mask.clone()
    diagonal_mask = torch.eye(MAX_U, device=device)[None, :, :] * matrix_mask
    offdiagonal_mask = identity_mask - diagonal_mask
    inequality_count = inequality_mask.sum(dim=1).clamp_min(1)[:, None]
    identity_count = identity_mask.sum(dim=(1, 2)).clamp_min(1)[:, None]

    started = time.perf_counter()
    for _ in range(steps):
        optimizer.zero_grad(set_to_none=True)
        parameters.data.mul_(parameter_mask[:, None, :])
        slack = (
            torch.einsum("bmd,brd->brm", affine_a, parameters)
            - affine_c[:, None, :]
        )
        hinge = torch.relu(1e-3 - slack) * inequality_mask[:, None, :]
        inequality_loss = (hinge.square().sum(dim=2) / inequality_count)

        hflat = origin[:, None, :] + torch.einsum(
            "bfd,brd->brf", basis, parameters
        )
        hmatrix = hflat.reshape(count, restarts, MAX_U, MAX_U)
        row_sum = hmatrix.sum(dim=3)
        total_sum = row_sum.sum(dim=2)
        product = row_sum[:, :, :, None] * row_sum[:, :, None, :]
        residual = (
            (1 - total_sum)[:, :, None, None] * hmatrix + product
        )
        diagonal_residual = (
            (1 - total_sum)[:, :, None]
            * torch.diagonal(hmatrix, dim1=2, dim2=3)
            - row_sum * (1 - row_sum)
        )
        residual = residual * offdiagonal_mask[:, None, :, :]
        residual = residual + torch.diag_embed(diagonal_residual)
        identity_loss = (
            residual.square().sum(dim=(2, 3)) / identity_count
        )
        regularizer = 1e-10 * parameters.square().mean(dim=2)
        loss = (identity_loss + 20 * inequality_loss + regularizer).mean()
        loss.backward()
        optimizer.step()
    torch.mps.synchronize()
    elapsed = time.perf_counter() - started

    with torch.no_grad():
        slack = (
            torch.einsum("bmd,brd->brm", affine_a, parameters)
            - affine_c[:, None, :]
        )
        masked_slack = torch.where(
            inequality_mask[:, None, :].bool(),
            slack,
            torch.full_like(slack, float("inf")),
        )
        minimum_slack = masked_slack.min(dim=2).values
        hflat = origin[:, None, :] + torch.einsum(
            "bfd,brd->brf", basis, parameters
        )
        hmatrix = hflat.reshape(count, restarts, MAX_U, MAX_U)
        row_sum = hmatrix.sum(dim=3)
        total_sum = row_sum.sum(dim=2)
        product = row_sum[:, :, :, None] * row_sum[:, :, None, :]
        residual = (
            (1 - total_sum)[:, :, None, None] * hmatrix + product
        ) * offdiagonal_mask[:, None, :, :]
        diagonal_residual = (
            (1 - total_sum)[:, :, None]
            * torch.diagonal(hmatrix, dim1=2, dim2=3)
            - row_sum * (1 - row_sum)
        )
        residual = residual + torch.diag_embed(diagonal_residual)
        rms = torch.sqrt(
            residual.square().sum(dim=(2, 3)) / identity_count
        )
        score = rms + 10 * torch.relu(-minimum_slack)
        best_restart = score.argmin(dim=1)
        rows = torch.arange(count, device=device)
        best_rms = rms[rows, best_restart].cpu().tolist()
        best_slack = minimum_slack[rows, best_restart].cpu().tolist()

    targets = []
    for index, problem in enumerate(problems):
        targets.append({
            "problem": index,
            "graph_index": problem["graph_index"],
            "seed": problem["seed"],
            "zmask": problem["zmask"],
            "nvertices": problem["nvertices"],
            "u": problem["u"],
            "clique": problem["clique"],
            "free_parameters": problem["free_parameters"],
            "best_special_identity_rms": best_rms[index],
            "best_normalized_minimum_strict_slack": best_slack[index],
        })
    targets.sort(
        key=lambda item: (
            item["best_special_identity_rms"]
            + 10 * max(0.0, -item["best_normalized_minimum_strict_slack"])
        ),
        reverse=True,
    )
    return {
        "systems": count,
        "device": str(device),
        "dtype": str(dtype),
        "batch_shape": [count, restarts],
        "max_free_parameters": max_parameters,
        "max_inequalities": max_inequalities,
        "steps": steps,
        "learning_rate": learning_rate,
        "seed": seed,
        "wall_seconds": elapsed,
        "targets": targets,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("selection_report", type=Path)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--restarts", type=int, default=32)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--learning-rate", type=float, default=0.03)
    parser.add_argument("--seed", type=int, default=600719)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    sample = json.loads(args.sample.read_text(encoding="utf-8"))
    selection = json.loads(args.selection_report.read_text(encoding="utf-8"))
    indices = set(selection["graph_decisions"]["combined_survivors"])
    graphs = [
        graph for graph in sample["graphs"] if graph.get("index") in indices
    ]
    extraction_started = time.perf_counter()
    problems = collect_problems(graphs, args.workers)
    extraction_seconds = time.perf_counter() - extraction_started
    triage = mps_optimize(
        problems,
        restarts=args.restarts,
        steps=args.steps,
        learning_rate=args.learning_rate,
        seed=args.seed,
    )
    report = {
        "schema": 1,
        "description": (
            "Heuristic batched Apple-MPS optimization of exact-strict-H "
            "survivors. No numerical result is a proof or rejection."
        ),
        "proof_status": "HEURISTIC_ONLY",
        "selected_graphs": len(graphs),
        "selected_indices": sorted(indices),
        "cpu_extraction_seconds": extraction_seconds,
        "cpu_workers": args.workers,
        "input_sample": {
            "file": args.sample.name,
            "sha256": file_sha256(args.sample),
        },
        "selection_report": {
            "file": args.selection_report.name,
            "sha256": file_sha256(args.selection_report),
        },
        "source": {
            "file": Path(__file__).name,
            "sha256": file_sha256(Path(__file__)),
        },
        "triage": triage,
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
