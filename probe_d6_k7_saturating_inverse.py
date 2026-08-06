#!/usr/bin/env python3
"""Bounded exact profile of special inverses on the final K7 cover residue.

This is a profiling/discovery program, not a production certificate checker.
It reconstructs the 189 survivors of the independently checked full joint
cover/support campaign from its committed inputs (never from ``.runs``
checkpoints), identifies every jointly surviving saturating-clique cover, and
reduces its affine system for ``H=K[C,C]^-1`` over ``fractions.Fraction``.

For a unique affine H it independently checks the affine equations, all
strict inequalities used by the strict-H layer, and the special-inverse
identities

    r = H 1,  delta = 1 - 1^T r,
    delta H_ij + r_i r_j = 0             (i != j),
    delta H_ii - r_i (1-r_i) = 0.

The current residue turns out to have no unique or one-parameter affine
systems.  As a deliberately bounded follow-up, the program sends only the
first three systems of each rank ``u`` through modular graded-order Groebner
screens in the direct Sherman--Morrison variables.  Each child has a hard
wall-clock timeout.  Modular output is discovery only; no modular failure or
timeout is treated as a rejection.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation
import run_d6_k7_support_capacity_pilot as joint_reference


Q = Fraction
ROOT = Path(__file__).resolve().parent
RANK_INPUT = ROOT / ".runs/d6_k7_rank_survivors.json"
JOINT_REPORT = ROOT / "d6_k7_joint_support_full_report.json"
JOINT_DECISIONS = ROOT / "d6_k7_joint_support_full_decisions.tsv.gz"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"

EXPECTED_JOINT_REPORT_SHA256 = (
    "b3fd2752a371f8909b79b50b8e2b6135f93423d0adf7286aee26d4a8c0643326"
)
EXPECTED_JOINT_DECISIONS_SHA256 = (
    "b6040e7796e3d7df511a4e460ad71082e63c074ac2cf2713c19d8f3b772bcc90"
)
EXPECTED_SURVIVORS = 189


def sha256(path: Path) -> str:
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


def fraction_json(value: Q) -> list[int]:
    value = Q(value)
    return [value.numerator, value.denominator]


def bits(mask: int) -> tuple[int, ...]:
    output = []
    while mask:
        bit = mask & -mask
        output.append(bit.bit_length() - 1)
        mask ^= bit
    return tuple(output)


def symmetric_coordinates(u: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(u) for j in range(i, u))


def bilinear_coefficients(left: Sequence[int], right: Sequence[int]) -> list[Q]:
    output = []
    for i, j in symmetric_coordinates(len(left)):
        if i == j:
            output.append(Q(left[i] * right[i]))
        else:
            output.append(Q(left[i] * right[j] + left[j] * right[i]))
    return output


def rref(
    coefficients: Sequence[Sequence[int | Q]],
    rhs: Sequence[int | Q],
    variables: int,
) -> dict:
    """Independent exact RREF; returns a rational affine parameterization."""

    if len(coefficients) != len(rhs):
        raise ValueError("coefficient/rhs count mismatch")
    rows = [
        [Q(value) for value in row] + [Q(target)]
        for row, target in zip(coefficients, rhs)
    ]
    if any(len(row) != variables + 1 for row in rows):
        raise ValueError("bad exact linear-system width")
    pivot_row = 0
    pivots: list[int] = []
    for column in range(variables):
        selected = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column]),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    value - factor * source
                    for value, source in zip(rows[row], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    for row in rows[pivot_row:]:
        if not any(row[:-1]) and row[-1]:
            return {
                "consistent": False,
                "rank": len(pivots),
                "pivots": tuple(pivots),
                "origin": None,
                "directions": None,
            }
    free = [column for column in range(variables) if column not in set(pivots)]
    origin = [Q(0)] * variables
    directions = [[Q(0)] * variables for _ in free]
    for parameter, column in enumerate(free):
        directions[parameter][column] = Q(1)
    for row, pivot in zip(rows[:pivot_row], pivots):
        origin[pivot] = row[-1]
        for parameter, column in enumerate(free):
            directions[parameter][pivot] = -row[column]
    return {
        "consistent": True,
        "rank": len(pivots),
        "pivots": tuple(pivots),
        "origin": tuple(origin),
        "directions": tuple(tuple(direction) for direction in directions),
    }


def rank_of_rows(rows: Sequence[Sequence[int | Q]], width: int) -> int:
    return int(rref(rows, [0] * len(rows), width)["rank"])


def symmetric_matrix(solution: Sequence[Q], u: int) -> list[list[Q]]:
    matrix = [[Q(0) for _ in range(u)] for _ in range(u)]
    for value, (i, j) in zip(solution, symmetric_coordinates(u)):
        matrix[i][j] = matrix[j][i] = Q(value)
    return matrix


def determinant(matrix: Sequence[Sequence[Q]]) -> Q:
    work = [[Q(value) for value in row] for row in matrix]
    n = len(work)
    if any(len(row) != n for row in work):
        raise ValueError("determinant needs a square matrix")
    result = Q(1)
    for column in range(n):
        selected = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if selected is None:
            return Q(0)
        if selected != column:
            work[column], work[selected] = work[selected], work[column]
            result = -result
        pivot = work[column][column]
        result *= pivot
        for row in range(column + 1, n):
            factor = work[row][column] / pivot
            for entry in range(column + 1, n):
                work[row][entry] -= factor * work[column][entry]
    return result


def inverse(matrix: Sequence[Sequence[Q]]) -> list[list[Q]] | None:
    n = len(matrix)
    work = [
        [Q(value) for value in row]
        + [Q(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    if any(len(row) != 2 * n for row in work):
        raise ValueError("inverse needs a square matrix")
    for column in range(n):
        selected = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if selected is None:
            return None
        work[column], work[selected] = work[selected], work[column]
        pivot = work[column][column]
        work[column] = [value / pivot for value in work[column]]
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    value - factor * source
                    for value, source in zip(work[row], work[column])
                ]
    return [row[n:] for row in work]


def quadratic(matrix: Sequence[Sequence[Q]], vector: Sequence[int]) -> Q:
    return sum(
        Q(vector[i]) * matrix[i][j] * Q(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def special_inverse_check(matrix: Sequence[Sequence[Q]]) -> dict:
    """Check the special-inverse characterization and its strict domain."""

    u = len(matrix)
    if any(len(row) != u for row in matrix):
        raise ValueError("H is not square")
    if any(matrix[i][j] != matrix[j][i] for i in range(u) for j in range(u)):
        return {"valid": False, "failures": ["not_symmetric"]}
    r = [sum(matrix[i]) for i in range(u)]
    delta = Q(1) - sum(r)
    failures = []
    if delta <= 0:
        failures.append("delta_not_positive")
    if any(value <= 0 for value in r):
        failures.append("row_sum_not_positive")
    for i in range(u):
        if delta * matrix[i][i] != r[i] * (1 - r[i]):
            failures.append(f"diagonal_identity_{i}")
    for i, j in combinations(range(u), 2):
        if delta * matrix[i][j] != -r[i] * r[j]:
            failures.append(f"offdiagonal_identity_{i}_{j}")
    recovered = inverse(matrix)
    if recovered is None:
        failures.append("singular_H")
    else:
        for i, j in combinations(range(u), 2):
            if recovered[i][j] != 1:
                failures.append(f"inverse_offdiagonal_{i}_{j}")
        for i in range(u):
            if recovered[i][i] <= 1:
                failures.append(f"inverse_diagonal_{i}")
    return {
        "valid": not failures,
        "failures": failures,
        "r": [fraction_json(value) for value in r],
        "delta": fraction_json(delta),
    }


def check_unique_solution(
    solution: Sequence[Q],
    columns: Sequence[Sequence[int]],
    equation_pairs: Sequence[tuple[int, int]],
    targets: Sequence[int],
    coefficients: Sequence[Sequence[Q]],
    u: int,
) -> dict:
    """Independently check all equations and strict conditions for unique H."""

    failures = []
    for row, target in zip(coefficients, targets):
        if sum(a * b for a, b in zip(row, solution)) != target:
            failures.append("affine_equation")
    matrix = symmetric_matrix(solution, u)
    for i in range(u):
        if matrix[i][i] <= 0:
            failures.append(f"H_diagonal_{i}")
    for i, j in combinations(range(u), 2):
        if matrix[i][j] >= 0:
            failures.append(f"H_offdiagonal_{i}_{j}")
    row_sums = [sum(row) for row in matrix]
    if any(value <= 0 for value in row_sums):
        failures.append("H_row_sum")
    total = sum(row_sums)
    if not (0 < total < 1):
        failures.append("H_total_sum")
    for vertex, column in enumerate(columns):
        if quadratic(matrix, column) <= 1:
            failures.append(f"outside_diagonal_{vertex}")
    for first, second in combinations(range(len(columns)), 2):
        difference = [
            columns[first][i] - columns[second][i] for i in range(u)
        ]
        if not any(difference) or quadratic(matrix, difference) <= 0:
            failures.append(f"column_difference_{first}_{second}")
    if not all(
        determinant([row[:size] for row in matrix[:size]]) > 0
        for size in range(1, u + 1)
    ):
        failures.append("not_positive_definite")
    special = special_inverse_check(matrix)
    failures.extend(f"special:{name}" for name in special["failures"])
    return {
        "valid": not failures,
        "failures": sorted(set(failures)),
        "special": special,
        "equations_checked": len(equation_pairs),
    }


def sherman_polynomial_rows(masks: Sequence[int], targets: Sequence[int], u: int) -> dict:
    """Linearize the ten degree-at-most-two equations in all monomials."""

    monomials = [(0,)]
    monomials.extend((1, i) for i in range(u))
    monomials.extend((2, i, j) for i in range(u) for j in range(i, u))
    position = {monomial: index for index, monomial in enumerate(monomials)}
    rows = []
    pair_targets = []
    for pair_index, (first, second) in enumerate(combinations(range(len(masks)), 2)):
        target = int(targets[pair_index])
        row = [0] * len(monomials)
        left = tuple(bits(masks[first]))
        right = tuple(bits(masks[second]))
        intersection = tuple(sorted(set(left) & set(right)))
        row[position[(0,)]] -= target
        for i in range(u):
            row[position[(1, i)]] -= target
        for i in intersection:
            row[position[(1, i)]] += 1
            for j in range(u):
                row[position[(2, min(i, j), max(i, j))]] += 1
        for i in left:
            for j in right:
                row[position[(2, min(i, j), max(i, j))]] -= 1
        rows.append(row)
        pair_targets.append([first, second, target])
    rank = rank_of_rows(rows, len(monomials))
    return {
        "monomial_count": len(monomials),
        "rank": rank,
        "free_dimension": len(monomials) - rank,
        "pair_targets": pair_targets,
    }


def affine_system_for_clique(adj: Sequence[int], clique_mask: int) -> dict:
    basis = list(bits(clique_mask))
    u = len(basis)
    remainder = [vertex for vertex in range(len(adj)) if vertex not in basis]
    columns = [
        tuple(int(bool(adj[vertex] & (1 << member))) for member in basis)
        for vertex in remainder
    ]
    coefficients = []
    targets = []
    equation_pairs = []
    for first, second in combinations(range(len(remainder)), 2):
        coefficients.append(bilinear_coefficients(columns[first], columns[second]))
        targets.append(int(bool(adj[remainder[first]] & (1 << remainder[second]))))
        equation_pairs.append((first, second))
    variables = u * (u + 1) // 2
    reduced = rref(coefficients, targets, variables)
    if not reduced["consistent"]:
        # Such a clique should already have failed the committed strict-H layer.
        unique = {"valid": False, "failures": ["affine_inconsistent"]}
    elif reduced["rank"] == variables:
        unique = check_unique_solution(
            reduced["origin"], columns, equation_pairs, targets, coefficients, u
        )
    else:
        unique = None
    masks = [sum(value << i for i, value in enumerate(column)) for column in columns]
    polynomial = sherman_polynomial_rows(masks, targets, u)
    return {
        "u": u,
        "basis": basis,
        "remainder": remainder,
        "masks": masks,
        "targets": targets,
        "equations": len(coefficients),
        "variables": variables,
        "rank": int(reduced["rank"]),
        "free_dimension": variables - int(reduced["rank"]),
        "consistent": bool(reduced["consistent"]),
        "unique_check": unique,
        "linearized_polynomial": polynomial,
    }


def joint_cover_survives(
    graph_n: Sequence[int],
    z_allowed: Sequence[int],
    n_allowed: Sequence[int],
    zero_forcing: reference.ZeroForcingSolver,
    clique_solver: reference.CliqueStructureSolver,
) -> bool:
    """Reconstruct the exact pre-capacity existential support conjunction."""

    for z_supports in propagation.labeled_support_families(z_allowed):
        propagated = propagation.analyze_support_assignment(
            graph_n, z_supports, n_allowed, zero_forcing, clique_solver
        )
        if propagated.failure is not None:
            continue
        if sparse_value.check_small_support_masks(
            graph_n, propagated.propagated_masks
        ).feasible:
            return True
    return False


def analyze_graph(payload: dict) -> dict:
    graph = payload["graph"]
    expected = payload["expected"]
    prior_failures = {
        (tuple(seed), int(zmask)) for seed, zmask in payload["prior_failures"]
    }
    tetrad_failures = {
        (tuple(seed), int(zmask)) for seed, zmask in payload["tetrad_failures"]
    }
    index = int(graph["index"])
    adj = tuple(int(row) for row in graph["adjacency"])
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    counts: Counter[str] = Counter()
    systems = []
    direct_rejected_seeds = 0
    seed_profiles = []
    for seed_mask in reference.clique_masks(adj, 7):
        counts["seeds"] += 1
        seed, outside, defects, ladj, eligible = reference.seed_instance(adj, seed_mask)
        total_term_rank = reference.matching_size(defects)
        jointly_surviving = 0
        direct_surviving = 0
        for zmask in reference.eligible_covers(ladj, eligible):
            counts["covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            key = (tuple(seed), zmask)
            status = joint_reference.current_cover_status(
                adj, outside, defects, zmask, baseline, key,
                prior_failures, tetrad_failures,
            )
            counts[f"current_{status}"] += 1
            if status != "passing":
                continue
            nlocal = tuple(
                vertex for vertex in range(len(outside))
                if not zmask & (1 << vertex)
            )
            zlocal = tuple(bits(zmask))
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nlocal]
            )
            if not joint_cover_survives(
                graph_n,
                [defects[vertex] for vertex in zlocal],
                [defects[vertex] for vertex in nlocal],
                zero_forcing,
                clique_solver,
            ):
                counts["joint_failing"] += 1
                continue
            counts["joint_passing"] += 1
            jointly_surviving += 1
            cover_failed = False
            cliques = tuple(reference.clique_masks(graph_n, baseline.k_rank_upper)) \
                if baseline.k_rank_upper else ()
            counts[f"saturating_cliques_per_cover:{len(cliques)}"] += 1
            if not cliques:
                counts["covers_without_saturating_clique"] += 1
            for clique_mask in cliques:
                system = affine_system_for_clique(graph_n, clique_mask)
                system.update({
                    "graph_index": index,
                    "seed_mask": int(seed_mask),
                    "seed": list(seed),
                    "zmask": int(zmask),
                    "nvertices": [outside[vertex] for vertex in nlocal],
                    "clique_mask": int(clique_mask),
                })
                systems.append(system)
                counts["saturating_cliques"] += 1
                if system["unique_check"] is not None and not system["unique_check"]["valid"]:
                    cover_failed = True
            if not cover_failed:
                direct_surviving += 1
            else:
                counts["direct_unique_rejected_covers"] += 1
        if not jointly_surviving:
            raise AssertionError(
                f"joint survivor graph {index} has a seed without a joint cover"
            )
        if not direct_surviving:
            direct_rejected_seeds += 1
        seed_profiles.append({
            "seed_mask": int(seed_mask),
            "jointly_surviving_covers": jointly_surviving,
            "directly_surviving_covers": direct_surviving,
        })
    if counts["current_passing"] != int(expected["current_passing_covers"]):
        raise ValueError(
            f"graph {index}: current cover mismatch "
            f"{counts['current_passing']} != {expected['current_passing_covers']}"
        )
    if counts["joint_passing"] != int(expected["pre_capacity_passing_covers"]):
        raise ValueError(
            f"graph {index}: joint cover mismatch "
            f"{counts['joint_passing']} != {expected['pre_capacity_passing_covers']}"
        )
    return {
        "index": index,
        "counts": dict(counts),
        "systems": systems,
        "seed_profiles": seed_profiles,
        "direct_unique_rejected": bool(direct_rejected_seeds),
    }


def synthetic_controls() -> dict:
    """Feasible and perturbed exact special-inverse controls."""

    e = [Q(2), Q(3), Q(5)]
    k = [[Q(1) + (e[i] if i == j else 0) for j in range(3)] for i in range(3)]
    h = inverse(k)
    assert h is not None
    positive = special_inverse_check(h)
    if not positive["valid"]:
        raise AssertionError(f"feasible special inverse failed: {positive}")
    coordinates = [h[i][j] for i, j in symmetric_coordinates(3)]
    identity = [[Q(int(i == j)) for j in range(6)] for i in range(6)]
    reduced = rref(identity, coordinates, 6)
    if reduced["origin"] != tuple(coordinates) or reduced["rank"] != 6:
        raise AssertionError("synthetic unique affine system failed")
    perturbed = [row[:] for row in h]
    perturbed[0][1] += Q(1, 100)
    perturbed[1][0] += Q(1, 100)
    negative = special_inverse_check(perturbed)
    if negative["valid"]:
        raise AssertionError("perturbed nonspecial H passed")
    return {
        "status": "PASS",
        "feasible_e": [fraction_json(value) for value in e],
        "feasible_delta": positive["delta"],
        "feasible_r": positive["r"],
        "perturbed_failure_count": len(negative["failures"]),
    }


def groebner_worker(encoded: str) -> None:
    """Child-process entry point; prints exactly one compact JSON object."""

    payload = json.loads(base64.urlsafe_b64decode(encoded.encode("ascii")))
    started_unix = time.time()
    started = time.perf_counter()
    try:
        import sympy as sp

        u = int(payload["u"])
        xs = sp.symbols(f"x0:{u}")
        total = 1 + sum(xs)
        masks = [int(value) for value in payload["masks"]]
        polynomials = []
        for target, (first, second) in zip(
            payload["targets"], combinations(range(len(masks)), 2)
        ):
            left = sum(xs[i] for i in bits(masks[first]))
            right = sum(xs[i] for i in bits(masks[second]))
            intersection = sum(
                xs[i] for i in bits(masks[first] & masks[second])
            )
            polynomials.append(
                sp.expand(total * intersection - left * right - int(target) * total)
            )
        modulus = payload.get("modulus")
        kwargs = {"order": "grlex"}
        if modulus is not None:
            kwargs["modulus"] = int(modulus)
        basis = sp.groebner(polynomials, *xs, **kwargs)
        polys = list(basis.polys)
        unit = len(polys) == 1 and polys[0].total_degree() == 0
        output = {
            "status": "COMPLETE",
            "pid": os.getpid(),
            "started_unix": started_unix,
            "finished_unix": time.time(),
            "modulus": modulus,
            "unit_ideal": bool(unit),
            "basis_polynomials": len(polys),
            "maximum_total_degree": max((poly.total_degree() for poly in polys), default=0),
            "maximum_terms": max((len(poly.terms()) for poly in polys), default=0),
            "total_terms": sum(len(poly.terms()) for poly in polys),
            "elapsed_seconds": time.perf_counter() - started,
        }
    except BaseException as error:  # noqa: BLE001 - isolated discovery child
        output = {
            "status": "ERROR",
            "pid": os.getpid(),
            "started_unix": started_unix,
            "finished_unix": time.time(),
            "error_type": type(error).__name__,
            "error": str(error),
            "elapsed_seconds": time.perf_counter() - started,
        }
    print(json.dumps(output, sort_keys=True, separators=(",", ":")))


def run_groebner_job(system: dict, modulus: int | None, timeout: float) -> dict:
    payload = {
        "u": system["u"],
        "masks": system["masks"],
        "targets": system["targets"],
        "modulus": modulus,
    }
    encoded = base64.urlsafe_b64encode(
        json.dumps(payload, separators=(",", ":")).encode("ascii")
    ).decode("ascii")
    started_unix = time.time()
    started = time.perf_counter()
    environment = dict(os.environ)
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = "1"
    process = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "--groebner-worker", encoded],
            cwd=ROOT,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    pid = process.pid
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            stdout, stderr = process.communicate(timeout=1.0)
            stop = "SIGTERM"
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            stop = "SIGKILL"
        return {
            "status": "TIMEOUT",
            "modulus": modulus,
            "pid": pid,
            "started_unix": started_unix,
            "finished_unix": time.time(),
            "stop": stop,
            "timeout_seconds": timeout,
            "elapsed_seconds": time.perf_counter() - started,
        }
    lines = [line for line in stdout.splitlines() if line.strip()]
    if process.returncode or len(lines) != 1:
        return {
            "status": "INFRA_ERROR",
            "modulus": modulus,
            "pid": pid,
            "started_unix": started_unix,
            "finished_unix": time.time(),
            "stop": "EXIT",
            "returncode": process.returncode,
            "stderr": stderr[-2000:],
            "stdout": stdout[-2000:],
            "elapsed_seconds": time.perf_counter() - started,
        }
    output = json.loads(lines[0])
    if int(output.get("pid", -1)) != pid:
        raise ValueError("Groebner child PID mismatch")
    output["stop"] = "EXIT"
    output["wall_seconds"] = time.perf_counter() - started
    return output


def bounded_groebner_pilot(
    systems: Sequence[dict],
    primes: Sequence[int],
    timeout: float,
    workers: int,
) -> list[dict]:
    selected = []
    for u in sorted({int(system["u"]) for system in systems}):
        selected.extend([system for system in systems if int(system["u"]) == u][:3])
    records = [{
        "identity": {
            key: system[key]
            for key in ("graph_index", "seed_mask", "zmask", "clique_mask", "u")
        },
        "masks": system["masks"],
        "targets": system["targets"],
        "linearized_polynomial": system["linearized_polynomial"],
        "modular": {},
        "rational": None,
    } for system in selected]
    jobs = {}
    with ThreadPoolExecutor(max_workers=max(1, min(workers, 6))) as executor:
        for position, system in enumerate(selected):
            for prime in primes:
                future = executor.submit(run_groebner_job, system, int(prime), timeout)
                jobs[future] = (position, int(prime))
        for future in as_completed(jobs):
            position, prime = jobs[future]
            records[position]["modular"][str(prime)] = future.result()
    rational_positions = [
        position for position, record in enumerate(records)
        if any(result.get("unit_ideal") for result in record["modular"].values())
    ]
    if rational_positions:
        with ThreadPoolExecutor(max_workers=max(1, min(workers, 3))) as executor:
            jobs = {
                executor.submit(run_groebner_job, selected[position], None, timeout): position
                for position in rational_positions
            }
            for future in as_completed(jobs):
                records[jobs[future]]["rational"] = future.result()
    return records


def load_boundary() -> tuple[list[dict], dict[int, dict], dict[str, str]]:
    inherited = joint_reference.verify_inputs()
    observed_report = sha256(JOINT_REPORT)
    observed_decisions = sha256(JOINT_DECISIONS)
    if observed_report != EXPECTED_JOINT_REPORT_SHA256:
        raise ValueError("joint report hash mismatch")
    if observed_decisions != EXPECTED_JOINT_DECISIONS_SHA256:
        raise ValueError("joint decisions hash mismatch")
    report = json.loads(JOINT_REPORT.read_text(encoding="utf-8"))
    if report.get("status") != "COMPLETE" or report["summary"]["status_counts"] != {
        "JOINT_PRE_CAPACITY_REJECTED": 69,
        "SURVIVOR": 189,
    }:
        raise ValueError("unexpected joint report boundary")
    rows = []
    with gzip.open(JOINT_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["status"] == "SURVIVOR":
                rows.append(row)
    rows.sort(key=lambda row: int(row["ordinal"]))
    if len(rows) != EXPECTED_SURVIVORS:
        raise ValueError("unexpected joint survivor count")
    payload = json.loads(RANK_INPUT.read_text(encoding="utf-8"))
    graphs = {int(graph["index"]): graph for graph in payload["graphs"]}
    selected = [int(row["index"]) for row in rows]
    if any(index not in graphs for index in selected):
        raise ValueError("joint survivor absent from pinned rank input")
    hashes = {
        **inherited,
        "d6_k7_joint_support_full_report.json": observed_report,
        "d6_k7_joint_support_full_decisions.tsv.gz": observed_decisions,
    }
    return rows, graphs, hashes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--groebner-timeout", type=float, default=8.0)
    parser.add_argument("--primes", default="101,103,107")
    parser.add_argument(
        "--output", type=Path,
        default=Path("/private/tmp/d6_k7_saturating_inverse_report.json"),
    )
    parser.add_argument("--groebner-worker", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.groebner_worker is not None:
        groebner_worker(args.groebner_worker)
        return
    if args.workers <= 0 or args.groebner_timeout <= 0:
        raise SystemExit("workers and timeout must be positive")
    primes = tuple(int(value) for value in args.primes.split(",") if value)
    if not (2 <= len(primes) <= 3):
        raise SystemExit("use two or three modular screening primes")

    controls = synthetic_controls()
    rows, graphs, input_hashes = load_boundary()
    selected = [int(row["index"]) for row in rows]
    selected_set = set(selected)
    prior = joint_reference.load_jsonl_witness_keys(
        PRIOR_CERTIFICATES, selected_set, "dual_failure_witnesses"
    )
    tetrad = joint_reference.load_jsonl_witness_keys(
        TETRAD_CERTIFICATES, selected_set, "tetrad_failure_witnesses"
    )
    joint_reference.load_sparse_survivors(selected_set)
    payloads = [{
        "graph": graphs[index],
        "expected": row,
        "prior_failures": [[list(seed), zmask] for seed, zmask in sorted(prior.get(index, set()))],
        "tetrad_failures": [[list(seed), zmask] for seed, zmask in sorted(tetrad.get(index, set()))],
    } for index, row in zip(selected, rows)]

    started = time.perf_counter()
    if args.workers == 1:
        graph_results = [analyze_graph(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            graph_results = list(executor.map(analyze_graph, payloads, chunksize=1))
    reconstruction_seconds = time.perf_counter() - started
    if [result["index"] for result in graph_results] != selected:
        raise ValueError("parallel reconstruction changed graph order")

    totals: Counter[str] = Counter()
    systems = []
    direct_rejected = []
    for result in graph_results:
        totals.update(result["counts"])
        systems.extend(result["systems"])
        if result["direct_unique_rejected"]:
            direct_rejected.append(result["index"])
    systems.sort(key=lambda system: (
        system["graph_index"], system["seed_mask"], system["zmask"], system["clique_mask"]
    ))
    affine_histogram = Counter(
        f"u{system['u']}_rank{system['rank']}_vars{system['variables']}_free{system['free_dimension']}"
        for system in systems
    )
    polynomial_histogram = Counter(
        f"u{system['u']}_rank{system['linearized_polynomial']['rank']}_"
        f"monomials{system['linearized_polynomial']['monomial_count']}_"
        f"free{system['linearized_polynomial']['free_dimension']}"
        for system in systems
    )
    unique_systems = sum(system["free_dimension"] == 0 for system in systems)
    dimension_one_systems = sum(system["free_dimension"] == 1 for system in systems)
    groebner = bounded_groebner_pilot(
        systems, primes, args.groebner_timeout, args.workers
    )
    modular_units = sum(
        any(result.get("unit_ideal") for result in record["modular"].values())
        for record in groebner
    )
    rational_units = sum(
        bool(record["rational"] and record["rational"].get("unit_ideal"))
        for record in groebner
    )

    report = {
        "schema": 1,
        "kind": "d6_k7_saturating_inverse_bounded_probe",
        "claim_scope": (
            "Exact reconstruction/profile and unique-H checks only. Modular Groebner "
            "screens and timeouts are discovery-only and reject nothing."
        ),
        "input": {
            "joint_survivors": len(selected),
            "ordered_indices_sha256": stable_hash(selected),
            "source_sha256": dict(sorted(input_hashes.items())),
        },
        "configuration": {
            "workers": args.workers,
            "groebner_primes": list(primes),
            "groebner_order": "grlex",
            "groebner_timeout_seconds_per_job": args.groebner_timeout,
            "groebner_systems_per_u": 3,
        },
        "controls": controls,
        "exact_profile": {
            "counts": dict(sorted(totals.items())),
            "affine_rank_histogram": dict(sorted(affine_histogram.items())),
            "linearized_degree_two_histogram": dict(sorted(polynomial_histogram.items())),
            "unique_affine_systems": unique_systems,
            "dimension_one_affine_systems": dimension_one_systems,
            "direct_unique_rejected_graphs": direct_rejected,
            "graphs_with_saturating_clique": sum(bool(result["systems"]) for result in graph_results),
            "graphs_without_saturating_clique": sum(not result["systems"] for result in graph_results),
        },
        "bounded_groebner_discovery": {
            "systems": groebner,
            "systems_screened": len(groebner),
            "systems_with_any_modular_unit_basis": modular_units,
            "systems_with_rational_unit_basis": rational_units,
            "theorem_rejections_credited": 0,
        },
        "timing": {
            "reconstruction_wall_seconds": reconstruction_seconds,
            "total_wall_seconds": time.perf_counter() - started,
        },
        "source_sha256": sha256(Path(__file__)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "COMPLETE",
        "output": str(args.output),
        "joint_covers": totals["joint_passing"],
        "saturating_cliques": totals["saturating_cliques"],
        "affine_histogram": dict(affine_histogram),
        "modular_units": modular_units,
        "rational_units": rational_units,
        "elapsed_seconds": report["timing"]["total_wall_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
