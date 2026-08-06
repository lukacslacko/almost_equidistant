#!/usr/bin/env python3
"""Exact prototype for the full K7 saturating-clique Schur equations.

This is deliberately a separate reference layer.  For a normalized K Gram
matrix of rank at most ``u``, a required ``u``-clique C is positive definite
and hence is a basis.  If ``p_y`` is the zero-one adjacency column from C to
an outside vertex y, then, with ``H = K[C,C]^{-1}``, Schur equality gives

    p_y.T H p_z = K[y,z] in {0,1}.

The off-diagonal equations are linear over Q in the symmetric entries of H.
This program solves them exactly.  It rejects an instance only for an exact
linear contradiction, a forced violation of a strict necessary inequality,
or (when H is unique) failure of

    H > 0,                 H^{-1} = J + positive diagonal.

It imports the existing K7 reference only to reproduce the candidate-cover
and prior-filter quantifiers.  The new algebra below is independent and uses
``fractions.Fraction`` throughout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import d6_k7_rank_reference as prior


Q = Fraction


def symmetric_coordinates(u: int) -> tuple[tuple[int, int], ...]:
    """The upper-triangular coordinates of a symmetric u by u matrix."""

    return tuple((i, j) for i in range(u) for j in range(i, u))


def bilinear_form_coefficients(
    left: Sequence[int], right: Sequence[int]
) -> tuple[Q, ...]:
    """Return coefficients of ``left.T H right`` for symmetric H."""

    if len(left) != len(right):
        raise ValueError("bilinear vectors have different dimensions")
    values = []
    for i, j in symmetric_coordinates(len(left)):
        if i == j:
            values.append(Q(left[i] * right[i]))
        else:
            values.append(Q(left[i] * right[j] + left[j] * right[i]))
    return tuple(values)


def quadratic_form_coefficients(vector: Sequence[int]) -> tuple[Q, ...]:
    return bilinear_form_coefficients(vector, vector)


def matrix_entry_coefficients(u: int, row: int, column: int) -> tuple[Q, ...]:
    if not (0 <= row < u and 0 <= column < u):
        raise ValueError("matrix entry outside H")
    if row > column:
        row, column = column, row
    return tuple(
        Q(int((i, j) == (row, column)))
        for i, j in symmetric_coordinates(u)
    )


def row_sum_coefficients(u: int, row: int) -> tuple[Q, ...]:
    return tuple(
        Q(1 if i == row or j == row else 0)
        for i, j in symmetric_coordinates(u)
    )


def total_sum_coefficients(u: int) -> tuple[Q, ...]:
    return tuple(
        Q(1 if i == j else 2) for i, j in symmetric_coordinates(u)
    )


def _fraction_json(value: Q) -> list[int]:
    return [value.numerator, value.denominator]


@dataclass(frozen=True)
class ExactAffineSystem:
    """Reduced exact affine system, retaining rational proof multipliers."""

    variables: int
    rows: tuple[tuple[Q, ...], ...]
    rhs: tuple[Q, ...]
    combinations: tuple[tuple[Q, ...], ...]
    pivots: tuple[int, ...]
    original_rows: int
    contradiction: tuple[Q, ...] | None = None
    contradiction_rhs: Q | None = None

    @property
    def consistent(self) -> bool:
        return self.contradiction is None

    @property
    def rank(self) -> int:
        return len(self.pivots)

    @property
    def unique(self) -> bool:
        return self.consistent and self.rank == self.variables

    def forced_value(
        self, functional: Sequence[Q]
    ) -> tuple[Q, tuple[Q, ...]] | None:
        """Return value and an equation combination if a form is forced."""

        if not self.consistent:
            raise ValueError("forced values are undefined for inconsistency")
        residual = [Q(value) for value in functional]
        if len(residual) != self.variables:
            raise ValueError("functional has the wrong dimension")
        value = Q(0)
        certificate = [Q(0)] * self.original_rows
        for row, rhs, combination, pivot in zip(
            self.rows, self.rhs, self.combinations, self.pivots
        ):
            factor = residual[pivot]
            if not factor:
                continue
            residual = [
                entry - factor * coefficient
                for entry, coefficient in zip(residual, row)
            ]
            value += factor * rhs
            certificate = [
                entry + factor * coefficient
                for entry, coefficient in zip(certificate, combination)
            ]
        if any(residual):
            return None
        return value, tuple(certificate)

    def unique_solution(self) -> tuple[Q, ...]:
        if not self.unique:
            raise ValueError("system does not have a unique solution")
        solution = [Q(0)] * self.variables
        for pivot, rhs in zip(self.pivots, self.rhs):
            solution[pivot] = rhs
        return tuple(solution)

    def parameterization(self) -> tuple[tuple[Q, ...], tuple[tuple[Q, ...], ...]]:
        """Return ``x = origin + sum(parameter_i * direction_i)``."""

        if not self.consistent:
            raise ValueError("inconsistent system has no parameterization")
        free = [
            column for column in range(self.variables)
            if column not in set(self.pivots)
        ]
        origin = [Q(0)] * self.variables
        directions = [[Q(0)] * self.variables for _ in free]
        for parameter, column in enumerate(free):
            directions[parameter][column] = Q(1)
        for row, rhs, pivot in zip(self.rows, self.rhs, self.pivots):
            origin[pivot] = rhs
            for parameter, column in enumerate(free):
                directions[parameter][pivot] = -row[column]
        return tuple(origin), tuple(tuple(row) for row in directions)


def exact_rref(
    coefficients: Sequence[Sequence[int | Q]],
    rhs: Sequence[int | Q],
    variables: int,
) -> ExactAffineSystem:
    """Reduce A x=b over Q and retain checkable row-combination witnesses."""

    if len(coefficients) != len(rhs):
        raise ValueError("coefficient/rhs row count mismatch")
    row_count = len(coefficients)
    augmented: list[list[Q]] = []
    transforms: list[list[Q]] = []
    for index, (row, value) in enumerate(zip(coefficients, rhs)):
        if len(row) != variables:
            raise ValueError("coefficient row has the wrong width")
        augmented.append([Q(entry) for entry in row] + [Q(value)])
        transforms.append([
            Q(int(index == column)) for column in range(row_count)
        ])

    pivot_row = 0
    pivots: list[int] = []
    for column in range(variables):
        selected = next(
            (row for row in range(pivot_row, row_count)
             if augmented[row][column]),
            None,
        )
        if selected is None:
            continue
        if selected != pivot_row:
            augmented[pivot_row], augmented[selected] = (
                augmented[selected], augmented[pivot_row]
            )
            transforms[pivot_row], transforms[selected] = (
                transforms[selected], transforms[pivot_row]
            )
        scale = augmented[pivot_row][column]
        augmented[pivot_row] = [entry / scale for entry in augmented[pivot_row]]
        transforms[pivot_row] = [entry / scale for entry in transforms[pivot_row]]
        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = augmented[row][column]
            if not factor:
                continue
            augmented[row] = [
                entry - factor * source
                for entry, source in zip(augmented[row], augmented[pivot_row])
            ]
            transforms[row] = [
                entry - factor * source
                for entry, source in zip(transforms[row], transforms[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break

    for row in range(pivot_row, row_count):
        if not any(augmented[row][:-1]) and augmented[row][-1]:
            return ExactAffineSystem(
                variables=variables,
                rows=tuple(tuple(item[:-1]) for item in augmented[:pivot_row]),
                rhs=tuple(item[-1] for item in augmented[:pivot_row]),
                combinations=tuple(tuple(item) for item in transforms[:pivot_row]),
                pivots=tuple(pivots),
                original_rows=row_count,
                contradiction=tuple(transforms[row]),
                contradiction_rhs=augmented[row][-1],
            )
    return ExactAffineSystem(
        variables=variables,
        rows=tuple(tuple(item[:-1]) for item in augmented[:pivot_row]),
        rhs=tuple(item[-1] for item in augmented[:pivot_row]),
        combinations=tuple(tuple(item) for item in transforms[:pivot_row]),
        pivots=tuple(pivots),
        original_rows=row_count,
    )


def _guided_rational_point(
    system: ExactAffineSystem,
    approximate: Sequence[float],
    max_denominator: int = 1_000_000_000,
) -> tuple[Q, ...] | None:
    """Choose rational free coordinates near a numerical affine solution."""

    if not system.consistent or len(approximate) != system.variables:
        return None
    free = [
        column for column in range(system.variables)
        if column not in set(system.pivots)
    ]
    point = [Q(0)] * system.variables
    for column in free:
        point[column] = Q(float(approximate[column])).limit_denominator(
            max_denominator
        )
    for row, rhs, pivot in zip(system.rows, system.rhs, system.pivots):
        point[pivot] = rhs - sum(
            row[column] * point[column] for column in free
        )
    return tuple(point)


def strict_motzkin_certificate(
    directions: Sequence[Sequence[Q]],
    thresholds: Sequence[Q],
    labels: Sequence[dict],
) -> dict | None:
    """Find and exactly verify a certificate that ``a_i*x>c_i`` is empty.

    Motzkin's strict alternative says infeasibility is equivalent to a
    nonzero lambda>=0 with

        sum_i lambda_i a_i = 0,    sum_i lambda_i c_i >= 0.

    We normalize ``sum lambda_i=1``.  HiGHS is used only to locate a sparse
    support; the returned certificate is reconstructed and checked over Q.
    A missed reconstruction merely leaves the instance unresolved.
    """

    if len(directions) != len(thresholds) or len(labels) != len(thresholds):
        raise ValueError("strict inequality arrays have different lengths")
    if not directions:
        return None
    dimension = len(directions[0])
    if any(len(row) != dimension for row in directions):
        raise ValueError("strict inequality directions have different widths")

    # With no affine parameter, a single false strict scalar inequality is
    # already its own rational Motzkin certificate.
    if dimension == 0:
        for index, threshold in enumerate(thresholds):
            if threshold >= 0:
                return {
                    "multipliers": [
                        {"inequality": index, "value": [1, 1]}
                    ],
                    "threshold_sum": _fraction_json(threshold),
                    "inequality_labels": list(labels),
                }
        return None

    try:
        from scipy.optimize import linprog
    except ImportError:
        return None

    count = len(directions)
    equality_rows = [
        [float(directions[index][coordinate]) for index in range(count)]
        for coordinate in range(dimension)
    ]
    equality_rows.append([1.0] * count)
    result = linprog(
        [-float(value) for value in thresholds],
        A_eq=equality_rows,
        b_eq=[0.0] * dimension + [1.0],
        bounds=[(0.0, None)] * count,
        method="highs-ds",
    )
    if not result.success or result.x is None:
        return None
    numeric_objective = sum(
        float(value) * multiplier
        for value, multiplier in zip(thresholds, result.x)
    )
    if numeric_objective < -1e-8:
        return None

    def verify(candidate: Sequence[Q]) -> dict | None:
        if len(candidate) != count or any(value < 0 for value in candidate):
            return None
        if sum(candidate) != 1:
            return None
        if any(
            sum(
                candidate[index] * directions[index][coordinate]
                for index in range(count)
            )
            for coordinate in range(dimension)
        ):
            return None
        objective = sum(
            multiplier * threshold
            for multiplier, threshold in zip(candidate, thresholds)
        )
        if objective < 0:
            return None
        return {
            "multipliers": [
                {
                    "inequality": index,
                    "value": _fraction_json(multiplier),
                }
                for index, multiplier in enumerate(candidate)
                if multiplier
            ],
            "threshold_sum": _fraction_json(objective),
            "inequality_labels": list(labels),
        }

    # HiGHS normally returns a basic solution.  Try several support cutoffs;
    # then try the full vector.  When the optimum is numerically zero, add
    # the exact objective-zero equation before rational reconstruction.
    supports: list[list[int]] = []
    for cutoff in (1e-7, 1e-9, 1e-11, 1e-13):
        support = [
            index for index, value in enumerate(result.x) if value > cutoff
        ]
        if support and support not in supports:
            supports.append(support)
    supports.append(list(range(count)))
    for support in supports:
        base_rows = [
            [directions[index][coordinate] for index in support]
            for coordinate in range(dimension)
        ]
        base_rhs = [Q(0)] * dimension
        base_rows.append([Q(1)] * len(support))
        base_rhs.append(Q(1))
        variants = [(base_rows, base_rhs)]
        if abs(numeric_objective) <= 1e-7:
            variants.insert(
                0,
                (
                    base_rows + [[thresholds[index] for index in support]],
                    base_rhs + [Q(0)],
                ),
            )
        for rows, rhs in variants:
            system = exact_rref(rows, rhs, len(support))
            if not system.consistent:
                continue
            approximate = [float(result.x[index]) for index in support]
            local = _guided_rational_point(system, approximate)
            if local is None:
                continue
            candidate = [Q(0)] * count
            for index, value in zip(support, local):
                candidate[index] = value
            certificate = verify(candidate)
            if certificate is not None:
                return certificate
    return None


def strict_affine_h_certificate(
    system: ExactAffineSystem,
    inequalities: Sequence[tuple[Sequence[Q], Q, dict]],
) -> dict | None:
    """Reduce strict H inequalities to free parameters and certify emptiness."""

    origin, basis = system.parameterization()
    directions: list[tuple[Q, ...]] = []
    thresholds: list[Q] = []
    labels: list[dict] = []
    for functional, lower_bound, label in inequalities:
        functional = tuple(Q(value) for value in functional)
        constant = sum(a * b for a, b in zip(functional, origin))
        directions.append(tuple(
            sum(a * b for a, b in zip(functional, direction))
            for direction in basis
        ))
        thresholds.append(Q(lower_bound) - constant)
        labels.append(label)
    return strict_motzkin_certificate(directions, thresholds, labels)


def symmetric_matrix(solution: Sequence[Q], u: int) -> list[list[Q]]:
    if len(solution) != u * (u + 1) // 2:
        raise ValueError("solution has the wrong symmetric dimension")
    matrix = [[Q(0) for _ in range(u)] for _ in range(u)]
    for value, (i, j) in zip(solution, symmetric_coordinates(u)):
        matrix[i][j] = matrix[j][i] = value
    return matrix


def determinant(matrix: Sequence[Sequence[Q]]) -> Q:
    """Exact determinant by fraction-preserving Gaussian elimination."""

    n = len(matrix)
    work = [[Q(value) for value in row] for row in matrix]
    if any(len(row) != n for row in work):
        raise ValueError("determinant requires a square matrix")
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
            if factor:
                for col in range(column + 1, n):
                    work[row][col] -= factor * work[column][col]
    return result


def positive_definite(matrix: Sequence[Sequence[Q]]) -> bool:
    """Exact Sylvester-criterion test for a symmetric matrix."""

    n = len(matrix)
    if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(n)):
        return False
    return all(
        determinant([list(row[:size]) for row in matrix[:size]]) > 0
        for size in range(1, n + 1)
    )


def inverse(matrix: Sequence[Sequence[Q]]) -> list[list[Q]] | None:
    n = len(matrix)
    work = [
        [Q(value) for value in row]
        + [Q(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    if any(len(row) != 2 * n for row in work):
        raise ValueError("inverse requires a square matrix")
    for column in range(n):
        selected = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if selected is None:
            return None
        work[column], work[selected] = work[selected], work[column]
        pivot = work[column][column]
        work[column] = [entry / pivot for entry in work[column]]
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    entry - factor * source
                    for entry, source in zip(work[row], work[column])
                ]
    return [row[n:] for row in work]


@dataclass(frozen=True)
class CliqueHResult:
    failed: bool
    reason: str
    u: int
    remainder_size: int
    equations: int
    rank: int
    unique: bool
    witness: dict | None = None


def _forced_witness(
    system: ExactAffineSystem,
    functional: Sequence[Q],
    kind: str,
    labels: dict,
) -> tuple[Q, dict] | None:
    forced = system.forced_value(functional)
    if forced is None:
        return None
    value, certificate = forced
    return value, {
        "failure_kind": kind,
        "forced_value": _fraction_json(value),
        "equation_combination": [_fraction_json(item) for item in certificate],
        **labels,
    }


def assess_saturating_clique(
    adj: Sequence[int], clique_mask: int
) -> CliqueHResult:
    """Apply the exact full-H screen to one saturating required clique."""

    clique = list(prior.bits(clique_mask))
    u = len(clique)
    if not u:
        raise ValueError("saturating clique must be nonempty")
    if any(
        not (adj[first] & (1 << second))
        for first, second in combinations(clique, 2)
    ):
        raise ValueError("specified basis is not a clique")
    remainder = [
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    ]
    columns = {
        vertex: tuple(int(bool(adj[vertex] & (1 << basis))) for basis in clique)
        for vertex in remainder
    }
    coefficients = []
    targets = []
    equation_labels = []
    for first, second in combinations(remainder, 2):
        coefficients.append(
            bilinear_form_coefficients(columns[first], columns[second])
        )
        targets.append(int(bool(adj[first] & (1 << second))))
        equation_labels.append([first, second])
    variables = u * (u + 1) // 2
    system = exact_rref(coefficients, targets, variables)
    witness_context = {
        "clique": clique,
        "equation_pairs": equation_labels,
    }
    common = {
        "u": u,
        "remainder_size": len(remainder),
        "equations": len(coefficients),
        "rank": system.rank,
        "unique": system.unique,
    }
    if not system.consistent:
        assert system.contradiction is not None
        assert system.contradiction_rhs is not None
        return CliqueHResult(
            failed=True,
            reason="linear_inconsistency",
            witness={
                "clique": clique,
                "equation_pairs": equation_labels,
                "equation_combination": [
                    _fraction_json(item) for item in system.contradiction
                ],
                "contradiction_rhs": _fraction_json(
                    system.contradiction_rhs
                ),
            },
            **common,
        )

    strict_inequalities: list[tuple[tuple[Q, ...], Q, dict]] = []
    for i in range(u):
        strict_inequalities.append((
            matrix_entry_coefficients(u, i, i),
            Q(0),
            {"kind": "H_diagonal_positive", "basis_indices": [i]},
        ))
    for i, j in combinations(range(u), 2):
        strict_inequalities.append((
            tuple(-value for value in matrix_entry_coefficients(u, i, j)),
            Q(0),
            {"kind": "H_offdiagonal_negative", "basis_indices": [i, j]},
        ))
    for i in range(u):
        strict_inequalities.append((
            row_sum_coefficients(u, i),
            Q(0),
            {"kind": "H_row_sum_positive", "basis_indices": [i]},
        ))
    total_form = total_sum_coefficients(u)
    strict_inequalities.extend((
        (
            total_form,
            Q(0),
            {"kind": "H_total_sum_positive"},
        ),
        (
            tuple(-value for value in total_form),
            Q(-1),
            {"kind": "H_total_sum_below_one"},
        ),
    ))
    for vertex in remainder:
        strict_inequalities.append((
            quadratic_form_coefficients(columns[vertex]),
            Q(1),
            {"kind": "outside_diagonal_above_one", "vertices": [vertex]},
        ))
    for first, second in combinations(remainder, 2):
        difference = tuple(
            a - b for a, b in zip(columns[first], columns[second])
        )
        if any(difference):
            strict_inequalities.append((
                quadratic_form_coefficients(difference),
                Q(0),
                {
                    "kind": "column_difference_form_positive",
                    "vertices": [first, second],
                },
            ))
    strict_certificate = strict_affine_h_certificate(
        system, strict_inequalities
    )
    if strict_certificate is not None:
        origin, basis = system.parameterization()
        return CliqueHResult(
            failed=True,
            reason="strict_affine_H_infeasible",
            witness={
                **witness_context,
                "certificate": strict_certificate,
                "affine_origin": [_fraction_json(value) for value in origin],
                "affine_directions": [
                    [_fraction_json(value) for value in direction]
                    for direction in basis
                ],
            },
            **common,
        )

    # H=(J+positive diagonal)^{-1} has positive diagonal, strictly negative
    # off-diagonal, positive row sums, and total entry sum strictly between
    # zero and one.  Test only a sign when its linear form is forced by the
    # Schur equations; otherwise leave the cover unresolved.
    for i in range(u):
        check = _forced_witness(
            system,
            matrix_entry_coefficients(u, i, i),
            "forced_nonpositive_H_diagonal",
            {**witness_context, "basis_indices": [i]},
        )
        if check is not None and check[0] <= 0:
            return CliqueHResult(True, check[1]["failure_kind"], witness=check[1], **common)
    for i, j in combinations(range(u), 2):
        check = _forced_witness(
            system,
            matrix_entry_coefficients(u, i, j),
            "forced_nonnegative_H_offdiagonal",
            {**witness_context, "basis_indices": [i, j]},
        )
        if check is not None and check[0] >= 0:
            return CliqueHResult(True, check[1]["failure_kind"], witness=check[1], **common)
    for i in range(u):
        check = _forced_witness(
            system,
            row_sum_coefficients(u, i),
            "forced_nonpositive_H_row_sum",
            {**witness_context, "basis_indices": [i]},
        )
        if check is not None and check[0] <= 0:
            return CliqueHResult(True, check[1]["failure_kind"], witness=check[1], **common)
    check = _forced_witness(
        system,
        total_sum_coefficients(u),
        "forced_invalid_H_total_sum",
        witness_context,
    )
    if check is not None and not (0 < check[0] < 1):
        return CliqueHResult(True, check[1]["failure_kind"], witness=check[1], **common)

    # The omitted Schur diagonal is K_yy=1+r_y^2>1.  Also, distinct p_y
    # columns give nonzero vector differences whose H-quadratic form is >0.
    for vertex in remainder:
        check = _forced_witness(
            system,
            quadratic_form_coefficients(columns[vertex]),
            "forced_outside_diagonal_at_most_one",
            {**witness_context, "vertices": [vertex]},
        )
        if check is not None and check[0] <= 1:
            return CliqueHResult(True, check[1]["failure_kind"], witness=check[1], **common)
    for first, second in combinations(remainder, 2):
        difference = tuple(
            a - b for a, b in zip(columns[first], columns[second])
        )
        if not any(difference):
            # The earlier mask rule rejects this case.  Retain a local sound
            # failure so this routine is safe to call independently.
            return CliqueHResult(
                True,
                "duplicate_basis_columns",
                witness={
                    **witness_context,
                    "vertices": [first, second],
                },
                **common,
            )
        check = _forced_witness(
            system,
            quadratic_form_coefficients(difference),
            "forced_nonpositive_column_difference_form",
            {**witness_context, "vertices": [first, second]},
        )
        if check is not None and check[0] <= 0:
            return CliqueHResult(True, check[1]["failure_kind"], witness=check[1], **common)

    if not system.unique:
        return CliqueHResult(False, "underdetermined", witness=None, **common)

    hmatrix = symmetric_matrix(system.unique_solution(), u)
    if not positive_definite(hmatrix):
        return CliqueHResult(
            True,
            "unique_H_not_positive_definite",
            witness={
                "H": [[_fraction_json(item) for item in row] for row in hmatrix]
            },
            **common,
        )
    kmatrix = inverse(hmatrix)
    assert kmatrix is not None
    for i, j in combinations(range(u), 2):
        if kmatrix[i][j] != 1:
            return CliqueHResult(
                True,
                "unique_H_inverse_offdiagonal_not_one",
                witness={
                    "basis_indices": [i, j],
                    "inverse_entry": _fraction_json(kmatrix[i][j]),
                },
                **common,
            )
    for i in range(u):
        if kmatrix[i][i] <= 1:
            return CliqueHResult(
                True,
                "unique_H_inverse_diagonal_at_most_one",
                witness={
                    "basis_indices": [i],
                    "inverse_entry": _fraction_json(kmatrix[i][i]),
                },
                **common,
            )
    return CliqueHResult(False, "unique_valid_special_H", witness=None, **common)


@dataclass
class GraphHResult:
    index: int | None
    stratum: str | None
    applicable: bool
    seeds: int
    covers: int
    prior_passing_covers: int
    h_passing_covers: int
    saturating_cliques: int
    prior_rejected: bool
    combined_rejected: bool
    marginal_rejected: bool
    clique_reasons: dict[str, int] = field(default_factory=dict)
    affine_rank_histogram: dict[str, int] = field(default_factory=dict)
    first_failure: dict | None = None
    elapsed_seconds: float = 0.0


def analyze_graph(graph: dict) -> GraphHResult:
    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    seed_count = cover_count = prior_passing = h_passing = 0
    saturating_count = 0
    prior_rejected = combined_rejected = False
    reasons: Counter[str] = Counter()
    ranks: Counter[str] = Counter()
    first_failure = None
    for seed_mask in prior.clique_masks(adj, 7):
        seed_count += 1
        _, outside, defects, ladj, eligible = prior.seed_instance(adj, seed_mask)
        covers = prior.eligible_covers(ladj, eligible)
        cover_count += len(covers)
        total_term_rank = prior.matching_size(defects)
        seed_prior_passes = 0
        seed_h_passes = 0
        for zmask in covers:
            analysis = prior.analyze_cover(
                adj,
                outside,
                defects,
                zmask,
                support_solver,
                zero_forcing,
                total_term_rank,
                clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            seed_prior_passes += 1
            prior_passing += 1
            nvertices = [
                outside[i] for i in range(len(outside))
                if not (zmask & (1 << i))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            clique_results = [
                assess_saturating_clique(graph_n, clique_mask)
                for clique_mask in prior.clique_masks(
                    graph_n, analysis.k_rank_upper
                )
            ] if analysis.k_rank_upper else []
            saturating_count += len(clique_results)
            for result in clique_results:
                reasons[result.reason] += 1
                ranks[f"u{result.u}_rank{result.rank}_vars{result.u*(result.u+1)//2}"] += 1
            failed = next((entry for entry in clique_results if entry.failed), None)
            if failed is not None:
                if first_failure is None:
                    first_failure = {
                        "seed": list(prior.bits(seed_mask)),
                        "zmask": zmask,
                        "reason": failed.reason,
                        "witness": failed.witness,
                    }
                continue
            seed_h_passes += 1
            h_passing += 1
        if not seed_prior_passes:
            prior_rejected = True
        if not seed_h_passes:
            combined_rejected = True
    return GraphHResult(
        index=graph.get("index"),
        stratum=graph.get("stratum"),
        applicable=bool(seed_count),
        seeds=seed_count,
        covers=cover_count,
        prior_passing_covers=prior_passing,
        h_passing_covers=h_passing,
        saturating_cliques=saturating_count,
        prior_rejected=prior_rejected,
        combined_rejected=combined_rejected,
        marginal_rejected=combined_rejected and not prior_rejected,
        clique_reasons=dict(reasons),
        affine_rank_histogram=dict(ranks),
        first_failure=first_failure,
        elapsed_seconds=time.perf_counter() - started,
    )


def _analyze_graph_worker(graph: dict) -> GraphHResult:
    return analyze_graph(graph)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_report(sample: dict, workers: int) -> dict:
    started = time.perf_counter()
    if workers == 1:
        results = [analyze_graph(graph) for graph in sample["graphs"]]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(
                executor.map(_analyze_graph_worker, sample["graphs"], chunksize=1)
            )
    reasons: Counter[str] = Counter()
    ranks: Counter[str] = Counter()
    for result in results:
        reasons.update(result.clique_reasons)
        ranks.update(result.affine_rank_histogram)
    indices = lambda predicate: sorted(
        result.index for result in results
        if result.index is not None and predicate(result)
    )
    return {
        "schema": 1,
        "description": (
            "Exact rational affine-H screen after the K7 enhanced joint "
            "cover filters; all candidate-nonedges remain optional zeros."
        ),
        "graphs": len(results),
        "workers": workers,
        "applicable_graphs": sum(result.applicable for result in results),
        "seeds": sum(result.seeds for result in results),
        "covers": sum(result.covers for result in results),
        "prior_passing_covers": sum(result.prior_passing_covers for result in results),
        "h_passing_covers": sum(result.h_passing_covers for result in results),
        "saturating_cliques": sum(result.saturating_cliques for result in results),
        "clique_result_reasons": dict(reasons),
        "affine_rank_histogram": dict(ranks),
        "graph_counts": {
            "prior_rejected": sum(result.prior_rejected for result in results),
            "combined_rejected": sum(result.combined_rejected for result in results),
            "marginal_h_rejected": sum(result.marginal_rejected for result in results),
            "combined_survivors": sum(
                result.applicable and not result.combined_rejected
                for result in results
            ),
        },
        "graph_decisions": {
            "prior_rejected": indices(lambda result: result.prior_rejected),
            "combined_rejected": indices(lambda result: result.combined_rejected),
            "marginal_h_rejected": indices(lambda result: result.marginal_rejected),
            "combined_survivors": indices(
                lambda result: result.applicable and not result.combined_rejected
            ),
        },
        "runtime": {
            "wall_seconds": time.perf_counter() - started,
            "sum_graph_seconds": sum(result.elapsed_seconds for result in results),
            "maximum_graph_seconds": max(
                (result.elapsed_seconds for result in results), default=0.0
            ),
        },
        "per_graph": [asdict(result) for result in results],
    }


def self_test() -> None:
    # Exact RREF: x+y=3, x-y=1, with checkable forced values.
    system = exact_rref([[1, 1], [1, -1]], [3, 1], 2)
    assert system.unique and system.unique_solution() == (Q(2), Q(1))
    forced = system.forced_value([2, 3])
    assert forced is not None and forced[0] == 7
    inconsistent = exact_rref([[1], [1]], [0, 1], 1)
    assert not inconsistent.consistent
    assert inconsistent.contradiction is not None
    assert sum(
        coefficient * value
        for coefficient, value in zip(inconsistent.contradiction, [0, 1])
    ) == inconsistent.contradiction_rhs

    # H=(J+diag(1,2))^-1 is positive definite and has the special inverse.
    k = [[Q(2), Q(1)], [Q(1), Q(3)]]
    h = inverse(k)
    assert h is not None and positive_definite(h)
    recovered = inverse(h)
    assert recovered == k
    coordinates = tuple(
        h[i][j] for i, j in symmetric_coordinates(2)
    )
    assert symmetric_matrix(coordinates, 2) == h

    # Exact strict alternative.  x>0 and -x>0 is impossible with the
    # rational certificate (1/2,1/2); x>0 and -x>-1 is feasible.
    labels = [{"kind": "lower"}, {"kind": "upper"}]
    certificate = strict_motzkin_certificate(
        [(Q(1),), (Q(-1),)], [Q(0), Q(0)], labels
    )
    assert certificate is not None
    assert certificate["threshold_sum"] == [0, 1]
    assert strict_motzkin_certificate(
        [(Q(1),), (Q(-1),)], [Q(0), Q(-1)], labels
    ) is None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", nargs="?", type=Path)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--prior-survivors-from",
        type=Path,
        help=(
            "select the union of combined_survivors and marginal_h_rejected "
            "from an earlier special-H report"
        ),
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        print("self-test passed")
        if args.sample is None:
            return
    if args.sample is None:
        parser.error("sample is required unless only --self-test is requested")
    sample = json.loads(args.sample.read_text(encoding="utf-8"))
    selected_indices = None
    if args.prior_survivors_from:
        selection_report = json.loads(
            args.prior_survivors_from.read_text(encoding="utf-8")
        )
        decisions = selection_report["graph_decisions"]
        selected_indices = set(decisions["combined_survivors"])
        selected_indices.update(decisions.get("marginal_h_rejected", []))
        sample = {
            **sample,
            "graphs": [
                graph for graph in sample["graphs"]
                if graph.get("index") in selected_indices
            ],
        }
    report = sample_report(sample, args.workers)
    report["input_sample"] = {
        "file": args.sample.name,
        "sha256": file_sha256(args.sample),
    }
    source = Path(__file__)
    report["source"] = {"file": source.name, "sha256": file_sha256(source)}
    if args.prior_survivors_from:
        report["selection"] = {
            "file": args.prior_survivors_from.name,
            "sha256": file_sha256(args.prior_survivors_from),
            "selected_indices": sorted(selected_indices or ()),
        }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
