#!/usr/bin/env python3
"""Pilot exact degree-four preordering duals for K7 Schur systems.

This is deliberately separate from the production degree-one campaign.  It
adds the strict outside-diagonal inequalities to a truncated positive
preordering.  Numerical linear programming only locates a rational identity;
``verify_raw_certificate`` checks the complete identity over ``Fraction``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import time
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import d6_k7_positive_polynomial_dual as degree_one
import d6_k7_rank_reference as prior
import d6_k7_special_h_reference as strict_h


Q = Fraction
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Q]


def fraction_json(value: Q) -> list[int]:
    return [value.numerator, value.denominator]


def parse_fraction(value: Sequence[int]) -> Q:
    if len(value) != 2:
        raise ValueError("fraction must have two entries")
    return Q(int(value[0]), int(value[1]))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def monomials(variables: int, maximum_degree: int) -> tuple[Exponent, ...]:
    if maximum_degree < 0:
        return ()
    return tuple(
        exponent
        for exponent in itertools.product(
            range(maximum_degree + 1), repeat=variables
        )
        if sum(exponent) <= maximum_degree
    )


def add_term(polynomial: Polynomial, exponent: Exponent, value: Q) -> None:
    if not value:
        return
    updated = polynomial.get(exponent, Q(0)) + value
    if updated:
        polynomial[exponent] = updated
    else:
        polynomial.pop(exponent, None)


def shifted(polynomial: Polynomial, multiplier: Exponent) -> Polynomial:
    return {
        tuple(a + b for a, b in zip(exponent, multiplier)): coefficient
        for exponent, coefficient in polynomial.items()
    }


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = tuple(
                a + b for a, b in zip(left_exponent, right_exponent)
            )
            add_term(output, exponent, left_value * right_value)
    return output


def outside_diagonal_polynomial(mask: Sequence[int]) -> Polynomial:
    """Numerator of ``p.T H p - 1`` for ``H=(J+diag(e))^-1``."""

    variables = len(mask)
    zero = (0,) * variables
    polynomial: Polynomial = {zero: Q(-1)}
    # (1 + sum x_i) w(A).
    for coordinate, present in enumerate(mask):
        if not present:
            continue
        linear = list(zero)
        linear[coordinate] = 1
        add_term(polynomial, tuple(linear), Q(1))
        for other in range(variables):
            quadratic = linear.copy()
            quadratic[other] += 1
            add_term(polynomial, tuple(quadratic), Q(1))
    # -w(A)^2.
    for first, first_present in enumerate(mask):
        if not first_present:
            continue
        for second, second_present in enumerate(mask):
            if not second_present:
                continue
            quadratic = [0] * variables
            quadratic[first] += 1
            quadratic[second] += 1
            add_term(polynomial, tuple(quadratic), Q(-1))
    # -(1 + sum x_i).
    for coordinate in range(variables):
        linear = list(zero)
        linear[coordinate] = 1
        add_term(polynomial, tuple(linear), Q(-1))
    return polynomial


@dataclass(frozen=True)
class Generator:
    vertices: tuple[int, ...]
    polynomial: Polynomial


def preordering_generators(
    system: degree_one.SchurPolynomialSystem,
    maximum_order: int,
) -> tuple[Generator, ...]:
    base = tuple(
        Generator(
            vertices=(vertex,),
            polynomial=outside_diagonal_polynomial(mask),
        )
        for vertex, mask in zip(system.remainder, system.masks)
    )
    output = list(base)
    for order in range(2, maximum_order + 1):
        for selected in itertools.combinations(base, order):
            polynomial: Polynomial = {(0,) * system.u: Q(1)}
            for generator in selected:
                polynomial = multiply(polynomial, generator.polynomial)
            output.append(Generator(
                vertices=tuple(item.vertices[0] for item in selected),
                polynomial=polynomial,
            ))
    return tuple(output)


def _guided_rational_point(
    system: strict_h.ExactAffineSystem,
    approximate: Sequence[float],
    maximum_denominator: int,
) -> tuple[Q, ...]:
    pivots = set(system.pivots)
    free = [
        column for column in range(system.variables) if column not in pivots
    ]
    point = [Q(0)] * system.variables
    for column in free:
        point[column] = Q(float(approximate[column])).limit_denominator(
            maximum_denominator
        )
    for row, rhs, pivot in zip(system.rows, system.rhs, system.pivots):
        point[pivot] = rhs - sum(
            row[column] * point[column] for column in free
        )
    return tuple(point)


def _coefficient_data(
    equations: Sequence[Polynomial],
    generators: Sequence[Generator],
    variables: int,
    output_degree: int,
) -> tuple[
    tuple[Exponent, ...],
    tuple[Exponent, ...],
    tuple[tuple[Exponent, ...], ...],
    tuple[tuple[int, ...], ...],
    int,
]:
    equation_monomials = monomials(variables, output_degree - 2)
    generator_monomials = tuple(
        monomials(
            variables,
            output_degree - max(sum(exponent) for exponent in generator.polynomial),
        )
        for generator in generators
    )
    outputs = monomials(variables, output_degree)
    output_index = {exponent: index for index, exponent in enumerate(outputs)}
    equation_columns = len(equations) * len(equation_monomials)
    columns = equation_columns + sum(map(len, generator_monomials))
    matrix = [[0] * columns for _ in outputs]
    for equation_index, equation in enumerate(equations):
        for monomial_index, monomial in enumerate(equation_monomials):
            column = equation_index * len(equation_monomials) + monomial_index
            for exponent, coefficient in shifted(equation, monomial).items():
                matrix[output_index[exponent]][column] += int(coefficient)
    column = equation_columns
    for generator, multiplier_monomials in zip(
        generators, generator_monomials
    ):
        for monomial in multiplier_monomials:
            # Identity convention: sum q_i f_i = P + sum r_j g_j.
            for exponent, coefficient in shifted(
                generator.polynomial, monomial
            ).items():
                matrix[output_index[exponent]][column] -= int(coefficient)
            column += 1
    return (
        equation_monomials,
        outputs,
        generator_monomials,
        tuple(tuple(row) for row in matrix),
        equation_columns,
    )


def verify_raw_certificate(
    equations: Sequence[Polynomial],
    generators: Sequence[Generator],
    variables: int,
    certificate: dict,
) -> Polynomial:
    """Exactly verify ``sum q f = P + sum r g`` and positivity data."""

    left: Polynomial = {}
    for item in certificate.get("equation_multipliers", []):
        equation_index = int(item["equation"])
        if not 0 <= equation_index < len(equations):
            raise ValueError("equation index out of range")
        monomial = tuple(int(value) for value in item["monomial"])
        if len(monomial) != variables or any(value < 0 for value in monomial):
            raise ValueError("invalid equation-multiplier monomial")
        coefficient = parse_fraction(item["coefficient"])
        for exponent, value in shifted(
            equations[equation_index], monomial
        ).items():
            add_term(left, exponent, coefficient * value)

    generator_by_vertices = {
        generator.vertices: generator for generator in generators
    }
    right: Polynomial = {}
    positive_mass = Q(0)
    for item in certificate.get("generator_multipliers", []):
        vertices = tuple(int(value) for value in item["vertices"])
        generator = generator_by_vertices.get(vertices)
        if generator is None:
            raise ValueError("unknown preordering generator")
        monomial = tuple(int(value) for value in item["monomial"])
        if len(monomial) != variables or any(value < 0 for value in monomial):
            raise ValueError("invalid generator-multiplier monomial")
        coefficient = parse_fraction(item["coefficient"])
        if coefficient <= 0:
            raise ValueError("generator multiplier is not positive")
        positive_mass += coefficient
        for exponent, value in shifted(generator.polynomial, monomial).items():
            add_term(right, exponent, coefficient * value)

    positive: Polynomial = {}
    for item in certificate.get("positive_polynomial", []):
        exponent = tuple(int(value) for value in item["monomial"])
        if len(exponent) != variables or any(value < 0 for value in exponent):
            raise ValueError("invalid positive-polynomial monomial")
        coefficient = parse_fraction(item["coefficient"])
        if coefficient <= 0:
            raise ValueError("positive-polynomial coefficient is not positive")
        positive_mass += coefficient
        add_term(positive, exponent, coefficient)
        add_term(right, exponent, coefficient)
    if positive_mass <= 0:
        raise ValueError("certificate has no strictly positive term")
    if left != right:
        raise ValueError("preordering polynomial identity does not check")
    return positive


def find_raw_certificate(
    equations: Sequence[Polynomial],
    generators: Sequence[Generator],
    variables: int,
    output_degree: int = 4,
) -> dict | None:
    """Locate and exactly reconstruct a truncated preordering identity."""

    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        return None
    (
        equation_monomials,
        outputs,
        generator_monomials,
        matrix,
        equation_columns,
    ) = _coefficient_data(
        equations, generators, variables, output_degree
    )
    if not matrix or not matrix[0]:
        return None
    numeric = np.asarray(matrix, dtype=float)
    generator_columns = numeric.shape[1] - equation_columns
    normalization = numeric.sum(axis=0)
    if generator_columns:
        normalization[equation_columns:] += 1.0
    bounds = (
        [(None, None)] * equation_columns
        + [(0.0, None)] * generator_columns
    )
    result = linprog(
        np.zeros(numeric.shape[1]),
        A_ub=-numeric,
        b_ub=np.zeros(numeric.shape[0]),
        A_eq=[normalization],
        b_eq=[1.0],
        bounds=bounds,
        method="highs-ds",
    )
    if not result.success or result.x is None:
        return None
    approximate_output = numeric @ result.x
    exact_normalization = [sum(row[column] for row in matrix) for column in range(len(matrix[0]))]
    for column in range(equation_columns, len(exact_normalization)):
        exact_normalization[column] += 1

    for cutoff in (1e-7, 1e-9, 1e-11, 1e-13):
        zero_rows = [
            row for row, value in zip(matrix, approximate_output)
            if value <= cutoff
        ]
        zero_generator_columns = [
            column
            for column in range(equation_columns, numeric.shape[1])
            if result.x[column] <= cutoff
        ]
        coordinate_rows = []
        for column in zero_generator_columns:
            row = [0] * numeric.shape[1]
            row[column] = 1
            coordinate_rows.append(row)
        affine = strict_h.exact_rref(
            zero_rows + coordinate_rows + [exact_normalization],
            [0] * (len(zero_rows) + len(coordinate_rows)) + [1],
            numeric.shape[1],
        )
        if not affine.consistent:
            continue
        for maximum_denominator in (1_000, 1_000_000, 1_000_000_000):
            rational = _guided_rational_point(
                affine, result.x, maximum_denominator
            )
            generator_values = rational[equation_columns:]
            if any(value < 0 for value in generator_values):
                continue
            positive = [
                sum(Q(coefficient) * value for coefficient, value in zip(row, rational))
                for row in matrix
            ]
            if any(value < 0 for value in positive):
                continue
            if sum(positive) + sum(generator_values) != 1:
                continue
            equation_items = []
            for column, value in enumerate(rational[:equation_columns]):
                if value:
                    equation_items.append({
                        "equation": column // len(equation_monomials),
                        "monomial": list(
                            equation_monomials[column % len(equation_monomials)]
                        ),
                        "coefficient": fraction_json(value),
                    })
            generator_items = []
            column = equation_columns
            for generator, multiplier_monomials in zip(
                generators, generator_monomials
            ):
                for monomial in multiplier_monomials:
                    value = rational[column]
                    if value:
                        generator_items.append({
                            "vertices": list(generator.vertices),
                            "monomial": list(monomial),
                            "coefficient": fraction_json(value),
                        })
                    column += 1
            certificate = {
                "schema": 1,
                "kind": "positive_orthant_outside_diagonal_preordering",
                "output_degree": output_degree,
                "maximum_generator_order": max(
                    (len(generator.vertices) for generator in generators),
                    default=0,
                ),
                "equation_multipliers": equation_items,
                "generator_multipliers": generator_items,
                "positive_polynomial": [
                    {
                        "monomial": list(outputs[row]),
                        "coefficient": fraction_json(value),
                    }
                    for row, value in enumerate(positive) if value
                ],
            }
            verify_raw_certificate(
                equations, generators, variables, certificate
            )
            return certificate
    return None


def find_clique_certificate(
    adj: Sequence[int],
    clique_mask: int,
    output_degree: int = 4,
    maximum_generator_order: int = 2,
) -> dict | None:
    system = degree_one.schur_polynomial_system(adj, clique_mask)
    generators = preordering_generators(system, maximum_generator_order)
    certificate = find_raw_certificate(
        system.equations, generators, system.u, output_degree
    )
    if certificate is None:
        return None
    certificate.update({
        "clique": list(system.clique),
        "remainder": list(system.remainder),
        "basis_neighbour_masks": [list(mask) for mask in system.masks],
        "equation_pairs": [list(pair) for pair in system.equation_pairs],
        "equation_targets": list(system.equation_targets),
    })
    verify_clique_certificate(adj, clique_mask, certificate)
    return certificate


def verify_clique_certificate(
    adj: Sequence[int], clique_mask: int, certificate: dict
) -> None:
    system = degree_one.schur_polynomial_system(adj, clique_mask)
    if certificate.get("clique") != list(system.clique):
        raise ValueError("certificate clique mismatch")
    if certificate.get("remainder") != list(system.remainder):
        raise ValueError("certificate remainder mismatch")
    if certificate.get("basis_neighbour_masks") != [
        list(mask) for mask in system.masks
    ]:
        raise ValueError("certificate basis masks mismatch")
    if certificate.get("equation_pairs") != [
        list(pair) for pair in system.equation_pairs
    ]:
        raise ValueError("certificate equation pairs mismatch")
    if certificate.get("equation_targets") != list(system.equation_targets):
        raise ValueError("certificate equation targets mismatch")
    generators = preordering_generators(
        system, int(certificate["maximum_generator_order"])
    )
    verify_raw_certificate(
        system.equations, generators, system.u, certificate
    )


@dataclass
class PilotGraphResult:
    index: int
    seeds: int = 0
    covers: int = 0
    enhanced_passing_covers: int = 0
    strict_h_passing_covers: int = 0
    covers_without_saturating_clique: int = 0
    baseline_degree1_failed_covers: int = 0
    pure_degree2_new_failed_covers: int = 0
    preordering_new_failed_covers: int = 0
    final_passing_covers: int = 0
    baseline_degree1_rejected: bool = False
    pure_degree2_rejected: bool = False
    preordering_rejected: bool = False
    certificates: list[dict] = field(default_factory=list)
    elapsed_seconds: float = 0.0


def analyze_graph(
    graph: dict,
    output_degree: int = 4,
    maximum_generator_order: int = 2,
) -> PilotGraphResult:
    """Pilot the hierarchy with complete seed/cover quantifiers for a graph."""

    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    result = PilotGraphResult(index=int(graph["index"]))
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    baseline_rejected = pure_rejected = preordering_rejected = False
    for seed_mask in prior.clique_masks(adj, 7):
        result.seeds += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(
            adj, seed_mask
        )
        total_term_rank = prior.matching_size(defects)
        baseline_seed_passes = pure_seed_passes = preordering_seed_passes = 0
        for zmask in prior.eligible_covers(ladj, eligible):
            result.covers += 1
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            result.enhanced_passing_covers += 1
            nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            clique_masks = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            if any(
                strict_h.assess_saturating_clique(graph_n, mask).failed
                for mask in clique_masks
            ):
                continue
            result.strict_h_passing_covers += 1
            if not clique_masks:
                result.covers_without_saturating_clique += 1
                baseline_seed_passes += 1
                pure_seed_passes += 1
                preordering_seed_passes += 1
                result.final_passing_covers += 1
                continue
            baseline_certificate = next((
                certificate
                for mask in clique_masks
                if (certificate := degree_one.find_clique_certificate(
                    graph_n, mask, multiplier_degree=1
                )) is not None
            ), None)
            if baseline_certificate is not None:
                result.baseline_degree1_failed_covers += 1
                continue
            baseline_seed_passes += 1
            pure_certificate = next((
                certificate
                for mask in clique_masks
                if (certificate := degree_one.find_clique_certificate(
                    graph_n, mask, multiplier_degree=2
                )) is not None
            ), None)
            if pure_certificate is not None:
                result.pure_degree2_new_failed_covers += 1
                result.certificates.append({
                    "kind": "pure_degree2",
                    "seed": seed,
                    "zmask": zmask,
                    "nvertices": nvertices,
                    "certificate": pure_certificate,
                })
                continue
            pure_seed_passes += 1
            preordering_certificate = next((
                certificate
                for mask in clique_masks
                if (certificate := find_clique_certificate(
                    graph_n, mask, output_degree, maximum_generator_order
                )) is not None
            ), None)
            if preordering_certificate is not None:
                result.preordering_new_failed_covers += 1
                result.certificates.append({
                    "kind": "outside_diagonal_preordering",
                    "seed": seed,
                    "zmask": zmask,
                    "nvertices": nvertices,
                    "certificate": preordering_certificate,
                })
                continue
            preordering_seed_passes += 1
            result.final_passing_covers += 1
        if not baseline_seed_passes:
            baseline_rejected = True
        if not pure_seed_passes:
            pure_rejected = True
        if not preordering_seed_passes:
            preordering_rejected = True
    result.baseline_degree1_rejected = baseline_rejected
    result.pure_degree2_rejected = pure_rejected
    result.preordering_rejected = preordering_rejected
    result.elapsed_seconds = time.perf_counter() - started
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=Path("d6_k7_rank_sample.json"))
    parser.add_argument("--indices", type=int, nargs="+", required=True)
    parser.add_argument("--output-degree", type=int, default=4)
    parser.add_argument("--maximum-generator-order", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.sample.read_text(encoding="utf-8"))
    requested = set(args.indices)
    graphs = [
        graph for graph in payload["graphs"] if int(graph["index"]) in requested
    ]
    if {int(graph["index"]) for graph in graphs} != requested:
        raise ValueError("one or more requested graph indices are absent")
    started = time.perf_counter()
    results = [
        analyze_graph(
            graph, args.output_degree, args.maximum_generator_order
        )
        for graph in sorted(graphs, key=lambda item: int(item["index"]))
    ]
    report = {
        "schema": 1,
        "kind": "d6_k7_dual_degree2_pilot",
        "input": {
            "path": args.sample.name,
            "sha256": file_sha256(args.sample),
            "indices": sorted(requested),
        },
        "parameters": {
            "pure_equation_multiplier_degree": 2,
            "preordering_output_degree": args.output_degree,
            "preordering_maximum_generator_order": args.maximum_generator_order,
        },
        "counts": {
            "graphs": len(results),
            "strict_h_passing_covers": sum(
                item.strict_h_passing_covers for item in results
            ),
            "covers_without_saturating_clique": sum(
                item.covers_without_saturating_clique for item in results
            ),
            "covers_with_saturating_clique": sum(
                item.strict_h_passing_covers
                - item.covers_without_saturating_clique
                for item in results
            ),
            "baseline_degree1_failed_covers": sum(
                item.baseline_degree1_failed_covers for item in results
            ),
            "baseline_degree1_surviving_saturating_clique_covers": sum(
                item.strict_h_passing_covers
                - item.covers_without_saturating_clique
                - item.baseline_degree1_failed_covers
                for item in results
            ),
            "pure_degree2_new_failed_covers": sum(
                item.pure_degree2_new_failed_covers for item in results
            ),
            "preordering_new_failed_covers": sum(
                item.preordering_new_failed_covers for item in results
            ),
            "pure_degree2_new_graph_rejections": sum(
                item.pure_degree2_rejected and not item.baseline_degree1_rejected
                for item in results
            ),
            "preordering_new_graph_rejections": sum(
                item.preordering_rejected and not item.pure_degree2_rejected
                for item in results
            ),
            "final_passing_covers": sum(
                item.final_passing_covers for item in results
            ),
        },
        "runtime": {"wall_seconds": time.perf_counter() - started},
        "per_graph": [asdict(item) for item in results],
    }
    report["source"] = {
        "path": Path(__file__).name,
        "sha256": file_sha256(Path(__file__)),
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
