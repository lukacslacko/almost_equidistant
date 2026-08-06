#!/usr/bin/env python3
"""Exact positive-polynomial dual certificates for K7 Schur systems.

For a saturating clique ``C`` in the normalized K matrix, write

    K[C,C] = J + diag(e_i),  x_i = 1/e_i > 0,
    T = 1 + sum_i x_i.

If an outside basis-neighbour mask is ``A``, put
``w(A)=sum_{i in A} x_i``.  The Schur equation for two outside masks is

    f_AB = T*w(A intersect B) - w(A)*w(B) - k_AB*T = 0,

where ``k_AB`` is zero or one.  A certificate consists of rational
polynomial multipliers ``q_AB`` for which

    P = sum_AB q_AB f_AB

is nonzero and every coefficient of ``P`` is nonnegative.  Since every
``x_i`` is strictly positive, ``P(x)>0``; the Schur equations instead give
``P(x)=0``.  This is an exact Farkas/positive-orthant certificate.

SciPy/HiGHS is used only to locate a certificate.  The returned sparse
rational identity is reconstructed and verified with ``Fraction``.  A
failed numerical search leaves the clique unresolved.
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

import d6_k7_rank_reference as prior
import d6_k7_special_h_reference as strict_h


Q = Fraction
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Q]


def fraction_json(value: Q) -> list[int]:
    return [value.numerator, value.denominator]


def parse_fraction(value: Sequence[int]) -> Q:
    if len(value) != 2:
        raise ValueError("fraction must have numerator and denominator")
    return Q(int(value[0]), int(value[1]))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def monomials(variables: int, maximum_degree: int) -> tuple[Exponent, ...]:
    """All monomials of total degree at most ``maximum_degree``."""

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
        tuple(left + right for left, right in zip(exponent, multiplier)): value
        for exponent, value in polynomial.items()
    }


@dataclass(frozen=True)
class SchurPolynomialSystem:
    u: int
    clique: tuple[int, ...]
    remainder: tuple[int, ...]
    masks: tuple[tuple[int, ...], ...]
    equation_pairs: tuple[tuple[int, int], ...]
    equation_targets: tuple[int, ...]
    equations: tuple[Polynomial, ...]


def schur_polynomial_system(
    adj: Sequence[int], clique_mask: int
) -> SchurPolynomialSystem:
    """Build the exact Sherman--Morrison polynomial equations."""

    clique = tuple(prior.bits(clique_mask))
    u = len(clique)
    if not u:
        raise ValueError("basis clique is empty")
    if any(
        not (adj[first] & (1 << second))
        for first, second in itertools.combinations(clique, 2)
    ):
        raise ValueError("basis mask is not a required clique")
    remainder = tuple(
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    )
    masks = tuple(
        tuple(int(bool(adj[vertex] & (1 << basis))) for basis in clique)
        for vertex in remainder
    )
    zero = (0,) * u
    pairs: list[tuple[int, int]] = []
    targets: list[int] = []
    equations: list[Polynomial] = []
    for left_index, right_index in itertools.combinations(
        range(len(remainder)), 2
    ):
        left = masks[left_index]
        right = masks[right_index]
        target = int(bool(
            adj[remainder[left_index]] & (1 << remainder[right_index])
        ))
        polynomial: Polynomial = {}

        # T*w(A intersect B).
        for coordinate in range(u):
            if not (left[coordinate] and right[coordinate]):
                continue
            linear = list(zero)
            linear[coordinate] = 1
            add_term(polynomial, tuple(linear), Q(1))
            for other in range(u):
                quadratic = linear.copy()
                quadratic[other] += 1
                add_term(polynomial, tuple(quadratic), Q(1))

        # -w(A)*w(B).
        for first in range(u):
            if not left[first]:
                continue
            for second in range(u):
                if not right[second]:
                    continue
                quadratic = list(zero)
                quadratic[first] += 1
                quadratic[second] += 1
                add_term(polynomial, tuple(quadratic), Q(-1))

        # -k_AB*T.
        if target:
            add_term(polynomial, zero, Q(-1))
            for coordinate in range(u):
                linear = list(zero)
                linear[coordinate] = 1
                add_term(polynomial, tuple(linear), Q(-1))

        pairs.append((remainder[left_index], remainder[right_index]))
        targets.append(target)
        equations.append(polynomial)
    return SchurPolynomialSystem(
        u=u,
        clique=clique,
        remainder=remainder,
        masks=masks,
        equation_pairs=tuple(pairs),
        equation_targets=tuple(targets),
        equations=tuple(equations),
    )


def coefficient_matrix(
    equations: Sequence[Polynomial],
    variables: int,
    multiplier_degree: int,
) -> tuple[
    tuple[Exponent, ...],
    tuple[Exponent, ...],
    tuple[tuple[int, ...], ...],
]:
    """Integer matrix from equation multipliers to output coefficients."""

    multipliers = monomials(variables, multiplier_degree)
    outputs = monomials(variables, multiplier_degree + 2)
    output_index = {exponent: index for index, exponent in enumerate(outputs)}
    columns = len(equations) * len(multipliers)
    matrix = [[0] * columns for _ in outputs]
    for equation_index, equation in enumerate(equations):
        for multiplier_index, multiplier in enumerate(multipliers):
            column = equation_index * len(multipliers) + multiplier_index
            for exponent, coefficient in equation.items():
                output = tuple(
                    left + right
                    for left, right in zip(exponent, multiplier)
                )
                matrix[output_index[output]][column] += int(coefficient)
    return multipliers, outputs, tuple(tuple(row) for row in matrix)


def _guided_rational_point(
    system: strict_h.ExactAffineSystem,
    approximate: Sequence[float],
    max_denominator: int,
) -> tuple[Q, ...]:
    """Rational point in an exact affine space, close to ``approximate``."""

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
        point[pivot] = rhs - sum(row[column] * point[column] for column in free)
    return tuple(point)


def verify_raw_certificate(
    equations: Sequence[Polynomial],
    variables: int,
    certificate: dict,
) -> Polynomial:
    """Verify a serialized rational identity and return its positive side."""

    combined: Polynomial = {}
    for item in certificate.get("multipliers", []):
        equation = int(item["equation"])
        if not 0 <= equation < len(equations):
            raise ValueError("certificate equation index is out of range")
        monomial = tuple(int(value) for value in item["monomial"])
        if len(monomial) != variables or any(value < 0 for value in monomial):
            raise ValueError("invalid multiplier monomial")
        coefficient = parse_fraction(item["coefficient"])
        for exponent, value in shifted(equations[equation], monomial).items():
            add_term(combined, exponent, coefficient * value)

    claimed: Polynomial = {}
    for item in certificate.get("positive_polynomial", []):
        exponent = tuple(int(value) for value in item["monomial"])
        if len(exponent) != variables or any(value < 0 for value in exponent):
            raise ValueError("invalid positive-polynomial monomial")
        coefficient = parse_fraction(item["coefficient"])
        if coefficient <= 0:
            raise ValueError("positive-polynomial coefficient is not positive")
        add_term(claimed, exponent, coefficient)
    if not claimed:
        raise ValueError("certificate polynomial is zero")
    if combined != claimed:
        raise ValueError("certificate polynomial identity does not check")
    return claimed


def verify_clique_certificate(
    adj: Sequence[int], clique_mask: int, certificate: dict
) -> Polynomial:
    system = schur_polynomial_system(adj, clique_mask)
    if list(system.clique) != certificate.get("clique"):
        raise ValueError("certificate clique does not match")
    if [list(pair) for pair in system.equation_pairs] != certificate.get(
        "equation_pairs"
    ):
        raise ValueError("certificate equation order does not match")
    return verify_raw_certificate(system.equations, system.u, certificate)


def find_raw_certificate(
    equations: Sequence[Polynomial],
    variables: int,
    multiplier_degree: int = 1,
) -> dict | None:
    """Locate, reconstruct, and exactly verify a polynomial certificate."""

    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        return None

    multipliers, outputs, matrix = coefficient_matrix(
        equations, variables, multiplier_degree
    )
    numeric = np.asarray(matrix, dtype=float)
    # Find q with Cq >= 0 and sum(Cq)=1.  The q multipliers are free.
    result = linprog(
        np.zeros(numeric.shape[1]),
        A_ub=-numeric,
        b_ub=np.zeros(numeric.shape[0]),
        A_eq=[numeric.sum(axis=0)],
        b_eq=[1.0],
        bounds=[(None, None)] * numeric.shape[1],
        method="highs-ds",
    )
    if not result.success or result.x is None:
        return None
    approximate_polynomial = numeric @ result.x
    normalization = [sum(row[column] for row in matrix) for column in range(len(matrix[0]))]

    # A basic HiGHS solution exposes the zero output coefficients.  Enforce
    # those zeros and the normalization exactly, then choose nearby rational
    # free coordinates.  Several cutoffs make reconstruction robust to LP
    # feasibility tolerances.
    for cutoff in (1e-7, 1e-9, 1e-11, 1e-13):
        zero_rows = [
            row for row, value in zip(matrix, approximate_polynomial)
            if value <= cutoff
        ]
        affine = strict_h.exact_rref(
            zero_rows + [normalization],
            [0] * len(zero_rows) + [1],
            numeric.shape[1],
        )
        if not affine.consistent:
            continue
        for maximum_denominator in (1_000, 1_000_000, 1_000_000_000):
            rational = _guided_rational_point(
                affine, result.x, maximum_denominator
            )
            positive = [
                sum(Q(coefficient) * value for coefficient, value in zip(row, rational))
                for row in matrix
            ]
            if any(value < 0 for value in positive) or sum(positive) != 1:
                continue
            certificate = {
                "schema": 1,
                "kind": "positive_orthant_polynomial_farkas",
                "multiplier_degree": multiplier_degree,
                "multipliers": [
                    {
                        "equation": column // len(multipliers),
                        "monomial": list(multipliers[column % len(multipliers)]),
                        "coefficient": fraction_json(value),
                    }
                    for column, value in enumerate(rational)
                    if value
                ],
                "positive_polynomial": [
                    {
                        "monomial": list(outputs[row]),
                        "coefficient": fraction_json(value),
                    }
                    for row, value in enumerate(positive)
                    if value
                ],
            }
            verify_raw_certificate(equations, variables, certificate)
            return certificate
    return None


def find_clique_certificate(
    adj: Sequence[int], clique_mask: int, multiplier_degree: int = 1
) -> dict | None:
    system = schur_polynomial_system(adj, clique_mask)
    certificate = find_raw_certificate(
        system.equations, system.u, multiplier_degree
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


@dataclass
class GraphDualResult:
    index: int | None
    applicable: bool
    seeds: int
    covers: int
    enhanced_passing_covers: int
    strict_h_passing_covers: int
    dual_passing_covers: int
    strict_h_passing_cliques: int
    dual_failing_cliques: int
    strict_h_rejected: bool
    dual_rejected: bool
    marginal_dual_rejected: bool
    first_certificate: dict | None = None
    dual_failure_witnesses: list[dict] = field(default_factory=list)
    elapsed_seconds: float = 0.0


def analyze_graph(graph: dict, multiplier_degree: int = 1) -> GraphDualResult:
    """Apply the exact dual with the same seed/cover quantifiers as prior."""

    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    seeds = covers = enhanced_passes = strict_passes = dual_passes = 0
    strict_cliques = dual_failures = 0
    strict_rejected = dual_rejected = False
    first_certificate = None
    dual_failure_witnesses: list[dict] = []
    for seed_mask in prior.clique_masks(adj, 7):
        seeds += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(adj, seed_mask)
        total_term_rank = prior.matching_size(defects)
        seed_strict_passes = seed_dual_passes = 0
        seed_covers = prior.eligible_covers(ladj, eligible)
        covers += len(seed_covers)
        for zmask in seed_covers:
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            enhanced_passes += 1
            nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            clique_masks = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            strict_results = [
                strict_h.assess_saturating_clique(graph_n, mask)
                for mask in clique_masks
            ]
            if any(result.failed for result in strict_results):
                continue
            seed_strict_passes += 1
            strict_passes += 1
            strict_cliques += len(clique_masks)
            failure = None
            for mask in clique_masks:
                certificate = find_clique_certificate(
                    graph_n, mask, multiplier_degree
                )
                if certificate is not None:
                    failure = certificate
                    dual_failures += 1
                    break
            if failure is not None:
                cover_witness = {
                    "seed": seed,
                    "zmask": zmask,
                    "nvertices": nvertices,
                    "certificate": failure,
                }
                dual_failure_witnesses.append(cover_witness)
                if first_certificate is None:
                    first_certificate = cover_witness
                continue
            seed_dual_passes += 1
            dual_passes += 1
        if not seed_strict_passes:
            strict_rejected = True
        if not seed_dual_passes:
            dual_rejected = True
    return GraphDualResult(
        index=graph.get("index"),
        applicable=bool(seeds),
        seeds=seeds,
        covers=covers,
        enhanced_passing_covers=enhanced_passes,
        strict_h_passing_covers=strict_passes,
        dual_passing_covers=dual_passes,
        strict_h_passing_cliques=strict_cliques,
        dual_failing_cliques=dual_failures,
        strict_h_rejected=strict_rejected,
        dual_rejected=dual_rejected,
        marginal_dual_rejected=dual_rejected and not strict_rejected,
        first_certificate=first_certificate,
        dual_failure_witnesses=dual_failure_witnesses,
        elapsed_seconds=time.perf_counter() - started,
    )


def sample_report(
    sample: dict, selected_indices: set[int], multiplier_degree: int
) -> dict:
    graphs = [
        graph for graph in sample["graphs"]
        if graph.get("index") in selected_indices
    ]
    started = time.perf_counter()
    results = [analyze_graph(graph, multiplier_degree) for graph in graphs]
    indices = lambda predicate: sorted(
        result.index for result in results
        if result.index is not None and predicate(result)
    )
    return {
        "schema": 1,
        "description": (
            "Exact degree-bounded positive-polynomial Farkas dual after "
            "the strict-affine-H K7 screen. Numerical LP only locates; every "
            "reported rejection includes a rational identity."
        ),
        "multiplier_degree": multiplier_degree,
        "graphs": len(results),
        "counts": {
            "seeds": sum(result.seeds for result in results),
            "covers": sum(result.covers for result in results),
            "enhanced_passing_covers": sum(
                result.enhanced_passing_covers for result in results
            ),
            "strict_H_passing_covers": sum(
                result.strict_h_passing_covers for result in results
            ),
            "dual_passing_covers": sum(
                result.dual_passing_covers for result in results
            ),
            "strict_H_passing_cliques": sum(
                result.strict_h_passing_cliques for result in results
            ),
            "dual_failing_cliques_first_per_cover": sum(
                result.dual_failing_cliques for result in results
            ),
            "marginal_dual_rejected": sum(
                result.marginal_dual_rejected for result in results
            ),
            "dual_survivors": sum(
                result.applicable and not result.dual_rejected
                for result in results
            ),
        },
        "decisions": {
            "marginal_dual_rejected": indices(
                lambda result: result.marginal_dual_rejected
            ),
            "dual_survivors": indices(
                lambda result: result.applicable and not result.dual_rejected
            ),
        },
        "runtime": {
            "wall_seconds": time.perf_counter() - started,
            "sum_graph_seconds": sum(result.elapsed_seconds for result in results),
        },
        "per_graph": [asdict(result) for result in results],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("--selection-report", type=Path, required=True)
    parser.add_argument("--multiplier-degree", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    sample = json.loads(args.sample.read_text(encoding="utf-8"))
    selection = json.loads(args.selection_report.read_text(encoding="utf-8"))
    selected = set(selection["graph_decisions"]["combined_survivors"])
    report = sample_report(sample, selected, args.multiplier_degree)
    report["input"] = {
        "sample": {"file": args.sample.name, "sha256": file_sha256(args.sample)},
        "selection_report": {
            "file": args.selection_report.name,
            "sha256": file_sha256(args.selection_report),
        },
        "selected_indices": sorted(selected),
    }
    source = Path(__file__)
    report["source"] = {"file": source.name, "sha256": file_sha256(source)}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
