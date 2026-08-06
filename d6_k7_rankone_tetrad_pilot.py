#!/usr/bin/env python3
"""Tiny exact rank-one Schur tetrad pilot for K7 residue covers.

For a required clique one smaller than the K-rank upper bound, the Schur
complement is positive semidefinite of rank at most one.  Its known
off-diagonal polynomial numerators therefore satisfy tetrad equalities and
triangle sign inequalities.  HiGHS only locates rational positive-polynomial
identities; every accepted identity is re-expanded exactly with ``Fraction``.
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

import d6_k7_positive_polynomial_dual as saturating_dual
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


def subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    output = dict(left)
    for exponent, value in right.items():
        add_term(output, exponent, -value)
    return output


@dataclass(frozen=True)
class RankOneSystem:
    clique: tuple[int, ...]
    remainder: tuple[int, ...]
    masks: tuple[tuple[int, ...], ...]
    pair_labels: tuple[tuple[int, int], ...]
    pair_polynomials: tuple[Polynomial, ...]
    tetrad_labels: tuple[tuple[tuple[int, int], ...], ...]
    tetrad_equations: tuple[Polynomial, ...]
    triangle_labels: tuple[tuple[int, int, int], ...]
    triangle_signs: tuple[Polynomial, ...]


def rank_one_system(adj: Sequence[int], clique_mask: int) -> RankOneSystem:
    """Build exact off-diagonal rank-one Schur consequences."""

    schur = saturating_dual.schur_polynomial_system(adj, clique_mask)
    # The saturating equation is f_yz=-g_yz.  Here g_yz=T*S_yz.
    pair_polynomials = tuple(
        {exponent: -value for exponent, value in equation.items()}
        for equation in schur.equations
    )
    pair_by_label = dict(zip(schur.equation_pairs, pair_polynomials))
    tetrad_labels = []
    tetrads = []
    for a, b, c, d in itertools.combinations(schur.remainder, 4):
        ab_cd = multiply(pair_by_label[(a, b)], pair_by_label[(c, d)])
        ac_bd = multiply(pair_by_label[(a, c)], pair_by_label[(b, d)])
        ad_bc = multiply(pair_by_label[(a, d)], pair_by_label[(b, c)])
        first = subtract(ab_cd, ac_bd)
        second = subtract(ab_cd, ad_bc)
        if first:
            tetrad_labels.append(((a, b), (c, d), (a, c), (b, d)))
            tetrads.append(first)
        if second:
            tetrad_labels.append(((a, b), (c, d), (a, d), (b, c)))
            tetrads.append(second)
    triangle_labels = []
    triangle_signs = []
    for a, b, c in itertools.combinations(schur.remainder, 3):
        sign = multiply(
            multiply(pair_by_label[(a, b)], pair_by_label[(a, c)]),
            pair_by_label[(b, c)],
        )
        if sign:
            triangle_labels.append((a, b, c))
            triangle_signs.append(sign)
    return RankOneSystem(
        clique=schur.clique,
        remainder=schur.remainder,
        masks=schur.masks,
        pair_labels=schur.equation_pairs,
        pair_polynomials=pair_polynomials,
        tetrad_labels=tuple(tetrad_labels),
        tetrad_equations=tuple(tetrads),
        triangle_labels=tuple(triangle_labels),
        triangle_signs=tuple(triangle_signs),
    )


def polynomial_degree(polynomial: Polynomial) -> int:
    return max(map(sum, polynomial), default=-1)


def verify_raw_certificate(
    equations: Sequence[Polynomial],
    signs: Sequence[Polynomial],
    variables: int,
    certificate: dict,
) -> Polynomial:
    """Check ``sum qE = P + sum r*t`` exactly, with P strictly positive."""

    left: Polynomial = {}
    for item in certificate.get("equation_multipliers", []):
        equation_index = int(item["equation"])
        if not 0 <= equation_index < len(equations):
            raise ValueError("tetrad equation index out of range")
        monomial = tuple(int(value) for value in item["monomial"])
        if len(monomial) != variables or any(value < 0 for value in monomial):
            raise ValueError("invalid tetrad multiplier monomial")
        coefficient = parse_fraction(item["coefficient"])
        for exponent, value in shifted(
            equations[equation_index], monomial
        ).items():
            add_term(left, exponent, coefficient * value)

    right: Polynomial = {}
    for item in certificate.get("sign_multipliers", []):
        sign_index = int(item["sign"])
        if not 0 <= sign_index < len(signs):
            raise ValueError("triangle sign index out of range")
        monomial = tuple(int(value) for value in item["monomial"])
        if len(monomial) != variables or any(value < 0 for value in monomial):
            raise ValueError("invalid sign multiplier monomial")
        coefficient = parse_fraction(item["coefficient"])
        if coefficient <= 0:
            raise ValueError("triangle sign multiplier is not positive")
        for exponent, value in shifted(signs[sign_index], monomial).items():
            add_term(right, exponent, coefficient * value)

    positive: Polynomial = {}
    for item in certificate.get("positive_polynomial", []):
        exponent = tuple(int(value) for value in item["monomial"])
        if len(exponent) != variables or any(value < 0 for value in exponent):
            raise ValueError("invalid positive-polynomial monomial")
        coefficient = parse_fraction(item["coefficient"])
        if coefficient <= 0:
            raise ValueError("positive-polynomial coefficient is not positive")
        add_term(positive, exponent, coefficient)
        add_term(right, exponent, coefficient)
    if not positive:
        raise ValueError("rank-one certificate needs a strict positive part")
    if left != right:
        raise ValueError("rank-one polynomial identity does not check")
    return positive


def find_raw_certificate(
    equations: Sequence[Polynomial],
    signs: Sequence[Polynomial],
    variables: int,
    output_degree: int,
) -> dict | None:
    """LP-locate and exactly check a degree-truncated rank-one identity."""

    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        return None
    equation_monomials = tuple(
        monomials(variables, output_degree - polynomial_degree(equation))
        for equation in equations
    )
    sign_monomials = tuple(
        monomials(variables, output_degree - polynomial_degree(sign))
        for sign in signs
    )
    outputs = monomials(variables, output_degree)
    output_index = {exponent: index for index, exponent in enumerate(outputs)}
    equation_columns = sum(map(len, equation_monomials))
    sign_columns = sum(map(len, sign_monomials))
    columns = equation_columns + sign_columns
    if not columns:
        return None
    matrix = [[0] * columns for _ in outputs]
    column = 0
    for equation, multipliers in zip(equations, equation_monomials):
        for monomial in multipliers:
            for exponent, coefficient in shifted(equation, monomial).items():
                matrix[output_index[exponent]][column] += int(coefficient)
            column += 1
    for sign, multipliers in zip(signs, sign_monomials):
        for monomial in multipliers:
            # sum qE = P + sum r*t, so output is qE-r*t.
            for exponent, coefficient in shifted(sign, monomial).items():
                matrix[output_index[exponent]][column] -= int(coefficient)
            column += 1
    numeric = np.asarray(matrix, dtype=float)
    normalization = numeric.sum(axis=0)
    result = linprog(
        np.zeros(columns),
        A_ub=-numeric,
        b_ub=np.zeros(len(outputs)),
        A_eq=[normalization],
        b_eq=[1.0],
        bounds=(
            [(None, None)] * equation_columns
            + [(0.0, None)] * sign_columns
        ),
        method="highs-ds",
    )
    if not result.success or result.x is None:
        return None

    # No equality is part of the mathematical certificate: normalization is
    # only an LP gauge.  Rationalize the multipliers directly, recompute the
    # output exactly, and accept only if every exact coefficient is >=0.
    for maximum_denominator in (1_000, 1_000_000, 1_000_000_000):
        rational = tuple(
            Q(float(value)).limit_denominator(maximum_denominator)
            for value in result.x
        )
        sign_values = rational[equation_columns:]
        if any(value < 0 for value in sign_values):
            continue
        positive = [
            sum(Q(coefficient) * value for coefficient, value in zip(row, rational))
            for row in matrix
        ]
        if any(value < 0 for value in positive) or not any(positive):
            continue
        equation_items = []
        column = 0
        for equation_index, multipliers in enumerate(equation_monomials):
            for monomial in multipliers:
                value = rational[column]
                if value:
                    equation_items.append({
                        "equation": equation_index,
                        "monomial": list(monomial),
                        "coefficient": fraction_json(value),
                    })
                column += 1
        sign_items = []
        for sign_index, multipliers in enumerate(sign_monomials):
            for monomial in multipliers:
                value = rational[column]
                if value:
                    sign_items.append({
                        "sign": sign_index,
                        "monomial": list(monomial),
                        "coefficient": fraction_json(value),
                    })
                column += 1
        certificate = {
            "schema": 1,
            "kind": "rank_one_schur_tetrad_positive_polynomial",
            "output_degree": output_degree,
            "equation_multipliers": equation_items,
            "sign_multipliers": sign_items,
            "positive_polynomial": [
                {
                    "monomial": list(outputs[row]),
                    "coefficient": fraction_json(value),
                }
                for row, value in enumerate(positive) if value
            ],
        }
        verify_raw_certificate(equations, signs, variables, certificate)
        return certificate
    return None


def find_clique_certificate(
    adj: Sequence[int], clique_mask: int, output_degree: int = 6
) -> dict | None:
    system = rank_one_system(adj, clique_mask)
    certificate = find_raw_certificate(
        system.tetrad_equations,
        system.triangle_signs if output_degree >= 6 else (),
        len(system.clique),
        output_degree,
    )
    if certificate is None:
        return None
    certificate.update({
        "clique": list(system.clique),
        "remainder": list(system.remainder),
        "basis_neighbour_masks": [list(mask) for mask in system.masks],
        "pair_labels": [list(pair) for pair in system.pair_labels],
        "tetrad_labels": [
            [list(pair) for pair in label] for label in system.tetrad_labels
        ],
        "triangle_labels": [list(label) for label in system.triangle_labels],
    })
    verify_clique_certificate(adj, clique_mask, certificate)
    return certificate


def verify_clique_certificate(
    adj: Sequence[int], clique_mask: int, certificate: dict
) -> None:
    system = rank_one_system(adj, clique_mask)
    if certificate.get("clique") != list(system.clique):
        raise ValueError("certificate clique mismatch")
    if certificate.get("remainder") != list(system.remainder):
        raise ValueError("certificate remainder mismatch")
    if certificate.get("basis_neighbour_masks") != [
        list(mask) for mask in system.masks
    ]:
        raise ValueError("certificate masks mismatch")
    if certificate.get("pair_labels") != [
        list(pair) for pair in system.pair_labels
    ]:
        raise ValueError("certificate pair ordering mismatch")
    if certificate.get("tetrad_labels") != [
        [list(pair) for pair in label] for label in system.tetrad_labels
    ]:
        raise ValueError("certificate tetrad ordering mismatch")
    if certificate.get("triangle_labels") != [
        list(label) for label in system.triangle_labels
    ]:
        raise ValueError("certificate triangle ordering mismatch")
    signs = system.triangle_signs if certificate["output_degree"] >= 6 else ()
    verify_raw_certificate(
        system.tetrad_equations, signs, len(system.clique), certificate
    )


@dataclass
class GraphResult:
    index: int
    seeds: int = 0
    covers: int = 0
    prior_passing_covers: int = 0
    no_saturating_clique_covers: int = 0
    covers_with_near_clique: int = 0
    covers_without_near_clique: int = 0
    near_cliques_tested: int = 0
    near_cliques_tested_degree4: int = 0
    near_cliques_tested_degree6: int = 0
    degree4_failed_covers: int = 0
    degree6_sign_new_failed_covers: int = 0
    tetrad_failed_covers: int = 0
    final_passing_covers: int = 0
    marginal_graph_rejected: bool = False
    certificates: list[dict] = field(default_factory=list)
    elapsed_seconds: float = 0.0


def analyze_graph(graph: dict, output_degree: int = 6) -> GraphResult:
    """Apply rank-one tetrads only to post-degree-one no-K_U covers."""

    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    result = GraphResult(index=int(graph["index"]))
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    marginal_rejected = False
    for seed_mask in prior.clique_masks(adj, 7):
        result.seeds += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(
            adj, seed_mask
        )
        total_term_rank = prior.matching_size(defects)
        seed_passes = 0
        for zmask in prior.eligible_covers(ladj, eligible):
            result.covers += 1
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            saturating = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            if any(
                strict_h.assess_saturating_clique(graph_n, mask).failed
                for mask in saturating
            ):
                continue
            baseline_failure = next((
                certificate
                for mask in saturating
                if (certificate := saturating_dual.find_clique_certificate(
                    graph_n, mask, multiplier_degree=1
                )) is not None
            ), None)
            if baseline_failure is not None:
                continue
            result.prior_passing_covers += 1
            if saturating:
                seed_passes += 1
                result.final_passing_covers += 1
                continue
            result.no_saturating_clique_covers += 1
            near_size = analysis.k_rank_upper - 1
            near_cliques = tuple(prior.clique_masks(
                graph_n, near_size
            )) if near_size > 0 else ()
            if not near_cliques:
                result.covers_without_near_clique += 1
                seed_passes += 1
                result.final_passing_covers += 1
                continue
            result.covers_with_near_clique += 1
            failure = None
            for near_mask in near_cliques:
                result.near_cliques_tested += 1
                result.near_cliques_tested_degree4 += 1
                certificate = find_clique_certificate(graph_n, near_mask, 4)
                if certificate is not None:
                    failure = {
                        "stage": "degree4_tetrads",
                        "seed": seed,
                        "zmask": zmask,
                        "nvertices": nvertices,
                        "k_rank_upper": analysis.k_rank_upper,
                        "certificate": certificate,
                    }
                    break
            if failure is None and output_degree >= 6:
                for near_mask in near_cliques:
                    result.near_cliques_tested += 1
                    result.near_cliques_tested_degree6 += 1
                    certificate = find_clique_certificate(
                        graph_n, near_mask, output_degree
                    )
                    if certificate is not None:
                        failure = {
                            "stage": "degree6_tetrads_and_signs",
                            "seed": seed,
                            "zmask": zmask,
                            "nvertices": nvertices,
                            "k_rank_upper": analysis.k_rank_upper,
                            "certificate": certificate,
                        }
                        break
            if failure is not None:
                result.tetrad_failed_covers += 1
                if failure["stage"] == "degree4_tetrads":
                    result.degree4_failed_covers += 1
                else:
                    result.degree6_sign_new_failed_covers += 1
                result.certificates.append(failure)
                continue
            seed_passes += 1
            result.final_passing_covers += 1
        if not seed_passes:
            marginal_rejected = True
    result.marginal_graph_rejected = marginal_rejected
    result.elapsed_seconds = time.perf_counter() - started
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample", type=Path,
        default=Path("d6_k7_rankone_tetrad_sample.json"),
    )
    parser.add_argument("--indices", type=int, nargs="+", required=True)
    parser.add_argument("--output-degree", type=int, default=6)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sample = json.loads(args.sample.read_text(encoding="utf-8"))
    requested = set(args.indices)
    graphs = [
        graph for graph in sample["graphs"] if int(graph["index"]) in requested
    ]
    if {int(graph["index"]) for graph in graphs} != requested:
        raise ValueError("one or more requested graph indices are absent")
    started = time.perf_counter()
    results = [
        analyze_graph(graph, args.output_degree)
        for graph in sorted(graphs, key=lambda item: int(item["index"]))
    ]
    report = {
        "schema": 1,
        "kind": "d6_k7_rankone_tetrad_pilot",
        "input": {
            "path": args.sample.name,
            "sha256": file_sha256(args.sample),
            "indices": sorted(requested),
        },
        "parameters": {"output_degree": args.output_degree},
        "counts": {
            "graphs": len(results),
            "prior_passing_covers": sum(
                item.prior_passing_covers for item in results
            ),
            "no_saturating_clique_covers": sum(
                item.no_saturating_clique_covers for item in results
            ),
            "covers_with_near_clique": sum(
                item.covers_with_near_clique for item in results
            ),
            "covers_without_near_clique": sum(
                item.covers_without_near_clique for item in results
            ),
            "near_cliques_tested": sum(
                item.near_cliques_tested for item in results
            ),
            "near_cliques_tested_degree4": sum(
                item.near_cliques_tested_degree4 for item in results
            ),
            "near_cliques_tested_degree6": sum(
                item.near_cliques_tested_degree6 for item in results
            ),
            "degree4_failed_covers": sum(
                item.degree4_failed_covers for item in results
            ),
            "degree6_sign_new_failed_covers": sum(
                item.degree6_sign_new_failed_covers for item in results
            ),
            "tetrad_failed_covers": sum(
                item.tetrad_failed_covers for item in results
            ),
            "marginal_graph_rejections": sum(
                item.marginal_graph_rejected for item in results
            ),
            "final_passing_covers": sum(
                item.final_passing_covers for item in results
            ),
        },
        "runtime": {"wall_seconds": time.perf_counter() - started},
        "per_graph": [asdict(item) for item in results],
    }
    source = Path(__file__)
    report["source"] = {
        "path": source.name, "sha256": file_sha256(source)
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
