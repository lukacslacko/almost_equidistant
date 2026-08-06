#!/usr/bin/env python3
"""Independent exact verifier for the K7 degree-two/preordering pilot.

This checker does not import the numerical locator.  It reconstructs the
Schur equations, strict outside-diagonal generators, and every retained
rational polynomial identity using arbitrary-precision ``Fraction``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as prior
import d6_k7_special_h_reference as strict_h
import verify_d6_k7_positive_polynomial_dual as verify_degree_one


Q = Fraction
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Q]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_fraction(value: object) -> Q:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or not all(isinstance(item, int) for item in value)
    ):
        raise ValueError("invalid serialized fraction")
    return Q(value[0], value[1])


def add_term(polynomial: Polynomial, exponent: Exponent, value: Q) -> None:
    if not value:
        return
    updated = polynomial.get(exponent, Q(0)) + value
    if updated:
        polynomial[exponent] = updated
    else:
        polynomial.pop(exponent, None)


def shift(polynomial: Polynomial, monomial: Exponent) -> Polynomial:
    return {
        tuple(a + b for a, b in zip(exponent, monomial)): coefficient
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


def bits(mask: int) -> tuple[int, ...]:
    output = []
    while mask:
        bit = mask & -mask
        output.append(bit.bit_length() - 1)
        mask ^= bit
    return tuple(output)


def reconstruct_system(
    adj: Sequence[int], clique_mask: int
) -> tuple[
    tuple[int, ...],
    tuple[int, ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, int], ...],
    tuple[int, ...],
    tuple[Polynomial, ...],
]:
    clique = bits(clique_mask)
    if not clique:
        raise ValueError("empty clique")
    if any(
        not (adj[first] & (1 << second))
        for first, second in itertools.combinations(clique, 2)
    ):
        raise ValueError("certificate basis is not a required clique")
    remainder = tuple(
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    )
    masks = tuple(
        tuple(int(bool(adj[vertex] & (1 << basis))) for basis in clique)
        for vertex in remainder
    )
    variables = len(clique)
    zero = (0,) * variables
    pairs = []
    targets = []
    equations = []
    for left_index, right_index in itertools.combinations(
        range(len(remainder)), 2
    ):
        left = masks[left_index]
        right = masks[right_index]
        target = int(bool(
            adj[remainder[left_index]] & (1 << remainder[right_index])
        ))
        polynomial: Polynomial = {}
        for coordinate in range(variables):
            if not (left[coordinate] and right[coordinate]):
                continue
            linear = list(zero)
            linear[coordinate] = 1
            add_term(polynomial, tuple(linear), Q(1))
            for other in range(variables):
                quadratic = linear.copy()
                quadratic[other] += 1
                add_term(polynomial, tuple(quadratic), Q(1))
        for first in range(variables):
            if not left[first]:
                continue
            for second in range(variables):
                if not right[second]:
                    continue
                quadratic = [0] * variables
                quadratic[first] += 1
                quadratic[second] += 1
                add_term(polynomial, tuple(quadratic), Q(-1))
        if target:
            add_term(polynomial, zero, Q(-1))
            for coordinate in range(variables):
                linear = list(zero)
                linear[coordinate] = 1
                add_term(polynomial, tuple(linear), Q(-1))
        pairs.append((remainder[left_index], remainder[right_index]))
        targets.append(target)
        equations.append(polynomial)
    return (
        clique,
        remainder,
        masks,
        tuple(pairs),
        tuple(targets),
        tuple(equations),
    )


def outside_diagonal(mask: Sequence[int]) -> Polynomial:
    variables = len(mask)
    zero = (0,) * variables
    polynomial: Polynomial = {zero: Q(-1)}
    for coordinate in range(variables):
        if not mask[coordinate]:
            continue
        linear = list(zero)
        linear[coordinate] = 1
        add_term(polynomial, tuple(linear), Q(1))
        for other in range(variables):
            quadratic = linear.copy()
            quadratic[other] += 1
            add_term(polynomial, tuple(quadratic), Q(1))
    for first in range(variables):
        if not mask[first]:
            continue
        for second in range(variables):
            if not mask[second]:
                continue
            quadratic = [0] * variables
            quadratic[first] += 1
            quadratic[second] += 1
            add_term(polynomial, tuple(quadratic), Q(-1))
    for coordinate in range(variables):
        linear = list(zero)
        linear[coordinate] = 1
        add_term(polynomial, tuple(linear), Q(-1))
    return polynomial


def reconstruct_generators(
    remainder: Sequence[int],
    masks: Sequence[Sequence[int]],
    maximum_order: int,
) -> dict[tuple[int, ...], Polynomial]:
    base = {
        (vertex,): outside_diagonal(mask)
        for vertex, mask in zip(remainder, masks)
    }
    output = dict(base)
    for order in range(2, maximum_order + 1):
        for vertices in itertools.combinations(remainder, order):
            polynomial: Polynomial = {(0,) * len(masks[0]): Q(1)}
            for vertex in vertices:
                polynomial = multiply(polynomial, base[(vertex,)])
            output[tuple(vertices)] = polynomial
    return output


def verify_clique_certificate(
    adj: Sequence[int], clique_mask: int, certificate: dict
) -> None:
    (
        clique,
        remainder,
        masks,
        pairs,
        targets,
        equations,
    ) = reconstruct_system(adj, clique_mask)
    if certificate.get("schema") != 1:
        raise ValueError("unexpected certificate schema")
    if certificate.get("kind") != (
        "positive_orthant_outside_diagonal_preordering"
    ):
        raise ValueError("unexpected certificate kind")
    if certificate.get("clique") != list(clique):
        raise ValueError("certificate clique mismatch")
    if certificate.get("remainder") != list(remainder):
        raise ValueError("certificate remainder mismatch")
    if certificate.get("basis_neighbour_masks") != [list(mask) for mask in masks]:
        raise ValueError("certificate basis-mask mismatch")
    if certificate.get("equation_pairs") != [list(pair) for pair in pairs]:
        raise ValueError("certificate equation-pair mismatch")
    if certificate.get("equation_targets") != list(targets):
        raise ValueError("certificate equation-target mismatch")
    output_degree = certificate.get("output_degree")
    maximum_order = certificate.get("maximum_generator_order")
    if not isinstance(output_degree, int) or output_degree < 2:
        raise ValueError("invalid output degree")
    if not isinstance(maximum_order, int) or maximum_order < 1:
        raise ValueError("invalid generator order")
    variables = len(clique)
    generators = reconstruct_generators(remainder, masks, maximum_order)

    left: Polynomial = {}
    for item in certificate.get("equation_multipliers", []):
        equation_index = item.get("equation")
        if not isinstance(equation_index, int) or not 0 <= equation_index < len(
            equations
        ):
            raise ValueError("equation index out of range")
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid equation multiplier")
        monomial = tuple(monomial_value)
        if (
            len(monomial) != variables
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) > output_degree - 2
        ):
            raise ValueError("equation multiplier exceeds degree bound")
        coefficient = parse_fraction(item.get("coefficient"))
        for exponent, value in shift(equations[equation_index], monomial).items():
            add_term(left, exponent, coefficient * value)

    right: Polynomial = {}
    positive_mass = Q(0)
    for item in certificate.get("generator_multipliers", []):
        vertices_value = item.get("vertices")
        if not isinstance(vertices_value, list):
            raise ValueError("invalid generator vertex list")
        vertices = tuple(vertices_value)
        generator = generators.get(vertices)
        if generator is None:
            raise ValueError("unknown generator")
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid generator multiplier")
        monomial = tuple(monomial_value)
        generator_degree = max(map(sum, generator))
        if (
            len(monomial) != variables
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) > output_degree - generator_degree
        ):
            raise ValueError("generator multiplier exceeds degree bound")
        coefficient = parse_fraction(item.get("coefficient"))
        if coefficient <= 0:
            raise ValueError("generator multiplier is not positive")
        positive_mass += coefficient
        for exponent, value in shift(generator, monomial).items():
            add_term(right, exponent, coefficient * value)

    for item in certificate.get("positive_polynomial", []):
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid positive polynomial")
        monomial = tuple(monomial_value)
        if (
            len(monomial) != variables
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) > output_degree
        ):
            raise ValueError("positive monomial exceeds degree bound")
        coefficient = parse_fraction(item.get("coefficient"))
        if coefficient <= 0:
            raise ValueError("positive polynomial has nonpositive coefficient")
        positive_mass += coefficient
        add_term(right, monomial, coefficient)
    if positive_mass <= 0:
        raise ValueError("zero certificate")
    if left != right:
        raise ValueError("exact preordering identity does not check")


def verify_report(sample_path: Path, report_path: Path) -> dict:
    sample = json.loads(sample_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("schema") != 1 or report.get("kind") != (
        "d6_k7_dual_degree2_pilot"
    ):
        raise ValueError("unexpected report schema")
    if report["input"]["sha256"] != file_sha256(sample_path):
        raise ValueError("sample hash mismatch")
    source_path = report_path.parent / report["source"]["path"]
    if report["source"]["sha256"] != file_sha256(source_path):
        raise ValueError("pilot source hash mismatch")
    graph_by_index = {
        int(graph["index"]): graph for graph in sample["graphs"]
    }
    recorded_indices = [int(item["index"]) for item in report["per_graph"]]
    if recorded_indices != report["input"]["indices"]:
        raise ValueError("report graph order differs from selected indices")
    if len(recorded_indices) != len(set(recorded_indices)):
        raise ValueError("report repeats a graph")
    checked = 0
    for recorded in report["per_graph"]:
        graph = graph_by_index[int(recorded["index"])]
        adj = tuple(graph["adjacency"])
        prior.validate_graph(adj)
        for witness in recorded["certificates"]:
            seed = tuple(int(vertex) for vertex in witness["seed"])
            seed_mask = sum(1 << vertex for vertex in seed)
            if seed_mask not in set(prior.clique_masks(adj, 7)):
                raise ValueError("witness seed is not a required K7")
            rebuilt_seed, outside, defects, ladj, eligible = prior.seed_instance(
                adj, seed_mask
            )
            if tuple(rebuilt_seed) != seed:
                raise ValueError("witness seed order mismatch")
            zmask = int(witness["zmask"])
            if zmask not in prior.eligible_covers(ladj, eligible):
                raise ValueError("witness zmask is not an eligible cover")
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                prior.SupportSolver(), prior.ZeroForcingSolver(),
                prior.matching_size(defects), prior.CliqueStructureSolver(),
            )
            if analysis.enhanced_joint_failed:
                raise ValueError("witness cover should have failed earlier")
            nvertices = witness["nvertices"]
            expected_nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            if nvertices != expected_nvertices:
                raise ValueError("witness induced vertex set mismatch")
            graph_n = tuple(prior.induced_graph(adj, nvertices))
            certificate = witness["certificate"]
            clique_mask = sum(1 << vertex for vertex in certificate["clique"])
            if clique_mask not in set(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )):
                raise ValueError("witness clique is not saturating")
            if any(
                strict_h.assess_saturating_clique(graph_n, mask).failed
                for mask in prior.clique_masks(graph_n, analysis.k_rank_upper)
            ):
                raise ValueError("witness cover failed the prior strict-H layer")
            if witness["kind"] == "pure_degree2":
                verify_degree_one.verify_clique_certificate(
                    graph_n, clique_mask, certificate
                )
            elif witness["kind"] == "outside_diagonal_preordering":
                verify_clique_certificate(graph_n, clique_mask, certificate)
            else:
                raise ValueError("unknown witness kind")
            checked += 1
    expected = (
        report["counts"]["pure_degree2_new_failed_covers"]
        + report["counts"]["preordering_new_failed_covers"]
    )
    if checked != expected:
        raise ValueError("certificate count does not match report totals")
    aggregate_fields = (
        "strict_h_passing_covers",
        "covers_without_saturating_clique",
        "covers_with_saturating_clique",
        "baseline_degree1_failed_covers",
        "baseline_degree1_surviving_saturating_clique_covers",
        "pure_degree2_new_failed_covers",
        "preordering_new_failed_covers",
        "final_passing_covers",
    )
    for field in aggregate_fields:
        if field == "covers_with_saturating_clique":
            rebuilt = sum(
                int(item["strict_h_passing_covers"])
                - int(item["covers_without_saturating_clique"])
                for item in report["per_graph"]
            )
        elif field == "baseline_degree1_surviving_saturating_clique_covers":
            rebuilt = sum(
                int(item["strict_h_passing_covers"])
                - int(item["covers_without_saturating_clique"])
                - int(item["baseline_degree1_failed_covers"])
                for item in report["per_graph"]
            )
        else:
            rebuilt = sum(int(item[field]) for item in report["per_graph"])
        if report["counts"][field] != rebuilt:
            raise ValueError(f"aggregate field {field} does not add up")
    return {"graphs": len(report["per_graph"]), "certificates": checked}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample", type=Path, default=Path("d6_k7_rank_sample.json")
    )
    parser.add_argument(
        "--report", type=Path, default=Path("d6_k7_dual_degree2_pilot_report.json")
    )
    args = parser.parse_args()
    summary = verify_report(args.sample, args.report)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
