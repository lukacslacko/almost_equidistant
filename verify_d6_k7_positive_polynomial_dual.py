#!/usr/bin/env python3
"""Independent exact checker for the K7 positive-polynomial dual report.

The checker never invokes the floating-point LP locator.  It rebuilds every
Schur polynomial over ``Fraction``, checks every stored polynomial identity,
and independently repeats the seed/cover quantifiers around the frozen prior
exact K7 and strict-H layers.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

import d6_k7_rank_reference as prior
import d6_k7_special_h_reference as strict_h


ROOT = Path(__file__).resolve().parent
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


def mask_bits(mask: int) -> tuple[int, ...]:
    output: list[int] = []
    while mask:
        bit = mask & -mask
        output.append(bit.bit_length() - 1)
        mask ^= bit
    return tuple(output)


def schur_polynomial_system(
    adj: tuple[int, ...], clique_mask: int
) -> tuple[
    tuple[int, ...],
    tuple[int, ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, int], ...],
    tuple[int, ...],
    tuple[Polynomial, ...],
]:
    """Independently reconstruct the Sherman--Morrison equations."""

    clique = mask_bits(clique_mask)
    if not clique:
        raise ValueError("empty certificate clique")
    if any(
        not (adj[first] & (1 << second))
        for first, second in itertools.combinations(clique, 2)
    ):
        raise ValueError("certificate core is not a required clique")
    remainder = tuple(
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    )
    masks = tuple(
        tuple(int(bool(adj[vertex] & (1 << basis))) for basis in clique)
        for vertex in remainder
    )
    zero = (0,) * len(clique)
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

        # (1 + sum x_i) * w(A intersect B).
        for coordinate in range(len(clique)):
            if not (left[coordinate] and right[coordinate]):
                continue
            linear = list(zero)
            linear[coordinate] = 1
            add_term(polynomial, tuple(linear), Q(1))
            for other in range(len(clique)):
                quadratic = linear.copy()
                quadratic[other] += 1
                add_term(polynomial, tuple(quadratic), Q(1))

        # -w(A)w(B).
        for first in range(len(clique)):
            if not left[first]:
                continue
            for second in range(len(clique)):
                if not right[second]:
                    continue
                quadratic = list(zero)
                quadratic[first] += 1
                quadratic[second] += 1
                add_term(polynomial, tuple(quadratic), Q(-1))

        # -k_AB(1 + sum x_i).
        if target:
            add_term(polynomial, zero, Q(-1))
            for coordinate in range(len(clique)):
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


def verify_clique_certificate(
    adj: tuple[int, ...], clique_mask: int, certificate: dict
) -> None:
    """Independently expand and check one rational polynomial identity."""

    clique, remainder, masks, pairs, targets, equations = (
        schur_polynomial_system(adj, clique_mask)
    )
    if certificate.get("schema") != 1:
        raise ValueError("unexpected certificate schema")
    if certificate.get("kind") != "positive_orthant_polynomial_farkas":
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
    multiplier_degree = certificate.get("multiplier_degree")
    if not isinstance(multiplier_degree, int) or multiplier_degree < 0:
        raise ValueError("invalid multiplier degree")

    combined: Polynomial = {}
    for item in certificate.get("multipliers", []):
        equation_index = item.get("equation")
        if not isinstance(equation_index, int) or not 0 <= equation_index < len(
            equations
        ):
            raise ValueError("certificate equation index out of range")
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid multiplier monomial")
        monomial = tuple(monomial_value)
        if (
            len(monomial) != len(clique)
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) > multiplier_degree
        ):
            raise ValueError("invalid multiplier monomial")
        coefficient = parse_fraction(item.get("coefficient"))
        for exponent, value in equations[equation_index].items():
            shifted = tuple(
                left + right for left, right in zip(exponent, monomial)
            )
            add_term(combined, shifted, coefficient * value)

    claimed: Polynomial = {}
    for item in certificate.get("positive_polynomial", []):
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid positive-polynomial monomial")
        monomial = tuple(monomial_value)
        if (
            len(monomial) != len(clique)
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) > multiplier_degree + 2
        ):
            raise ValueError("invalid positive-polynomial monomial")
        coefficient = parse_fraction(item.get("coefficient"))
        if coefficient <= 0:
            raise ValueError("nonpositive claimed coefficient")
        add_term(claimed, monomial, coefficient)
    if not claimed:
        raise ValueError("claimed positive polynomial is zero")
    if combined != claimed:
        raise ValueError("rational polynomial identity does not check")


def verify_graph(graph: dict, recorded: dict) -> dict:
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    witnesses = {
        (tuple(item["seed"]), int(item["zmask"])): item
        for item in recorded["dual_failure_witnesses"]
    }
    if len(witnesses) != len(recorded["dual_failure_witnesses"]):
        raise AssertionError("duplicate cover witness")
    consumed: set[tuple[tuple[int, ...], int]] = set()
    counts = {
        "seeds": 0,
        "covers": 0,
        "enhanced_passing_covers": 0,
        "strict_h_passing_covers": 0,
        "dual_passing_covers": 0,
        "strict_h_passing_cliques": 0,
        "dual_failing_cliques": 0,
    }
    strict_rejected = dual_rejected = False
    for seed_mask in prior.clique_masks(adj, 7):
        counts["seeds"] += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(
            adj, seed_mask
        )
        total_term_rank = prior.matching_size(defects)
        seed_strict_passes = seed_dual_passes = 0
        covers = prior.eligible_covers(ladj, eligible)
        counts["covers"] += len(covers)
        for zmask in covers:
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            counts["enhanced_passing_covers"] += 1
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
            counts["strict_h_passing_covers"] += 1
            counts["strict_h_passing_cliques"] += len(clique_masks)
            seed_strict_passes += 1
            key = (tuple(seed), zmask)
            witness = witnesses.get(key)
            if witness is None:
                counts["dual_passing_covers"] += 1
                seed_dual_passes += 1
                continue
            if witness["nvertices"] != nvertices:
                raise AssertionError("witness induced-vertex list differs")
            certificate = witness["certificate"]
            clique_mask = sum(1 << vertex for vertex in certificate["clique"])
            if clique_mask not in clique_masks:
                raise AssertionError("witness clique is not saturating")
            verify_clique_certificate(tuple(graph_n), clique_mask, certificate)
            counts["dual_failing_cliques"] += 1
            consumed.add(key)
        if not seed_strict_passes:
            strict_rejected = True
        if not seed_dual_passes:
            dual_rejected = True
    if consumed != set(witnesses):
        raise AssertionError("orphan or inapplicable dual witness")
    expected = {
        **counts,
        "strict_h_rejected": strict_rejected,
        "dual_rejected": dual_rejected,
        "marginal_dual_rejected": dual_rejected and not strict_rejected,
    }
    for key, value in expected.items():
        if recorded[key] != value:
            raise AssertionError(
                f"graph {graph.get('index')} field {key}: "
                f"recorded {recorded[key]!r}, rebuilt {value!r}"
            )
    return expected


def verify_report(sample_path: Path, selection_path: Path, report_path: Path) -> None:
    sample = json.loads(sample_path.read_text(encoding="utf-8"))
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report["input"]["sample"]["sha256"] != file_sha256(sample_path):
        raise AssertionError("sample hash mismatch")
    if report["input"]["selection_report"]["sha256"] != file_sha256(
        selection_path
    ):
        raise AssertionError("selection hash mismatch")
    source = ROOT / report["source"]["file"]
    if report["source"]["sha256"] != file_sha256(source):
        raise AssertionError("dual source hash mismatch")
    selected = set(selection["graph_decisions"]["combined_survivors"])
    if selected != set(report["input"]["selected_indices"]):
        raise AssertionError("selected-index set mismatch")
    graph_by_index = {graph["index"]: graph for graph in sample["graphs"]}
    record_by_index = {item["index"]: item for item in report["per_graph"]}
    if set(record_by_index) != selected:
        raise AssertionError("per-graph report coverage mismatch")

    rebuilt = [
        verify_graph(graph_by_index[index], record_by_index[index])
        for index in sorted(selected)
    ]
    rejected = sorted(
        index for index in selected
        if record_by_index[index]["marginal_dual_rejected"]
    )
    survivors = sorted(selected - set(rejected))
    if report["decisions"]["marginal_dual_rejected"] != rejected:
        raise AssertionError("marginal rejection list mismatch")
    if report["decisions"]["dual_survivors"] != survivors:
        raise AssertionError("survivor list mismatch")
    aggregate_mapping = {
        "seeds": "seeds",
        "covers": "covers",
        "enhanced_passing_covers": "enhanced_passing_covers",
        "strict_H_passing_covers": "strict_h_passing_covers",
        "dual_passing_covers": "dual_passing_covers",
        "strict_H_passing_cliques": "strict_h_passing_cliques",
        "dual_failing_cliques_first_per_cover": "dual_failing_cliques",
    }
    for report_key, rebuilt_key in aggregate_mapping.items():
        value = sum(item[rebuilt_key] for item in rebuilt)
        if report["counts"][report_key] != value:
            raise AssertionError(f"aggregate count {report_key} mismatch")
    if report["counts"]["marginal_dual_rejected"] != len(rejected):
        raise AssertionError("aggregate marginal rejection count mismatch")
    if report["counts"]["dual_survivors"] != len(survivors):
        raise AssertionError("aggregate survivor count mismatch")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample", type=Path, default=ROOT / "d6_k7_rank_sample.json"
    )
    parser.add_argument(
        "--selection", type=Path, default=ROOT / "d6_k7_strict_h_report.json"
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k7_positive_polynomial_dual_report.json",
    )
    args = parser.parse_args()
    verify_report(args.sample, args.selection, args.report)
    print("positive-polynomial dual report verified")


if __name__ == "__main__":
    main()
