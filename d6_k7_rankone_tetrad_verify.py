#!/usr/bin/env python3
"""Independent exact checker for the K7 rank-one Schur tetrad pilot.

The checker never imports the new LP locator.  It rebuilds the off-diagonal
Schur numerators from the frozen independent degree-one checker, constructs
all tetrads and triangle signs, verifies every rational identity, and checks
that each witness belongs to the claimed seed, cover, and near-saturating
clique.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as prior
import verify_d6_k7_positive_polynomial_dual as schur_checker


Q = Fraction
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Q]
ROOT = Path(__file__).resolve().parent


def resolve_source(base: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    local = base / path
    return local if local.exists() else ROOT / path


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


def subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    output = dict(left)
    for exponent, value in right.items():
        add_term(output, exponent, -value)
    return output


def reconstruct_rank_one_system(
    adj: Sequence[int], clique_mask: int
) -> tuple[
    tuple[int, ...],
    tuple[int, ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, int], ...],
    tuple[tuple[tuple[int, int], ...], ...],
    tuple[Polynomial, ...],
    tuple[tuple[int, int, int], ...],
    tuple[Polynomial, ...],
]:
    (
        clique,
        remainder,
        masks,
        pairs,
        _targets,
        schur_equations,
    ) = schur_checker.schur_polynomial_system(tuple(adj), clique_mask)
    pair_polynomials = tuple(
        {exponent: -value for exponent, value in equation.items()}
        for equation in schur_equations
    )
    pair_by_label = dict(zip(pairs, pair_polynomials))
    tetrad_labels = []
    tetrads = []
    for a, b, c, d in itertools.combinations(remainder, 4):
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
    signs = []
    for a, b, c in itertools.combinations(remainder, 3):
        sign = multiply(
            multiply(pair_by_label[(a, b)], pair_by_label[(a, c)]),
            pair_by_label[(b, c)],
        )
        if sign:
            triangle_labels.append((a, b, c))
            signs.append(sign)
    return (
        clique,
        remainder,
        masks,
        pairs,
        tuple(tetrad_labels),
        tuple(tetrads),
        tuple(triangle_labels),
        tuple(signs),
    )


def polynomial_degree(polynomial: Polynomial) -> int:
    return max(map(sum, polynomial), default=-1)


def verify_clique_certificate(
    adj: Sequence[int], clique_mask: int, certificate: dict
) -> None:
    (
        clique,
        remainder,
        masks,
        pairs,
        tetrad_labels,
        equations,
        triangle_labels,
        all_signs,
    ) = reconstruct_rank_one_system(adj, clique_mask)
    if certificate.get("schema") != 1:
        raise ValueError("unexpected certificate schema")
    if certificate.get("kind") != (
        "rank_one_schur_tetrad_positive_polynomial"
    ):
        raise ValueError("unexpected certificate kind")
    if certificate.get("clique") != list(clique):
        raise ValueError("certificate clique mismatch")
    if certificate.get("remainder") != list(remainder):
        raise ValueError("certificate remainder mismatch")
    if certificate.get("basis_neighbour_masks") != [list(mask) for mask in masks]:
        raise ValueError("certificate basis masks mismatch")
    if certificate.get("pair_labels") != [list(pair) for pair in pairs]:
        raise ValueError("certificate pair labels mismatch")
    if certificate.get("tetrad_labels") != [
        [list(pair) for pair in label] for label in tetrad_labels
    ]:
        raise ValueError("certificate tetrad labels mismatch")
    if certificate.get("triangle_labels") != [
        list(label) for label in triangle_labels
    ]:
        raise ValueError("certificate triangle labels mismatch")
    output_degree = certificate.get("output_degree")
    if not isinstance(output_degree, int) or output_degree < 4:
        raise ValueError("invalid output degree")
    signs = all_signs if output_degree >= 6 else ()
    variables = len(clique)

    left: Polynomial = {}
    for item in certificate.get("equation_multipliers", []):
        equation_index = item.get("equation")
        if not isinstance(equation_index, int) or not 0 <= equation_index < len(
            equations
        ):
            raise ValueError("tetrad equation index out of range")
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid tetrad multiplier")
        monomial = tuple(monomial_value)
        if (
            len(monomial) != variables
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) + polynomial_degree(equations[equation_index])
            > output_degree
        ):
            raise ValueError("tetrad multiplier exceeds degree bound")
        coefficient = parse_fraction(item.get("coefficient"))
        for exponent, value in shift(equations[equation_index], monomial).items():
            add_term(left, exponent, coefficient * value)

    right: Polynomial = {}
    for item in certificate.get("sign_multipliers", []):
        sign_index = item.get("sign")
        if not isinstance(sign_index, int) or not 0 <= sign_index < len(signs):
            raise ValueError("triangle sign index out of range")
        monomial_value = item.get("monomial")
        if not isinstance(monomial_value, list):
            raise ValueError("invalid triangle-sign multiplier")
        monomial = tuple(monomial_value)
        if (
            len(monomial) != variables
            or not all(isinstance(value, int) and value >= 0 for value in monomial)
            or sum(monomial) + polynomial_degree(signs[sign_index])
            > output_degree
        ):
            raise ValueError("triangle-sign multiplier exceeds degree bound")
        coefficient = parse_fraction(item.get("coefficient"))
        if coefficient <= 0:
            raise ValueError("triangle-sign multiplier is not positive")
        for exponent, value in shift(signs[sign_index], monomial).items():
            add_term(right, exponent, coefficient * value)

    positive: Polynomial = {}
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
            raise ValueError("positive-polynomial coefficient is not positive")
        add_term(positive, monomial, coefficient)
        add_term(right, monomial, coefficient)
    if not positive:
        raise ValueError("certificate has no strict positive polynomial")
    if left != right:
        raise ValueError("exact rank-one identity does not check")


def _verify_current_sample_sources(
    sample_path: Path, sample: dict
) -> dict[int, dict]:
    sample_sources = sample.get("source")
    if not isinstance(sample_sources, dict):
        raise ValueError("sample has no source map")
    resolved_sources = {}
    for name, record in sample_sources.items():
        path = resolve_source(sample_path.parent, record["path"])
        if file_sha256(path) != record["sha256"]:
            raise ValueError(f"sample source hash mismatch: {name}")
        resolved_sources[name] = path
    selection = json.loads(
        resolved_sources["positive_dual_selection"].read_text(
            encoding="utf-8"
        )
    )
    selected = {int(value) for value in selection["selected_indices"]}
    rank_payload = json.loads(
        resolved_sources["rank_input"].read_text(encoding="utf-8")
    )
    rank_graph_by_index = {
        int(graph["index"]): graph for graph in rank_payload["graphs"]
    }
    candidates = []
    decision_by_index = {}
    with gzip.open(
        resolved_sources["positive_dual_decisions"],
        "rt", encoding="utf-8", newline="",
    ) as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            decision_by_index[index] = row
            if (
                row["status"] == "SURVIVOR"
                and index in selected
                and int(row["dual_passing_covers"]) == 1
                and int(row["strict_h_passing_cliques"]) == 0
            ):
                candidates.append((int(row["covers"]), index))
    candidates.sort()
    expected_sample_indices = [
        index for _covers, index in candidates[:int(sample["count"])]
    ]
    if sample["indices"] != expected_sample_indices:
        raise ValueError("sample does not follow its deterministic rule")
    compact_graph_by_index = {
        int(graph["index"]): graph for graph in sample["graphs"]
    }
    if set(compact_graph_by_index) != set(expected_sample_indices):
        raise ValueError("compact sample graph set mismatch")
    for index, graph in compact_graph_by_index.items():
        if graph != rank_graph_by_index[index]:
            raise ValueError("compact sample graph differs from rank input")
    return decision_by_index


def verify_report(sample_path: Path, report_path: Path) -> dict:
    sample = json.loads(sample_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    current_sample = sample.get("kind") == (
        "d6_k7_rankone_tetrad_current_residue_sample"
    )
    if sample.get("schema") != 1 or not isinstance(sample.get("graphs"), list):
        raise ValueError("unexpected sample schema")
    decision_by_index = (
        _verify_current_sample_sources(sample_path, sample)
        if current_sample else {}
    )
    if report.get("schema") != 1 or report.get("kind") != (
        "d6_k7_rankone_tetrad_pilot"
    ):
        raise ValueError("unexpected report schema")
    if report["input"]["sha256"] != file_sha256(sample_path):
        raise ValueError("sample hash mismatch")
    source = report_path.parent / report["source"]["path"]
    if not source.exists():
        source = ROOT / report["source"]["path"]
    if report["source"]["sha256"] != file_sha256(source):
        raise ValueError("pilot source hash mismatch")
    graph_by_index = {
        int(graph["index"]): graph for graph in sample["graphs"]
    }
    recorded_indices = [int(item["index"]) for item in report["per_graph"]]
    if recorded_indices != report["input"]["indices"]:
        raise ValueError("selected graph order mismatch")
    checked = 0
    for recorded in report["per_graph"]:
        graph = graph_by_index[int(recorded["index"])]
        if current_sample:
            decision = decision_by_index[int(recorded["index"])]
            if not (
                decision["status"] == "SURVIVOR"
                and int(decision["dual_passing_covers"]) == 1
                and int(decision["strict_h_passing_cliques"]) == 0
            ):
                raise ValueError("reported graph is outside the residue stratum")
            if not (
                int(recorded["prior_passing_covers"]) == 1
                and int(recorded["no_saturating_clique_covers"]) == 1
                and int(recorded["tetrad_failed_covers"]) == 1
                and int(recorded["final_passing_covers"]) == 0
                and recorded["marginal_graph_rejected"] is True
                and len(recorded["certificates"]) == 1
            ):
                raise ValueError(
                    "pilot did not eliminate the unique prior-passing cover"
                )
        adj = tuple(graph["adjacency"])
        prior.validate_graph(adj)
        for witness in recorded["certificates"]:
            if witness.get("stage") not in (
                "degree4_tetrads", "degree6_tetrads_and_signs"
            ):
                raise ValueError("unknown rank-one certificate stage")
            seed = tuple(int(vertex) for vertex in witness["seed"])
            seed_mask = sum(1 << vertex for vertex in seed)
            if seed_mask not in set(prior.clique_masks(adj, 7)):
                raise ValueError("witness seed is not a K7")
            rebuilt_seed, outside, defects, ladj, eligible = prior.seed_instance(
                adj, seed_mask
            )
            if tuple(rebuilt_seed) != seed:
                raise ValueError("witness seed order mismatch")
            zmask = int(witness["zmask"])
            if zmask not in prior.eligible_covers(ladj, eligible):
                raise ValueError("witness is not an eligible cover")
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                prior.SupportSolver(), prior.ZeroForcingSolver(),
                prior.matching_size(defects), prior.CliqueStructureSolver(),
            )
            if analysis.enhanced_joint_failed:
                raise ValueError("witness cover failed a prior rank layer")
            if int(witness["k_rank_upper"]) != analysis.k_rank_upper:
                raise ValueError("witness K-rank upper bound mismatch")
            expected_nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            if witness["nvertices"] != expected_nvertices:
                raise ValueError("witness induced vertex set mismatch")
            graph_n = tuple(prior.induced_graph(adj, expected_nvertices))
            saturating = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            if saturating:
                raise ValueError("rank-one witness cover has a saturating clique")
            certificate = witness["certificate"]
            if (
                witness["stage"] == "degree4_tetrads"
                and certificate.get("output_degree") != 4
            ):
                raise ValueError("degree-four stage has the wrong certificate")
            if (
                witness["stage"] == "degree6_tetrads_and_signs"
                and certificate.get("output_degree", 0) < 6
            ):
                raise ValueError("degree-six stage has the wrong certificate")
            clique_mask = sum(1 << vertex for vertex in certificate["clique"])
            if clique_mask not in set(prior.clique_masks(
                graph_n, analysis.k_rank_upper - 1
            )):
                raise ValueError("certificate clique is not near-saturating")
            verify_clique_certificate(graph_n, clique_mask, certificate)
            checked += 1
    if checked != report["counts"]["tetrad_failed_covers"]:
        raise ValueError("certificate count differs from cover-failure count")
    aggregate_fields = (
        "prior_passing_covers",
        "no_saturating_clique_covers",
        "covers_with_near_clique",
        "covers_without_near_clique",
        "near_cliques_tested",
        "near_cliques_tested_degree4",
        "near_cliques_tested_degree6",
        "degree4_failed_covers",
        "degree6_sign_new_failed_covers",
        "tetrad_failed_covers",
        "final_passing_covers",
    )
    for field in aggregate_fields:
        rebuilt = sum(int(item[field]) for item in report["per_graph"])
        if report["counts"][field] != rebuilt:
            raise ValueError(f"aggregate field {field} does not add up")
    if report["counts"]["tetrad_failed_covers"] != (
        report["counts"]["degree4_failed_covers"]
        + report["counts"]["degree6_sign_new_failed_covers"]
    ):
        raise ValueError("rank-one stage counts do not partition failures")
    return {"graphs": len(recorded_indices), "certificates": checked}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample", type=Path,
        default=Path("d6_k7_rankone_tetrad_sample.json"),
    )
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_rankone_tetrad_pilot_report.json"),
    )
    args = parser.parse_args()
    print(json.dumps(
        verify_report(args.sample, args.report), sort_keys=True
    ))


if __name__ == "__main__":
    main()
