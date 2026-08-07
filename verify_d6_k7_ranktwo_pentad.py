#!/usr/bin/env python3
"""Independent symbolic and artifact checker for rank-two Schur pentads."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence

import d6_k7_rank_reference as prior
import verify_d6_k7_positive_polynomial_dual as schur_checker
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
Q = Fraction
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Q]
Edge = tuple[int, int]

EXPECTED_RANK_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_SELECTION_SHA256 = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
EXPECTED_PRIOR_CERTIFICATES_SHA256 = (
    "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
)
EXPECTED_UNION_SHA256 = (
    "1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599"
)


# Independently transcribed expanded expression.  This verifier never imports
# the production pentad module or its polynomial helpers.
TERMS: tuple[tuple[int, tuple[Edge, ...]], ...] = (
    (+1, ((0, 1), (0, 2), (1, 3), (2, 4), (3, 4))),
    (-1, ((0, 1), (0, 2), (1, 4), (2, 3), (3, 4))),
    (-1, ((0, 1), (0, 3), (1, 2), (2, 4), (3, 4))),
    (+1, ((0, 1), (0, 3), (1, 4), (2, 3), (2, 4))),
    (+1, ((0, 1), (0, 4), (1, 2), (2, 3), (3, 4))),
    (-1, ((0, 1), (0, 4), (1, 3), (2, 3), (2, 4))),
    (+1, ((0, 2), (0, 3), (1, 2), (1, 4), (3, 4))),
    (-1, ((0, 2), (0, 3), (1, 3), (1, 4), (2, 4))),
    (-1, ((0, 2), (0, 4), (1, 2), (1, 3), (3, 4))),
    (+1, ((0, 2), (0, 4), (1, 3), (1, 4), (2, 3))),
    (+1, ((0, 3), (0, 4), (1, 2), (1, 3), (2, 4))),
    (-1, ((0, 3), (0, 4), (1, 2), (1, 4), (2, 3))),
)


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


def add_term(target: Polynomial, exponent: Exponent, value: Q) -> None:
    if not value:
        return
    updated = target.get(exponent, Q(0)) + value
    if updated:
        target[exponent] = updated
    else:
        target.pop(exponent, None)


def add(left: Polynomial, right: Polynomial, scale: Q = Q(1)) -> Polynomial:
    output = dict(left)
    for exponent, value in right.items():
        add_term(output, exponent, scale * value)
    return output


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = tuple(
                a + b for a, b in zip(left_exponent, right_exponent)
            )
            add_term(output, exponent, left_value * right_value)
    return output


def product(polynomials: Sequence[Polynomial], variables: int) -> Polynomial:
    output = {(0,) * variables: Q(1)}
    for item in polynomials:
        output = multiply(output, item)
    return output


def variable(index: int, variables: int) -> Polynomial:
    exponent = [0] * variables
    exponent[index] = 1
    return {tuple(exponent): Q(1)}


def pair(first: int, second: int) -> Edge:
    return (first, second) if first < second else (second, first)


def independent_pentad(
    vertices: Sequence[int],
    pairs: Mapping[Edge, Polynomial],
    variables: int,
) -> Polynomial:
    vertices = tuple(int(vertex) for vertex in vertices)
    if len(vertices) != 5 or len(set(vertices)) != 5:
        raise ValueError("pentad needs five distinct vertices")
    output: Polynomial = {}
    for sign, local_edges in TERMS:
        factors = [
            pairs[pair(vertices[left], vertices[right])]
            for left, right in local_edges
        ]
        for exponent, value in product(factors, variables).items():
            add_term(output, exponent, Q(sign) * value)
    return output


def verify_elimination_derivation() -> dict:
    """Derive the pentad by eliminating t from two residual tetrads."""

    labels = tuple(itertools.combinations(range(5), 2))
    pair_index = {edge: index for index, edge in enumerate(labels)}
    g = {edge: variable(index, len(labels)) for edge, index in pair_index.items()}
    mul = lambda *items: product(items, len(labels))
    # With t=||v_0||^2, residual rank-one tetrads reduce (after removing
    # their common factor t) to t*A+B=0 and t*C+D=0.
    A = add(mul(g[(1, 2)], g[(3, 4)]), mul(g[(1, 3)], g[(2, 4)]), Q(-1))
    B: Polynomial = {}
    for sign, edges in (
        (-1, ((0, 1), (0, 2), (3, 4))),
        (+1, ((0, 1), (0, 3), (2, 4))),
        (+1, ((0, 2), (0, 4), (1, 3))),
        (-1, ((0, 3), (0, 4), (1, 2))),
    ):
        B = add(B, mul(*(g[edge] for edge in edges)), Q(sign))
    C = add(mul(g[(1, 2)], g[(3, 4)]), mul(g[(1, 4)], g[(2, 3)]), Q(-1))
    D: Polynomial = {}
    for sign, edges in (
        (-1, ((0, 1), (0, 2), (3, 4))),
        (+1, ((0, 1), (0, 4), (2, 3))),
        (+1, ((0, 2), (0, 3), (1, 4))),
        (-1, ((0, 3), (0, 4), (1, 2))),
    ):
        D = add(D, mul(*(g[edge] for edge in edges)), Q(sign))
    eliminated = add(multiply(A, D), multiply(C, B), Q(-1))
    expanded = independent_pentad(range(5), g, len(labels))
    if eliminated != expanded:
        raise ValueError("elimination expansion differs from the pentad")
    if len(expanded) != 12 or any(sum(exponent) != 5 for exponent in expanded):
        raise ValueError("pentad does not have twelve degree-five terms")
    return {"formal_variables": 10, "terms": len(expanded), "degree": 5}


def verify_generic_rank_two_gram_identity() -> dict:
    """Substitute five generic two-dimensional vectors and expand to zero."""

    variables = 10  # x_0,y_0,...,x_4,y_4
    coordinates = [
        (variable(2 * index, variables), variable(2 * index + 1, variables))
        for index in range(5)
    ]
    gram: dict[Edge, Polynomial] = {}
    for first, second in itertools.combinations(range(5), 2):
        gram[(first, second)] = add(
            multiply(coordinates[first][0], coordinates[second][0]),
            multiply(coordinates[first][1], coordinates[second][1]),
        )
    expanded = independent_pentad(range(5), gram, variables)
    if expanded:
        raise ValueError("generic two-dimensional Gram pentad is nonzero")
    return {
        "coordinate_variables": variables,
        "formal_products_before_cancellation": 12 * (2 ** 5),
        "residual_terms": 0,
    }


def verify_pentad_permutation_orbit() -> dict:
    """Check all 120 labelings give only the same pentad up to sign."""

    labels = tuple(itertools.combinations(range(5), 2))
    formal = {
        edge: variable(index, len(labels))
        for index, edge in enumerate(labels)
    }
    base = independent_pentad(range(5), formal, len(labels))
    negative = {exponent: -value for exponent, value in base.items()}
    counts = Counter()
    for permutation in itertools.permutations(range(5)):
        image = independent_pentad(permutation, formal, len(labels))
        if image == base:
            counts["same"] += 1
        elif image == negative:
            counts["opposite"] += 1
        else:
            raise ValueError(
                f"pentad labeling {permutation} is not equal up to sign"
            )
    if sum(counts.values()) != 120:
        raise ValueError("pentad permutation orbit is incomplete")
    return {
        "labelings_checked": 120,
        "same_sign": counts["same"],
        "opposite_sign": counts["opposite"],
        "genuinely_distinct_up_to_sign": 1,
    }


def parse_fraction(value: object) -> Q:
    if (
        not isinstance(value, list) or len(value) != 2
        or not all(isinstance(item, int) for item in value)
    ):
        raise ValueError("invalid serialized fraction")
    return Q(value[0], value[1])


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def reconstruct_full_targets(
    report_path: Path, report: dict
) -> tuple[list[dict], dict[int, dict]]:
    """Independently replay the frozen input/seed/cover target quantifiers."""

    provenance = report["provenance"]
    rank_path = resolve(report_path.parent, provenance["rank_input"]["path"])
    selection_path = resolve(
        report_path.parent, provenance["selection"]["path"]
    )
    union_path = resolve(report_path.parent, provenance["union"]["path"])
    observed = {
        "rank_input": sha256(rank_path),
        "selection": sha256(selection_path),
        "union": sha256(union_path),
    }
    expected = {
        "rank_input": EXPECTED_RANK_INPUT_SHA256,
        "selection": EXPECTED_SELECTION_SHA256,
        "union": EXPECTED_UNION_SHA256,
    }
    if observed != expected:
        raise ValueError("full pentad frozen input hash mismatch")
    if any(provenance[name]["sha256"] != observed[name] for name in observed):
        raise ValueError("full pentad report provenance mismatch")

    rank_payload = json.loads(rank_path.read_text(encoding="utf-8"))
    graphs = rank_payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != 17_764:
        raise ValueError("rank input graph count mismatch")
    graph_by_index: dict[int, dict] = {}
    for graph in graphs:
        index = int(graph["index"])
        if index in graph_by_index:
            raise ValueError("rank input repeats a graph")
        adjacency = tuple(int(row) for row in graph["adjacency"])
        prior.validate_graph(adjacency)
        graph_by_index[index] = graph

    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    selected = [int(index) for index in selection["selected_indices"]]
    if (
        len(selected) != 12_839 or len(set(selected)) != len(selected)
        or selection["selected_indices_sha256"] != stable_hash(selected)
        or selection["rank_input"]["sha256"] != observed["rank_input"]
    ):
        raise ValueError("tetrad selection mismatch")
    certificate_path = resolve(
        selection_path.parent,
        selection["prior_degree_one"]["certificates"]["path"],
    )
    if sha256(certificate_path) != EXPECTED_PRIOR_CERTIFICATES_SHA256:
        raise ValueError("prior certificate archive hash mismatch")
    if provenance["prior_certificates"]["sha256"] != sha256(certificate_path):
        raise ValueError("prior certificate provenance mismatch")
    witness_keys: dict[int, set[tuple[tuple[int, ...], int]]] = {}
    with gzip.open(certificate_path, "rt", encoding="utf-8") as stream:
        for line in stream:
            entry = json.loads(line)
            witness_keys[int(entry["index"])] = {
                (tuple(int(vertex) for vertex in item["seed"]), int(item["zmask"]))
                for item in entry["dual_failure_witnesses"]
            }

    union = json.loads(union_path.read_text(encoding="utf-8"))
    if (
        union.get("schema") != "d6-k7-rankone-pattern-union-v1"
        or union.get("status") != "COMPLETE"
    ):
        raise ValueError("union schema/status mismatch")
    residue = [int(index) for index in union["sets"]["exact_residue"]["indices"]]
    if union["sets"]["exact_residue"]["indices_sha256"] != stable_hash(residue):
        raise ValueError("union residue hash mismatch")
    profile_list = union["cover_structure"]["residue_profiles"]
    if [int(item["index"]) for item in profile_list] != residue:
        raise ValueError("union residue profile ordering mismatch")
    profiles = {int(item["index"]): item for item in profile_list}
    graph_indices = [
        int(item["index"])
        for item in profile_list if item["has_no_near_clique_cover"]
    ]
    if len(graph_indices) != 50 or any(index not in set(selected) for index in graph_indices):
        raise ValueError("no-near residue graph population mismatch")

    targets: list[dict] = []
    for graph_index in graph_indices:
        adjacency = tuple(int(row) for row in graph_by_index[graph_index]["adjacency"])
        support_solver = prior.SupportSolver()
        zero_forcing = prior.ZeroForcingSolver()
        clique_solver = prior.CliqueStructureSolver()
        for seed_mask in prior.clique_masks(adjacency, 7):
            seed, outside, defects, ladj, eligible = prior.seed_instance(
                adjacency, seed_mask
            )
            total_term_rank = prior.matching_size(defects)
            for zmask in prior.eligible_covers(ladj, eligible):
                analysis = prior.analyze_cover(
                    adjacency, outside, defects, zmask,
                    support_solver, zero_forcing, total_term_rank, clique_solver,
                )
                if analysis.enhanced_joint_failed:
                    continue
                nvertices = [
                    outside[position] for position in range(len(outside))
                    if not (zmask & (1 << position))
                ]
                graph_n = prior.induced_graph(adjacency, nvertices)
                upper = int(analysis.k_rank_upper)
                if upper and next(prior.clique_masks(graph_n, upper), None) is not None:
                    continue
                key = (tuple(seed), int(zmask))
                if key in witness_keys.get(graph_index, set()):
                    raise ValueError("no-saturating target has a prior dual witness")
                if upper > 1 and next(
                    prior.clique_masks(graph_n, upper - 1), None
                ) is not None:
                    continue
                omega = prior.clique_number(graph_n)
                clique_mask = next(prior.clique_masks(graph_n, omega))
                targets.append({
                    "graph_index": graph_index,
                    "seed": [int(vertex) for vertex in seed],
                    "zmask": int(zmask),
                    "nvertices": [int(vertex) for vertex in nvertices],
                    "graph_n": [int(row) for row in graph_n],
                    "rank_upper": upper,
                    "maximum_clique_size": omega,
                    "fixed_maximum_clique": list(prior.bits(clique_mask)),
                })
    if len(targets) != 53 or len({item["graph_index"] for item in targets}) != 50:
        raise ValueError("independent target replay is not 53 covers/50 graphs")
    return targets, profiles


def independent_schur_pentads(
    graph_n: Sequence[int], clique: Sequence[int]
) -> tuple[
    tuple[int, ...], tuple[tuple[int, ...], ...], tuple[Edge, ...],
    tuple[tuple[int, ...], ...], tuple[Polynomial, ...]
]:
    clique = tuple(int(vertex) for vertex in clique)
    clique_mask = sum(1 << vertex for vertex in clique)
    (
        checked_clique,
        remainder,
        masks,
        pair_labels,
        _targets,
        schur_equations,
    ) = schur_checker.schur_polynomial_system(tuple(graph_n), clique_mask)
    if checked_clique != clique or len(clique) != 5:
        raise ValueError("independent full Schur clique mismatch")
    pair_polynomials = tuple(
        {exponent: -value for exponent, value in equation.items()}
        for equation in schur_equations
    )
    pair_by_label = dict(zip(pair_labels, pair_polynomials))
    labels = tuple(itertools.combinations(remainder, 5))
    equations = tuple(
        independent_pentad(label, pair_by_label, len(clique))
        for label in labels
    )
    return remainder, masks, pair_labels, labels, equations


def verify_report(report_path: Path, expected_sha256: str) -> dict:
    observed_hash = sha256(report_path)
    if observed_hash != expected_sha256:
        raise ValueError(
            f"report hash {observed_hash} != explicit pin {expected_sha256}"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind") != "d6_k7_ranktwo_pentad_first_core_pilot"
        or report.get("status") != "REJECTED"
    ):
        raise ValueError("unexpected report schema/kind/status")

    elimination = verify_elimination_derivation()
    generic_gram = verify_generic_rank_two_gram_identity()
    permutation_orbit = verify_pentad_permutation_orbit()

    upstream = report["upstream"]
    basis_path = resolve(report_path.parent, upstream["basis_pilot"]["path"])
    basis_verification_path = resolve(
        report_path.parent, upstream["basis_verification"]["path"]
    )
    if sha256(basis_path) != upstream["basis_pilot"]["sha256"]:
        raise ValueError("upstream basis pilot hash mismatch")
    if sha256(basis_verification_path) != upstream[
        "basis_verification"
    ]["sha256"]:
        raise ValueError("upstream basis verification hash mismatch")
    basis_report = json.loads(basis_path.read_text(encoding="utf-8"))
    basis_verification = json.loads(
        basis_verification_path.read_text(encoding="utf-8")
    )
    if (
        basis_verification.get("status") != "PASS"
        or basis_verification["report"]["sha256"] != sha256(basis_path)
    ):
        raise ValueError("upstream independent verification gate failed")
    cover = basis_report["covers"][0]
    scope = report["scope"]
    expected_scope = {
        "graph_index": int(cover["graph_index"]),
        "seed": cover["seed"],
        "zmask": int(cover["zmask"]),
        "nvertices": cover["nvertices"],
        "graph_n": cover["graph_n"],
        "rank_upper": int(cover["rank_upper"]),
        "clique": cover["fixed_maximum_clique"],
    }
    for name, value in expected_scope.items():
        if scope[name] != value:
            raise ValueError(f"scope mismatch in {name}")
    if int(scope["cores_tested"]) != 1 or int(scope["covers_tested"]) != 1:
        raise ValueError("pilot is not restricted to exactly one core/cover")

    graph_n = tuple(int(row) for row in scope["graph_n"])
    clique = tuple(int(vertex) for vertex in scope["clique"])
    clique_mask = sum(1 << vertex for vertex in clique)
    (
        checked_clique,
        remainder,
        masks,
        pair_labels,
        _targets,
        schur_equations,
    ) = schur_checker.schur_polynomial_system(graph_n, clique_mask)
    if checked_clique != clique or len(clique) != 5:
        raise ValueError("independent Schur clique mismatch")
    pair_polynomials = tuple(
        {exponent: -value for exponent, value in equation.items()}
        for equation in schur_equations
    )
    pair_by_label = dict(zip(pair_labels, pair_polynomials))
    pentad_labels = tuple(itertools.combinations(remainder, 5))
    equations = tuple(
        independent_pentad(label, pair_by_label, len(clique))
        for label in pentad_labels
    )

    certificate = report["certificate"]
    if (
        certificate.get("schema") != 1
        or certificate.get("kind")
        != "rank_two_schur_pentad_coefficientwise_positive"
    ):
        raise ValueError("unexpected pentad certificate schema/kind")
    if certificate.get("clique") != list(clique):
        raise ValueError("certificate clique mismatch")
    if certificate.get("remainder") != list(remainder):
        raise ValueError("certificate remainder mismatch")
    if certificate.get("basis_neighbour_masks") != [list(mask) for mask in masks]:
        raise ValueError("certificate basis masks mismatch")
    if certificate.get("pair_labels") != [list(edge) for edge in pair_labels]:
        raise ValueError("certificate pair ordering mismatch")
    index = int(certificate["pentad_index"])
    if not 0 <= index < len(equations):
        raise ValueError("certificate pentad index out of range")
    if certificate.get("pentad_vertices") != list(pentad_labels[index]):
        raise ValueError("certificate pentad label mismatch")
    sign = int(certificate["sign"])
    if sign not in (-1, 1):
        raise ValueError("certificate sign is not plus/minus one")
    positive = {
        exponent: Q(sign) * value for exponent, value in equations[index].items()
    }
    if not positive or any(value <= 0 for value in positive.values()):
        raise ValueError("signed pentad is not strictly coefficientwise positive")
    serialized: Polynomial = {}
    for item in certificate["positive_polynomial"]:
        exponent = tuple(int(value) for value in item["monomial"])
        coefficient = parse_fraction(item["coefficient"])
        if (
            len(exponent) != len(clique) or any(value < 0 for value in exponent)
            or coefficient <= 0 or exponent in serialized
        ):
            raise ValueError("invalid serialized positive monomial")
        serialized[exponent] = coefficient
    if serialized != positive:
        raise ValueError("serialized coefficientwise identity mismatch")

    summary = report["certificate_summary"]
    observed_summary = {
        "pentad_index": index,
        "pentad_vertices": list(pentad_labels[index]),
        "sign": sign,
        "positive_terms": len(positive),
        "polynomial_degree": max(map(sum, positive), default=-1),
        "coefficient_min": [min(positive.values()).numerator,
                            min(positive.values()).denominator],
        "coefficient_max": [max(positive.values()).numerator,
                            max(positive.values()).denominator],
        "monomial_degree_counts": {
            str(degree): count for degree, count in sorted(Counter(
                map(sum, positive)
            ).items())
        },
    }
    if summary != observed_summary:
        raise ValueError("certificate summary mismatch")

    return {
        "schema": 1,
        "kind": "d6_k7_ranktwo_pentad_first_core_verification",
        "status": "PASS",
        "report": {"path": str(report_path), "sha256": observed_hash},
        "symbolic_checks": {
            "elimination_derivation": elimination,
            "generic_rank_two_gram_substitution": generic_gram,
            "permutation_orbit": permutation_orbit,
        },
        "checked": {
            "cores": 1,
            "covers": 1,
            "pentads_reconstructed": len(equations),
            "certificate_pentad_index": index,
            "positive_terms": len(positive),
        },
        "claim": (
            "The 12-term pentad was independently derived and expanded to "
            "zero on five generic 2D vectors. The stored exact Sherman "
            "substitution is coefficientwise strictly positive, so this one "
            "rank-at-most-two Schur cover is impossible."
        ),
    }


def classify_equation(equation: Polynomial) -> str:
    if not equation:
        return "ZERO"
    if all(value > 0 for value in equation.values()):
        return "POSITIVE"
    if all(value < 0 for value in equation.values()):
        return "NEGATIVE"
    return "MIXED"


def polynomial_payload(equation: Polynomial) -> list[dict]:
    return [
        {
            "monomial": list(exponent),
            "coefficient": [coefficient.numerator, coefficient.denominator],
        }
        for exponent, coefficient in sorted(equation.items())
    ]


def verify_full_cover_record(target: dict, record: dict, ordinal: int) -> dict:
    """Independently verify one full-scan cover record with bounded lifetime."""

    expected_context = {
        "ordinal": ordinal,
        "graph_index": int(target["graph_index"]),
        "seed": target["seed"],
        "zmask": int(target["zmask"]),
        "nvertices": target["nvertices"],
        "graph_n": target["graph_n"],
        "rank_upper": int(target["rank_upper"]),
        "fixed_maximum_clique": target["fixed_maximum_clique"],
    }
    for name, value in expected_context.items():
        if record[name] != value:
            raise ValueError(f"cover {ordinal} context mismatch in {name}")
    graph_n = tuple(int(row) for row in target["graph_n"])
    clique = tuple(int(vertex) for vertex in target["fixed_maximum_clique"])
    remainder, masks, pairs, labels, equations = independent_schur_pentads(
        graph_n, clique
    )
    if len(equations) != 21 or int(record["pentads_tested"]) != 21:
        raise ValueError("cover does not exhaust all 21 pentads")
    summaries = record["equation_summaries"]
    if len(summaries) != len(equations):
        raise ValueError("pentad summary count mismatch")
    one_sign_indices = []
    local_classes = Counter()
    for index, (vertices, equation, item) in enumerate(zip(
        labels, equations, summaries
    )):
        classification = classify_equation(equation)
        expected_summary = {
            "pentad_index": index,
            "vertices": list(vertices),
            "classification": classification,
            "terms": len(equation),
            "degree": max(map(sum, equation), default=-1),
            "polynomial_sha256": stable_hash(polynomial_payload(equation)),
        }
        if item != expected_summary:
            raise ValueError(
                f"cover {ordinal} pentad {index} summary mismatch"
            )
        local_classes[classification] += 1
        if classification in ("POSITIVE", "NEGATIVE"):
            one_sign_indices.append(index)
    if record["classification_counts"] != dict(local_classes):
        raise ValueError("cover pentad classification counts mismatch")
    if record["one_sign_pentad_indices"] != one_sign_indices:
        raise ValueError("cover one-sign pentad index list mismatch")

    certificate = record["certificate"]
    certificate_checked = False
    if one_sign_indices:
        if record["status"] != "REJECTED" or certificate is None:
            raise ValueError("one-sign cover is not certified rejected")
        index = one_sign_indices[0]
        sign = 1 if classify_equation(equations[index]) == "POSITIVE" else -1
        if (
            certificate.get("schema") != 1
            or certificate.get("kind")
            != "rank_two_schur_pentad_coefficientwise_positive"
            or certificate.get("clique") != list(clique)
            or certificate.get("remainder") != list(remainder)
            or certificate.get("basis_neighbour_masks")
            != [list(mask) for mask in masks]
            or certificate.get("pair_labels") != [list(edge) for edge in pairs]
            or int(certificate["pentad_index"]) != index
            or certificate.get("pentad_vertices") != list(labels[index])
            or int(certificate["sign"]) != sign
        ):
            raise ValueError("full pentad certificate context mismatch")
        positive = {
            exponent: Q(sign) * value
            for exponent, value in equations[index].items()
        }
        if not positive or any(value <= 0 for value in positive.values()):
            raise ValueError("full signed pentad is not strictly positive")
        serialized: Polynomial = {}
        for item in certificate["positive_polynomial"]:
            exponent = tuple(int(value) for value in item["monomial"])
            coefficient = parse_fraction(item["coefficient"])
            if (
                len(exponent) != 5 or any(value < 0 for value in exponent)
                or coefficient <= 0 or exponent in serialized
            ):
                raise ValueError("invalid full positive-polynomial term")
            serialized[exponent] = coefficient
        if serialized != positive:
            raise ValueError("full positive-polynomial identity mismatch")
        if int(certificate["polynomial_degree"]) != max(map(sum, positive)):
            raise ValueError("full certificate degree mismatch")
        certificate_checked = True
    else:
        if record["status"] != "UNRESOLVED" or certificate is not None:
            raise ValueError("no-one-sign cover is not marked unresolved")
    return {
        "ordinal": ordinal,
        "graph_index": int(record["graph_index"]),
        "status": record["status"],
        "pentads_checked": len(equations),
        "classification_counts": dict(local_classes),
        "one_sign_count": len(one_sign_indices),
        "certificate_checked": certificate_checked,
    }


def verify_full_report(report_path: Path, expected_sha256: str) -> dict:
    """Independently reconstruct all 53 covers and every stored identity."""

    observed_hash = sha256(report_path)
    if observed_hash != expected_sha256:
        raise ValueError(
            f"full report hash {observed_hash} != explicit pin {expected_sha256}"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind")
        != "d6_k7_ranktwo_pentad_full_coefficientwise_scan"
        or report.get("status") != "COMPLETE"
    ):
        raise ValueError("unexpected full report schema/kind/status")

    symbolic = {
        "elimination_derivation": verify_elimination_derivation(),
        "generic_rank_two_gram_substitution": (
            verify_generic_rank_two_gram_identity()
        ),
        "permutation_orbit": verify_pentad_permutation_orbit(),
    }
    for name, expected_hash in report["source_hashes"].items():
        path = ROOT / name
        if sha256(path) != expected_hash:
            raise ValueError(f"full report source hash mismatch: {name}")

    targets, profiles = reconstruct_full_targets(report_path, report)
    records = report["covers"]
    if len(records) != len(targets):
        raise ValueError("full report cover record count mismatch")
    status_counts = Counter()
    classification_counts = Counter()
    rejected_by_graph: Counter[int] = Counter()
    cover_counts: Counter[int] = Counter()
    certificates_checked = 0
    pentads_checked = 0
    multiple = 0
    for ordinal, (target, record) in enumerate(zip(targets, records)):
        expected_context = {
            "ordinal": ordinal,
            "graph_index": target["graph_index"],
            "seed": target["seed"],
            "zmask": target["zmask"],
            "nvertices": target["nvertices"],
            "graph_n": target["graph_n"],
            "rank_upper": target["rank_upper"],
            "fixed_maximum_clique": target["fixed_maximum_clique"],
        }
        for name, value in expected_context.items():
            if record[name] != value:
                raise ValueError(
                    f"cover {ordinal} context mismatch in {name}"
                )
        graph_n = tuple(int(row) for row in target["graph_n"])
        clique = tuple(int(vertex) for vertex in target["fixed_maximum_clique"])
        remainder, masks, pairs, labels, equations = independent_schur_pentads(
            graph_n, clique
        )
        if len(equations) != 21 or int(record["pentads_tested"]) != 21:
            raise ValueError("cover does not exhaust all 21 pentads")
        summaries = record["equation_summaries"]
        if len(summaries) != len(equations):
            raise ValueError("pentad summary count mismatch")
        one_sign_indices = []
        local_classes = Counter()
        for index, (vertices, equation, item) in enumerate(zip(
            labels, equations, summaries
        )):
            classification = classify_equation(equation)
            expected_summary = {
                "pentad_index": index,
                "vertices": list(vertices),
                "classification": classification,
                "terms": len(equation),
                "degree": max(map(sum, equation), default=-1),
                "polynomial_sha256": stable_hash(polynomial_payload(equation)),
            }
            if item != expected_summary:
                raise ValueError(
                    f"cover {ordinal} pentad {index} summary mismatch"
                )
            local_classes[classification] += 1
            classification_counts[classification] += 1
            pentads_checked += 1
            if classification in ("POSITIVE", "NEGATIVE"):
                one_sign_indices.append(index)
        if record["classification_counts"] != dict(local_classes):
            raise ValueError("cover pentad classification counts mismatch")
        if record["one_sign_pentad_indices"] != one_sign_indices:
            raise ValueError("cover one-sign pentad index list mismatch")
        if len(one_sign_indices) > 1:
            multiple += 1

        certificate = record["certificate"]
        if one_sign_indices:
            if record["status"] != "REJECTED" or certificate is None:
                raise ValueError("one-sign cover is not certified rejected")
            index = one_sign_indices[0]
            sign = 1 if classify_equation(equations[index]) == "POSITIVE" else -1
            if (
                certificate.get("schema") != 1
                or certificate.get("kind")
                != "rank_two_schur_pentad_coefficientwise_positive"
                or certificate.get("clique") != list(clique)
                or certificate.get("remainder") != list(remainder)
                or certificate.get("basis_neighbour_masks")
                != [list(mask) for mask in masks]
                or certificate.get("pair_labels") != [list(edge) for edge in pairs]
                or int(certificate["pentad_index"]) != index
                or certificate.get("pentad_vertices") != list(labels[index])
                or int(certificate["sign"]) != sign
            ):
                raise ValueError("full pentad certificate context mismatch")
            positive = {
                exponent: Q(sign) * value
                for exponent, value in equations[index].items()
            }
            if not positive or any(value <= 0 for value in positive.values()):
                raise ValueError("full signed pentad is not strictly positive")
            serialized: Polynomial = {}
            for item in certificate["positive_polynomial"]:
                exponent = tuple(int(value) for value in item["monomial"])
                coefficient = parse_fraction(item["coefficient"])
                if (
                    len(exponent) != 5 or any(value < 0 for value in exponent)
                    or coefficient <= 0 or exponent in serialized
                ):
                    raise ValueError("invalid full positive-polynomial term")
                serialized[exponent] = coefficient
            if serialized != positive:
                raise ValueError("full positive-polynomial identity mismatch")
            if int(certificate["polynomial_degree"]) != max(map(sum, positive)):
                raise ValueError("full certificate degree mismatch")
            certificates_checked += 1
            rejected_by_graph[int(record["graph_index"])] += 1
        else:
            if record["status"] != "UNRESOLVED" or certificate is not None:
                raise ValueError("no-one-sign cover is not marked unresolved")
        status_counts[record["status"]] += 1
        cover_counts[int(record["graph_index"])] += 1

    all_no_near_rejected = sorted(
        index for index, count in cover_counts.items()
        if rejected_by_graph[index] == count
    )
    marginal = sorted(
        index for index in all_no_near_rejected
        if not profiles[index]["has_saturating_cover"]
        and not profiles[index]["has_tetrad_resistant_near_clique_cover"]
    )
    summary = report["summary"]
    exact_summary = {
        "cover_status_counts": dict(status_counts),
        "pentads_tested": pentads_checked,
        "pentad_classification_counts": dict(classification_counts),
        "covers_with_multiple_one_sign_pentads": multiple,
        "graphs_with_all_no_near_covers_rejected": len(all_no_near_rejected),
        "all_no_near_covers_rejected_indices": all_no_near_rejected,
        "marginal_graph_rejections": len(marginal),
        "marginal_graph_rejected_indices": marginal,
    }
    for name, value in exact_summary.items():
        if summary[name] != value:
            raise ValueError(f"full summary mismatch in {name}")

    positive = tuple(int(row) for row in lower_bound_18_graph())
    prior.validate_graph(positive)
    positive_k7 = sum(1 for _ in prior.clique_masks(positive, 7))
    expected_positive = {
        "vertices": 18,
        "K7_seeds": positive_k7,
        "status": "PASS_NOT_APPLICABLE_NO_K7",
    }
    if positive_k7 or report["positive_18_control"] != expected_positive:
        raise ValueError("known realizable 18-point positive control failed")

    return {
        "schema": 1,
        "kind": "d6_k7_ranktwo_pentad_full_verification",
        "status": "PASS",
        "report": {"path": str(report_path), "sha256": observed_hash},
        "symbolic_checks": symbolic,
        "checked": {
            "graphs": len(cover_counts),
            "covers": len(records),
            "pentads": pentads_checked,
            "certificates": certificates_checked,
            "marginal_graph_rejections": len(marginal),
            "positive_18_control": True,
        },
        "claim": (
            "All target quantifiers, 21 canonical pentads per cover, stored "
            "positive identities, graph-level marginal logic, pentad labeling "
            "symmetry, and the realizable 18-point control were reconstructed "
            "independently with exact arithmetic."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_ranktwo_pentad_pilot_report.json"),
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_ranktwo_pentad_pilot_verification.json"),
    )
    args = parser.parse_args()
    result = verify_report(args.report, args.report_sha256)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "status": result["status"],
        "checked": result["checked"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
