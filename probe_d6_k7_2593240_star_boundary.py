#!/usr/bin/env python3
"""Exact exploratory rank-one Schur contradiction for K7 graph 2593240.

The independently checked K7 support boundary leaves one ``zmask=0`` family
for each of the graph's two K7 seeds.  A required K6 in the normalized Gram
matrix leaves a positive-semidefinite Schur complement of rank at most one.
This probe reconstructs the exact zero/one entries (zeros only from disjoint
propagated support supersets), checks that the two systems are isomorphic, and
performs exact polynomial elimination on the rank-one tetrads.

The tetrads have one positive solution for the six inverse diagonal weights,
namely ``(2/3, 5, 1, 1, 1, 1)`` up to the displayed isomorphism.  At that
point the known Schur off-diagonal support is a nonzero star.  Such a support
cannot be ``h h^T``: two nonzero centre-leaf products force the corresponding
leaf-leaf product to be nonzero.  No floating-point arithmetic enters the
decision.

This is an exploratory exact certificate, not a production residue update.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 2_593_240
SEEDS = (
    (2, 4, 7, 10, 13, 15, 18),
    (6, 7, 8, 9, 12, 15, 18),
)
SEED_MASKS = tuple(sum(1 << vertex for vertex in seed) for seed in SEEDS)
CLIQUES = (
    (0, 1, 2, 7, 9, 10),
    (0, 1, 3, 7, 9, 10),
)
REMAINDERS = tuple(
    tuple(vertex for vertex in range(12) if vertex not in set(clique))
    for clique in CLIQUES
)
EXPECTED_PROPAGATED_MASKS = (
    (100, 127, 5, 98, 25, 19, 26, 108, 11, 20, 38, 72),
    (127, 98, 21, 10, 28, 97, 25, 102, 13, 18, 35, 68),
)
EXPECTED_CURRENT_COVERS = ((0, 8), (0, 32))
EXPECTED_RAW_COVER_COUNTS = (502, 502)
EXPECTED_RAW_COVERS_SHA256 = (
    "fb44017cc67555b5954484a5bd012ba2c8b4ff80e5ad878eb517a14ac791d769",
    "f915ae089f9885b1d800c698306a8180f57ebe9d8818481094ca07996930a6bd",
)
EXPECTED_BASIS_MASKS = (
    (43, 30, 54, 58, 46, 11),
    (57, 29, 43, 53, 45, 11),
)
PAIR_ORDER = tuple(combinations(range(6), 2))
EXPECTED_PAIR_TARGETS = (
    tuple(map(int, "011111111110111")),
    tuple(map(int, "111110111111101")),
)
# Maps canonical first-system labels to second-system labels.
REMAINDER_ISOMORPHISM = (2, 1, 3, 0, 4, 5)
COORDINATE_ISOMORPHISM = (1, 0, 2, 3, 4, 5)
EXPECTED_ADJACENCY_SHA256 = (
    "56f6a6a8bacf5a67e8fd0a0688b86b610d2d8c0842060b6423f8d484f340ed82"
)
UPSTREAM = {
    "d6_current_residue_manifest_v6.json": (
        "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
    ),
    "d6_current_residue_manifest_v6_verification.json": (
        "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
    ),
    "d6_k7_one_two_star_increment_report.json": (
        "12b3d18b1ea81961f58d831d4c7c322fbceb2300e7533b64161e8011fbe9d1ec"
    ),
    "d6_k7_one_two_star_increment_verification.json": (
        "c6f40578685036cb1156cb0ea8bf06d6004a70d5f376aff5ac1e16a1c0ca09eb"
    ),
}
EXPECTED_SEMANTICS = {
    "allowed_unpinned_coordinates_may_be_zero": True,
    "candidate_nonedges_optional": True,
    "cap500000_campaign_used": False,
    "floating_point_enters_rejection": False,
    "graph_rejected_if_any_k7_seed_is_infeasible": True,
    "normalization_factor_t_nonzero_on_cover_complement_N": True,
    "only_required_edges_enter_star_equations": True,
    "positive_sqrt7_embedding_checked_exactly": True,
    "propagated_masks_are_support_supersets": True,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


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


def bits(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(mask.bit_length()) if mask & (1 << index))


def validate_adjacency(adjacency: Sequence[int]) -> None:
    n = len(adjacency)
    full = (1 << n) - 1
    require(n == 19, "target graph order")
    for vertex, row in enumerate(adjacency):
        require(not row & ~full, "adjacency bit outside graph")
        require(not row & (1 << vertex), "adjacency loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(adjacency[other] & (1 << vertex)),
                "asymmetric adjacency",
            )


def is_clique(adjacency: Sequence[int], vertices: Sequence[int]) -> bool:
    return all(
        bool(adjacency[first] & (1 << second))
        for first, second in combinations(vertices, 2)
    )


def induced_graph(
    adjacency: Sequence[int], vertices: Sequence[int]
) -> tuple[int, ...]:
    position = {vertex: index for index, vertex in enumerate(vertices)}
    return tuple(
        sum(
            1 << position[other]
            for other in vertices
            if adjacency[vertex] & (1 << other)
        )
        for vertex in vertices
    )


def reconstruct_seed(
    adjacency: Sequence[int], seed: Sequence[int]
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...], int]:
    """Independently reconstruct outside labels, defects, L, and eligibility."""

    require(is_clique(adjacency, seed), "frozen seed is not a required K7")
    seed_position = {vertex: index for index, vertex in enumerate(seed)}
    outside = tuple(vertex for vertex in range(19) if vertex not in set(seed))
    defects = tuple(
        sum(
            1 << seed_position[seed_vertex]
            for seed_vertex in seed
            if not adjacency[vertex] & (1 << seed_vertex)
        )
        for vertex in outside
    )
    ladj = [0] * len(outside)
    for first, second in combinations(range(len(outside)), 2):
        if (
            adjacency[outside[first]] & (1 << outside[second])
            and not defects[first] & defects[second]
        ):
            ladj[first] |= 1 << second
            ladj[second] |= 1 << first
    eligible = sum(
        1 << index
        for index, defect in enumerate(defects)
        if defect.bit_count() >= 3
    )
    return outside, defects, tuple(ladj), eligible


def is_vertex_cover(ladj: Sequence[int], zmask: int) -> bool:
    remaining = ((1 << len(ladj)) - 1) & ~zmask
    return all(not (ladj[vertex] & remaining) for vertex in bits(remaining))


def enumerate_raw_covers(ladj: Sequence[int], eligible: int) -> tuple[int, ...]:
    return tuple(
        zmask
        for zmask in range(1 << len(ladj))
        if not zmask & ~eligible
        and zmask.bit_count() <= 7
        and is_vertex_cover(ladj, zmask)
    )


def load_boundary() -> dict:
    for name, expected in UPSTREAM.items():
        require(sha256(ROOT / name) == expected, f"upstream hash: {name}")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest_v6.json").read_text(encoding="utf-8")
    )
    manifest_check = json.loads(
        (ROOT / "d6_current_residue_manifest_v6_verification.json").read_text(
            encoding="utf-8"
        )
    )
    star = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_report.json").read_text(
            encoding="utf-8"
        )
    )
    star_check = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v6"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "v6 manifest boundary",
    )
    require(
        manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v6.json"],
        "v6 manifest verification boundary",
    )
    require(
        star.get("kind")
        == "d6_k7_one_free_neighbour_two_free_center_increment"
        and star.get("status") == "COMPLETE"
        and star.get("semantics") == EXPECTED_SEMANTICS,
        "star report boundary",
    )
    require(
        star_check.get("kind") == "d6_k7_one_two_star_increment_verification"
        and star_check.get("status") == "PASS"
        and star_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_one_two_star_increment_report.json"]
        and star_check.get("semantics") == EXPECTED_SEMANTICS
        and star_check.get("checked", {}).get("eligible_covers") == 19_932
        and star_check.get("checked", {}).get("current_passing_families") == 88,
        "independent star verification boundary",
    )
    graph = next(
        (
            record
            for record in manifest["classes"]["K7"]["graphs"]
            if int(record["index"]) == TARGET_INDEX
        ),
        None,
    )
    require(graph is not None, "target absent from v6 K7 class")
    adjacency = tuple(map(int, graph["adjacency"]))
    validate_adjacency(adjacency)
    require(stable_hash(list(adjacency)) == EXPECTED_ADJACENCY_SHA256, "adjacency hash")
    graph_record = next(
        (record for record in star["records"] if int(record["index"]) == TARGET_INDEX),
        None,
    )
    require(
        graph_record is not None and graph_record["decision"] == "SURVIVOR",
        "star target",
    )
    return {
        "adjacency": adjacency,
        "graph_record": graph_record,
        "independent_checked": star_check["checked"],
    }


def quantifier_audit(boundary: dict) -> tuple[dict, ...]:
    adjacency = boundary["adjacency"]
    k7s = tuple(
        chosen
        for chosen in combinations(range(19), 7)
        if is_clique(adjacency, chosen)
    )
    require(k7s == SEEDS, "unexpected K7 seed orbit/list")
    archived_by_seed = {
        int(record["seed_mask"]): record for record in boundary["graph_record"]["seeds"]
    }
    audits = []
    for position, seed in enumerate(SEEDS):
        outside, defects, ladj, eligible = reconstruct_seed(adjacency, seed)
        raw_covers = enumerate_raw_covers(ladj, eligible)
        require(
            len(raw_covers) == EXPECTED_RAW_COVER_COUNTS[position],
            f"raw eligible-cover count at seed {position}",
        )
        require(
            stable_hash(list(raw_covers)) == EXPECTED_RAW_COVERS_SHA256[position],
            f"raw cover hash at seed {position}",
        )
        archived = archived_by_seed[SEED_MASKS[position]]["current_covers"]
        expected_current = EXPECTED_CURRENT_COVERS[position]
        require(
            tuple(int(row["zmask"]) for row in archived) == expected_current,
            f"current covers at seed {position}",
        )
        by_z = {int(row["zmask"]): row for row in archived}
        nonzero = by_z[expected_current[1]]
        require(
            nonzero["status"] == "INFEASIBLE"
            and int(nonzero["current_passing_families"]) == 0,
            f"nonzero current cover at seed {position}",
        )
        zero = by_z[0]
        require(
            zero["status"] == "PASSING"
            and int(zero["current_passing_families"]) == 1
            and int(zero["star_passing_families"]) == 1,
            f"zmask-zero family multiplicity at seed {position}",
        )
        witness = zero["first_passing_witness"]
        require(witness.get("z_supports") == [], "zmask-zero support assignment")
        propagated = tuple(map(int, witness["propagated_masks"]))
        require(
            propagated == EXPECTED_PROPAGATED_MASKS[position],
            f"propagated mask boundary at seed {position}",
        )
        audits.append(
            {
                "seed": seed,
                "outside": outside,
                "defects": defects,
                "ladj": ladj,
                "eligible_mask": eligible,
                "raw_covers": raw_covers,
                "current_covers": expected_current,
                "prior_layer_eliminated_raw_covers": len(raw_covers) - len(archived),
                "current_infeasible_covers": (expected_current[1],),
                "sole_new_branch": 0,
                "propagated_masks": propagated,
            }
        )
    return tuple(audits)


def exact_k_entry(
    graph_n: Sequence[int], propagated: Sequence[int], first: int, second: int
) -> tuple[int, str]:
    """Return a proved K off-diagonal entry, never guessing a nonedge."""

    if graph_n[first] & (1 << second):
        require(
            propagated[first] & propagated[second],
            "required edge disjoint after propagation",
        )
        return 1, "required_unit_edge"
    if not propagated[first] & propagated[second]:
        return 0, "disjoint_propagated_support_supersets"
    raise ValueError(f"optional K entry remains at pair {first},{second}")


def reconstruct_exact_patterns(
    boundary: dict, audits: Sequence[dict]
) -> tuple[dict, ...]:
    patterns = []
    for position, audit in enumerate(audits):
        graph_n = induced_graph(boundary["adjacency"], audit["outside"])
        propagated = audit["propagated_masks"]
        clique = CLIQUES[position]
        remainder = REMAINDERS[position]
        require(is_clique(graph_n, clique), "local rank-six basis is not a clique")

        basis_masks = []
        basis_reasons = []
        for vertex in remainder:
            mask = 0
            reasons = []
            for coordinate, basis_vertex in enumerate(clique):
                target, reason = exact_k_entry(
                    graph_n, propagated, vertex, basis_vertex
                )
                mask |= target << coordinate
                reasons.append(reason)
            basis_masks.append(mask)
            basis_reasons.append(reasons)
        require(
            tuple(basis_masks) == EXPECTED_BASIS_MASKS[position],
            f"basis masks at seed {position}",
        )

        pair_targets = []
        pair_reasons = []
        for first, second in PAIR_ORDER:
            target, reason = exact_k_entry(
                graph_n, propagated, remainder[first], remainder[second]
            )
            pair_targets.append(target)
            pair_reasons.append(reason)
        require(
            tuple(pair_targets) == EXPECTED_PAIR_TARGETS[position],
            f"pair targets at seed {position}",
        )
        patterns.append(
            {
                "graph_n": graph_n,
                "basis_masks": tuple(basis_masks),
                "basis_entry_reasons": tuple(tuple(row) for row in basis_reasons),
                "pair_targets": tuple(pair_targets),
                "pair_entry_reasons": tuple(pair_reasons),
            }
        )
    return tuple(patterns)


def mapped_mask(mask: int, coordinate_map: Sequence[int]) -> int:
    return sum(
        1 << coordinate_map[coordinate]
        for coordinate in bits(mask)
    )


def verify_system_isomorphism(patterns: Sequence[dict]) -> None:
    first, second = patterns
    for old_remainder, new_remainder in enumerate(REMAINDER_ISOMORPHISM):
        require(
            mapped_mask(
                first["basis_masks"][old_remainder], COORDINATE_ISOMORPHISM
            )
            == second["basis_masks"][new_remainder],
            "basis-mask isomorphism",
        )
    first_targets = dict(zip(PAIR_ORDER, first["pair_targets"], strict=True))
    second_targets = dict(zip(PAIR_ORDER, second["pair_targets"], strict=True))
    for (old_first, old_second), target in first_targets.items():
        new_pair = tuple(
            sorted(
                (
                    REMAINDER_ISOMORPHISM[old_first],
                    REMAINDER_ISOMORPHISM[old_second],
                )
            )
        )
        require(second_targets[new_pair] == target, "remainder-pair isomorphism")


def mask_sum(mask: int, variables: Sequence[sp.Symbol]) -> sp.Expr:
    return sum((variables[index] for index in bits(mask)), sp.Integer(0))


def schur_pair_numerators(
    basis_masks: Sequence[int], pair_targets: Sequence[int]
) -> tuple[tuple[sp.Symbol, ...], dict[tuple[int, int], sp.Expr]]:
    """Return ``g_ij=T*S_ij`` for the known rank-one Schur entries."""

    variables = sp.symbols("a b c d e f")
    total = 1 + sum(variables)
    sums = [mask_sum(mask, variables) for mask in basis_masks]
    numerators = {}
    for pair, target in zip(PAIR_ORDER, pair_targets, strict=True):
        first, second = pair
        intersection = mask_sum(
            basis_masks[first] & basis_masks[second], variables
        )
        numerators[pair] = sp.expand(
            target * total - total * intersection + sums[first] * sums[second]
        )
    return variables, numerators


def tetrad_equations(
    numerators: dict[tuple[int, int], sp.Expr]
) -> dict[tuple[tuple[int, ...], int], sp.Expr]:
    equations = {}
    for vertices in combinations(range(6), 4):
        first, second, third, fourth = vertices
        common = numerators[(first, second)] * numerators[(third, fourth)]
        equations[(vertices, 0)] = sp.expand(
            common - numerators[(first, third)] * numerators[(second, fourth)]
        )
        equations[(vertices, 1)] = sp.expand(
            common - numerators[(first, fourth)] * numerators[(second, third)]
        )
    return equations


def require_polynomial_zero(
    expression: sp.Expr, variables: Iterable[sp.Symbol], label: str
) -> None:
    require(sp.Poly(sp.together(expression), *variables, domain=sp.QQ).is_zero, label)


def groebner_membership(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    targets: Sequence[sp.Expr],
    label: str,
) -> list[str]:
    basis = sp.groebner(equations, *variables, order="grevlex", domain=sp.QQ)
    for index, target in enumerate(targets):
        _quotients, remainder = basis.reduce(target)
        require(
            sp.Poly(remainder, *variables, domain=sp.QQ).is_zero,
            f"{label} target {index}",
        )
    return [str(sp.factor(polynomial.as_expr())) for polynomial in basis.polys]


def exact_algebra_certificate(pattern: dict) -> dict:
    variables, numerators = schur_pair_numerators(
        pattern["basis_masks"], pattern["pair_targets"]
    )
    a, b, c, d, e, f = variables
    total = 1 + sum(variables)
    tetrads = tetrad_equations(numerators)

    first_key = ((0, 3, 4, 5), 0)
    first_factor = -f * (a + 1) * (c - e) * total
    require_polynomial_zero(
        tetrads[first_key] - first_factor,
        variables,
        "first equality factor",
    )

    # Every variable and T are positive.  The first equality therefore gives
    # e=c.  Exact ideal membership supplies three small consequences after
    # that substitution.  They imply f=c and d=c in the positive orthant.
    stage_one_variables = (a, b, d, f, c)
    stage_one_equations = tuple(
        sp.cancel(expression.subs(e, c) / total.subs(e, c))
        for expression in tetrads.values()
    )
    stage_one_targets = (
        (c - f) * (2 * c + 1),
        (c - d) * (c - 1),
        (c - d) * (c + d - 2),
    )
    stage_one_basis = groebner_membership(
        stage_one_equations,
        stage_one_variables,
        stage_one_targets,
        "stage-one elimination",
    )

    # With d=e=f=c, a second exact elimination leaves one positive root and
    # then fixes a and b rationally.
    equal_tail = {d: c, e: c, f: c}
    stage_two_equations = tuple(
        sp.cancel(expression.subs(equal_tail) / total.subs(equal_tail))
        for expression in tetrads.values()
    )
    stage_two_targets = (
        (c - 1) * (c + 1) * (2 * c + 1) * (3 * c + 1),
        9 * a - 6 * c**3 - 11 * c**2 + 3 * c + 8,
        3 * b - 6 * c**3 - 11 * c**2 + 3 * c - 1,
    )
    stage_two_basis = groebner_membership(
        stage_two_equations,
        (a, b, c),
        stage_two_targets,
        "stage-two elimination",
    )

    unique_tetrad_point = {
        a: sp.Rational(2, 3),
        b: 5,
        c: 1,
        d: 1,
        e: 1,
        f: 1,
    }
    require(
        all(sp.expand(expression.subs(unique_tetrad_point)) == 0 for expression in tetrads.values()),
        "unique point does not satisfy all tetrads",
    )
    evaluated = {
        pair: sp.factor(expression.subs(unique_tetrad_point))
        for pair, expression in numerators.items()
    }
    expected_nonzero = {
        (0, 1): sp.Rational(-8, 3),
        (0, 2): sp.Integer(8),
        (0, 3): sp.Rational(-8, 3),
        (0, 4): sp.Rational(-8, 3),
        (0, 5): sp.Rational(-28, 3),
    }
    require(
        {pair: value for pair, value in evaluated.items() if value}
        == expected_nonzero,
        "star-support evaluation",
    )
    require(
        evaluated[(0, 1)] != 0
        and evaluated[(0, 2)] != 0
        and evaluated[(1, 2)] == 0,
        "rank-one support-clique contradiction",
    )

    return {
        "variables": [str(variable) for variable in variables],
        "pair_order": [list(pair) for pair in PAIR_ORDER],
        "pair_numerators": [str(numerators[pair]) for pair in PAIR_ORDER],
        "tetrads_checked": len(tetrads),
        "first_factor_identity": str(sp.factor(first_factor)),
        "stage_one_targets": [str(sp.factor(target)) for target in stage_one_targets],
        "stage_one_groebner_basis": stage_one_basis,
        "stage_one_deduction": (
            "e=c; (c-f)(2c+1)=0 and c>0 give f=c; if d!=c, "
            "the other two products give c=1 and c+d=2, hence d=c, "
            "a contradiction. Thus d=e=f=c."
        ),
        "stage_two_targets": [str(sp.factor(target)) for target in stage_two_targets],
        "stage_two_groebner_basis": stage_two_basis,
        "unique_positive_tetrad_point": ["2/3", "5", "1", "1", "1", "1"],
        "star_nonzero_entries": {
            f"{first},{second}": str(value)
            for (first, second), value in expected_nonzero.items()
        },
        "rank_one_contradiction": (
            "If g_ij=h_i*h_j, then g_01 and g_02 nonzero imply h_0,h_1,h_2 "
            "are all nonzero, so g_12=h_1*h_2 is nonzero; exactly g_12=0."
        ),
    }


def build_report() -> dict:
    boundary = load_boundary()
    audits = quantifier_audit(boundary)
    patterns = reconstruct_exact_patterns(boundary, audits)
    verify_system_isomorphism(patterns)
    algebra = exact_algebra_certificate(patterns[0])
    return {
        "schema": 1,
        "kind": "d6_k7_rankone_star_2593240_exact_probe",
        "status": "EXACT_CONTRADICTION",
        "claim": (
            "At the pinned independently verified K7 boundary, both required "
            "K7 seeds of graph 2593240 have no realizable cover/support branch."
        ),
        "theorem_credit": "EXPLORATORY_NOT_YET_PRODUCTION_UNION",
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "source_sha256": sha256(Path(__file__).resolve()),
        "semantics": {
            "candidate_nonedges_optional": True,
            "zero_K_entries_require_disjoint_propagated_support_supersets": True,
            "required_edges_give_unit_K_entries": True,
            "only_one_infeasible_required_K7_seed_needed": True,
            "both_required_K7_seeds_checked": True,
            "floating_point_enters_contradiction": False,
        },
        "target": {
            "index": TARGET_INDEX,
            "adjacency_sha256": EXPECTED_ADJACENCY_SHA256,
            "required_k7_seeds": [list(seed) for seed in SEEDS],
            "required_k7_seed_masks": list(SEED_MASKS),
        },
        "quantifier": [
            {
                "seed": list(audit["seed"]),
                "outside": list(audit["outside"]),
                "raw_eligible_covers": len(audit["raw_covers"]),
                "raw_eligible_covers_sha256": stable_hash(list(audit["raw_covers"])),
                "current_covers": list(audit["current_covers"]),
                "prior_layer_eliminated_raw_covers": audit[
                    "prior_layer_eliminated_raw_covers"
                ],
                "current_infeasible_covers": list(audit["current_infeasible_covers"]),
                "sole_new_branch": audit["sole_new_branch"],
                "sole_new_branch_passing_families": 1,
            }
            for audit in audits
        ],
        "exact_patterns": [
            {
                "local_rank6_clique": list(CLIQUES[position]),
                "global_rank6_clique": [
                    audits[position]["outside"][index]
                    for index in CLIQUES[position]
                ],
                "local_remainder": list(REMAINDERS[position]),
                "global_remainder": [
                    audits[position]["outside"][index]
                    for index in REMAINDERS[position]
                ],
                "propagated_support_supersets": list(
                    audits[position]["propagated_masks"]
                ),
                "basis_masks": list(pattern["basis_masks"]),
                "pair_order": [list(pair) for pair in PAIR_ORDER],
                "pair_targets": list(pattern["pair_targets"]),
                "pair_target_string": "".join(map(str, pattern["pair_targets"])),
                "zero_entry_justification": (
                    "Every zero in the displayed basis/pair pattern has disjoint "
                    "propagated support supersets; it is not inferred from a "
                    "candidate nonedge alone."
                ),
            }
            for position, pattern in enumerate(patterns)
        ],
        "system_isomorphism": {
            "first_remainder_to_second": list(REMAINDER_ISOMORPHISM),
            "first_coordinate_to_second": list(COORDINATE_ISOMORPHISM),
            "checked_exactly": True,
        },
        "algebra": algebra,
        "rank_argument": {
            "normalized_gram_rank_upper_bound": 7,
            "required_clique_order": 6,
            "clique_principal_block_positive_definite": True,
            "schur_complement_rank_upper_bound": 1,
            "schur_complement_positive_semidefinite": True,
            "strict_domain": "a,b,c,d,e,f > 0",
        },
    }


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/tmp/d6_k7_2593240_star_boundary_probe.json"),
    )
    args = parser.parse_args()
    report = build_report()
    atomic_json(args.output.resolve(), report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "index": TARGET_INDEX,
                "required_k7_seeds": len(SEEDS),
                "raw_eligible_covers": sum(
                    item["raw_eligible_covers"] for item in report["quantifier"]
                ),
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
