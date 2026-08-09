#!/usr/bin/env python3
"""Exact exploratory Schur contradiction for K7 graph 3949382.

This standalone probe pins the independently checked v6/star boundary,
reconstructs one required K7 seed and its full zero-factor quantifier, and
checks a strict-positivity contradiction for the sole surviving support
family.  Candidate nonedges are never assumed non-unit: a zero Schur entry is
used only when the independently propagated support *supersets* are disjoint.

The output is an exploratory exact certificate, not a production residue
update.  No floating-point arithmetic enters the decision.
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
TARGET_INDEX = 3_949_382
SEED = (2, 5, 7, 12, 14, 16, 18)
SEED_MASK = sum(1 << vertex for vertex in SEED)
CLIQUE = (0, 2, 6, 7, 9, 10, 11)
REMAINDER = (1, 3, 4, 5, 8)
EXPECTED_BASIS_MASKS = (83, 47, 98, 86, 90)
PAIR_ORDER = tuple((i, j) for j in range(1, 5) for i in range(j))
EXPECTED_PAIR_TARGETS = tuple(map(int, "0111011011"))
EXPECTED_PROPAGATED_MASKS = (
    84, 36, 127, 66, 34, 33, 81, 88, 40, 29, 82, 61,
)
EXPECTED_CURRENT_COVERS = (0, 2048, 3072)
EXPECTED_RAW_COVERS_SHA256 = (
    "5a47d6b62e4b141b1f024d577c3aaec0533cff7480b1ca45de921a55af43d879"
)
EXPECTED_ADJACENCY_SHA256 = (
    "ede388747bbaa2de587db5986806c5805e7961ab278a8f5162a7d656863e36b8"
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
    adjacency: Sequence[int],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...], int]:
    """Independently reconstruct outside labels, defects, L, and eligibility."""

    require(is_clique(adjacency, SEED), "frozen seed is not a required K7")
    seed_position = {vertex: index for index, vertex in enumerate(SEED)}
    outside = tuple(vertex for vertex in range(19) if vertex not in set(SEED))
    defects = tuple(
        sum(
            1 << seed_position[seed_vertex]
            for seed_vertex in SEED
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
    require(graph_record is not None and graph_record["decision"] == "SURVIVOR", "star target")
    seed_record = next(
        (
            record
            for record in graph_record["seeds"]
            if int(record["seed_mask"]) == SEED_MASK
        ),
        None,
    )
    require(seed_record is not None and tuple(seed_record["seed"]) == SEED, "seed record")
    return {
        "adjacency": adjacency,
        "graph_record": graph_record,
        "seed_record": seed_record,
        "independent_checked": star_check["checked"],
    }


def quantifier_audit(boundary: dict) -> dict:
    adjacency = boundary["adjacency"]
    outside, defects, ladj, eligible = reconstruct_seed(adjacency)
    raw_covers = enumerate_raw_covers(ladj, eligible)
    require(len(raw_covers) == 128, "raw eligible-cover count")
    require(stable_hash(list(raw_covers)) == EXPECTED_RAW_COVERS_SHA256, "raw cover hash")
    archived = boundary["seed_record"]["current_covers"]
    require(tuple(int(row["zmask"]) for row in archived) == EXPECTED_CURRENT_COVERS, "current covers")
    by_z = {int(row["zmask"]): row for row in archived}
    require(
        by_z[2048]["status"] == "INFEASIBLE"
        and int(by_z[2048]["current_passing_families"]) == 0
        and by_z[3072]["status"] == "INFEASIBLE"
        and int(by_z[3072]["current_passing_families"]) == 0,
        "nonzero current covers were not eliminated",
    )
    zero = by_z[0]
    require(
        zero["status"] == "PASSING"
        and int(zero["current_passing_families"]) == 1
        and int(zero["star_passing_families"]) == 1,
        "zmask-zero family multiplicity",
    )
    witness = zero["first_passing_witness"]
    require(witness.get("z_supports") == [], "zmask-zero support assignment")
    propagated = tuple(map(int, witness["propagated_masks"]))
    require(propagated == EXPECTED_PROPAGATED_MASKS, "propagated mask boundary")
    return {
        "outside": outside,
        "defects": defects,
        "ladj": ladj,
        "eligible_mask": eligible,
        "raw_covers": raw_covers,
        "current_covers": EXPECTED_CURRENT_COVERS,
        "prior_layer_eliminated_raw_covers": len(raw_covers) - len(archived),
        "current_infeasible_covers": (2048, 3072),
        "sole_new_branch": 0,
        "propagated_masks": propagated,
    }


def exact_k_entry(
    graph_n: Sequence[int], propagated: Sequence[int], first: int, second: int
) -> tuple[int, str]:
    """Return a proved K off-diagonal entry, never guessing a nonedge."""

    if graph_n[first] & (1 << second):
        require(propagated[first] & propagated[second], "required edge disjoint after propagation")
        return 1, "required_unit_edge"
    if not propagated[first] & propagated[second]:
        return 0, "disjoint_propagated_support_supersets"
    raise ValueError(f"optional K entry remains at pair {first},{second}")


def reconstruct_exact_pattern(boundary: dict, audit: dict) -> dict:
    outside = audit["outside"]
    graph_n = induced_graph(boundary["adjacency"], outside)
    propagated = audit["propagated_masks"]
    require(is_clique(graph_n, CLIQUE), "local rank-seven basis is not a clique")
    seven_cliques = tuple(
        chosen
        for chosen in combinations(range(12), 7)
        if is_clique(graph_n, chosen)
    )
    require(seven_cliques == (CLIQUE,), "rank-seven clique is not unique")

    basis_masks = []
    basis_reasons = []
    for vertex in REMAINDER:
        mask = 0
        reasons = []
        for coordinate, basis_vertex in enumerate(CLIQUE):
            target, reason = exact_k_entry(
                graph_n, propagated, vertex, basis_vertex
            )
            mask |= target << coordinate
            reasons.append(reason)
        basis_masks.append(mask)
        basis_reasons.append(reasons)
    require(tuple(basis_masks) == EXPECTED_BASIS_MASKS, "basis masks")

    pair_targets = []
    pair_reasons = []
    for first, second in PAIR_ORDER:
        target, reason = exact_k_entry(
            graph_n, propagated, REMAINDER[first], REMAINDER[second]
        )
        pair_targets.append(target)
        pair_reasons.append(reason)
    require(tuple(pair_targets) == EXPECTED_PAIR_TARGETS, "remainder pair targets")
    return {
        "graph_n": graph_n,
        "basis_masks": tuple(basis_masks),
        "basis_entry_reasons": tuple(tuple(row) for row in basis_reasons),
        "pair_targets": tuple(pair_targets),
        "pair_entry_reasons": tuple(pair_reasons),
    }


def mask_sum(mask: int, variables: Sequence[sp.Symbol]) -> sp.Expr:
    return sum((variables[index] for index in bits(mask)), sp.Integer(0))


def schur_pair_equations(
    basis_masks: Sequence[int], pair_targets: Sequence[int]
) -> tuple[tuple[sp.Symbol, ...], dict[tuple[int, int], sp.Expr]]:
    variables = sp.symbols("a b c d e f g")
    total = 1 + sum(variables)
    sums = [mask_sum(mask, variables) for mask in basis_masks]
    equations = {}
    for (first, second), target in zip(PAIR_ORDER, pair_targets, strict=True):
        intersection = mask_sum(
            basis_masks[first] & basis_masks[second], variables
        )
        equations[(first, second)] = sp.expand(
            total * intersection - sums[first] * sums[second] - target * total
        )
    return variables, equations


def require_polynomial_zero(expression: sp.Expr, variables: Iterable[sp.Symbol], label: str) -> None:
    require(sp.Poly(sp.together(expression), *variables, domain=sp.QQ).is_zero, label)


def exact_algebra_certificate(pattern: dict) -> dict:
    variables, equations = schur_pair_equations(
        pattern["basis_masks"], pattern["pair_targets"]
    )
    a, b, c, d, e, f, g = variables
    rsum = b + f + g
    require_polynomial_zero(
        equations[(0, 2)] - equations[(2, 3)] + (a - c) * rsum,
        variables,
        "a=c difference identity",
    )
    require_polynomial_zero(
        equations[(2, 3)] - equations[(2, 4)] + (c - d) * rsum,
        variables,
        "c=d difference identity",
    )

    t = sp.symbols("t")
    substitution = {a: t, c: t, d: t}
    total = 1 + 3 * t + b + e + f + g
    p = t + b + e + g
    q = 3 * t + b + f
    r = b + f + g
    reduced = (
        sp.expand((t + b) * total - p * q),
        sp.expand((p - t - 1) * total - p**2),
        sp.expand((b + g - 1) * total - p * r),
        sp.expand((b + f - 1) * total - q * r),
    )
    mapping = (0, 2, 3, 1, 0, 2, 1, 0, 2, 1)
    for pair, expected in zip(PAIR_ORDER, mapping, strict=True):
        require_polynomial_zero(
            equations[pair].subs(substitution) - reduced[expected],
            (t, b, e, f, g),
            f"ten-to-four reduction {pair}",
        )

    h = t + b
    ell = 2 * t + f
    u = b + g - 1
    require_polynomial_zero(total - (p + ell + 1), (t, b, e, f, g), "T decomposition")
    require_polynomial_zero(q - (h + ell), (t, b, f), "Q decomposition")
    require_polynomial_zero(
        reduced[0] - (h * (ell + 1) - p * ell),
        (t, b, e, f, g),
        "E0 elimination identity",
    )

    p_solution = sp.cancel(h * (ell + 1) / ell)
    g_solution = sp.cancel(1 - b + h * (f + 1) / ell)
    e_solution = sp.cancel(h / ell - g_solution)
    e_from_p = sp.solve(sp.Eq(p, p_solution), e, dict=False)[0]
    require_polynomial_zero(
        e_from_p.subs(g, g_solution) - e_solution,
        (t, b, f),
        "P and E2 support reconstruction",
    )
    e2_after_p = sp.factor(reduced[2].subs(e, e_from_p))
    require_polynomial_zero(
        e2_after_p
        - (ell + 1) * (u * ell - h * (f + 1)) / ell,
        (t, b, f, g),
        "E2 elimination identity",
    )

    e1_after = sp.factor(
        reduced[1].subs({e: e_solution, g: g_solution}, simultaneous=True)
    )
    a_polynomial = t**2 + 2 * t + f - b * (t + f)
    require_polynomial_zero(
        e1_after + (ell + 1) * a_polynomial / ell,
        (t, b, f),
        "E1 factor identity",
    )
    b_solution = sp.cancel((t**2 + 2 * t + f) / (t + f))
    e_final = sp.factor(e_solution.subs(b, b_solution))
    g_final = sp.factor(g_solution.subs(b, b_solution))
    require_polynomial_zero(
        e_final - (t + 1) * (t - f) / (t + f),
        (t, f),
        "final e formula",
    )
    require_polynomial_zero(
        g_final + (t + 1) * (t - f - 1) / (t + f),
        (t, f),
        "final g formula",
    )

    e3_after = sp.factor(
        reduced[3].subs(
            {b: b_solution, e: e_final, g: g_final}, simultaneous=True
        )
    )
    final_polynomial = t**2 - t * f - t - f
    require_polynomial_zero(
        e3_after - (ell + 1) ** 2 * final_polynomial / (t + f) ** 2,
        (t, f),
        "E3 final factor identity",
    )

    boundary = {a: 1, b: 3, c: 1, d: 1, e: 2, f: 0, g: 0}
    require(
        all(sp.expand(equation.subs(boundary)) == 0 for equation in equations.values()),
        "nonnegative boundary control",
    )
    return {
        "variables": [str(variable) for variable in variables],
        "pair_order": [list(pair) for pair in PAIR_ORDER],
        "pair_equations": [str(equations[pair]) for pair in PAIR_ORDER],
        "equal_sum_deductions": [
            "E_AC-E_CD=-(a-c)(b+f+g)=0 and b+f+g>0, hence a=c",
            "E_CD-E_CE=-(c-d)(b+f+g)=0 and b+f+g>0, hence c=d",
        ],
        "reduced_variables": ["t", "b", "e", "f", "g"],
        "reduced_equations": [str(expression) for expression in reduced],
        "positive_factors": ["t", "f", "2*t+f", "t+f", "2*t+f+1"],
        "solved_forms": {
            "b": str(b_solution),
            "e": str(e_final),
            "g": str(g_final),
        },
        "final_factor": str(final_polynomial),
        "strict_contradiction": (
            "E3=0 gives t*(t-f-1)=f>0, hence t-f-1>0; "
            "but g=-(t+1)*(t-f-1)/(t+f)<0, contradicting g>0."
        ),
        "boundary_control": {
            "point": [1, 3, 1, 1, 2, 0, 0],
            "all_ten_equations_zero": True,
            "strict_domain": False,
        },
    }


def build_report() -> dict:
    boundary = load_boundary()
    audit = quantifier_audit(boundary)
    pattern = reconstruct_exact_pattern(boundary, audit)
    algebra = exact_algebra_certificate(pattern)
    zero_pairs = [
        list(pair)
        for pair, target in zip(PAIR_ORDER, pattern["pair_targets"], strict=True)
        if target == 0
    ]
    return {
        "schema": 1,
        "kind": "d6_k7_saturating_clique_3949382_exact_probe",
        "status": "EXACT_CONTRADICTION",
        "claim": (
            "At the pinned independently verified K7 boundary, one required "
            "K7 seed of graph 3949382 has no realizable cover/support branch."
        ),
        "theorem_credit": "EXPLORATORY_NOT_YET_PRODUCTION_UNION",
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "source_sha256": sha256(Path(__file__).resolve()),
        "semantics": {
            "candidate_nonedges_optional": True,
            "zero_K_entries_require_disjoint_propagated_support_supersets": True,
            "required_edges_give_unit_K_entries": True,
            "only_one_infeasible_required_K7_seed_needed": True,
            "floating_point_enters_contradiction": False,
        },
        "target": {
            "index": TARGET_INDEX,
            "adjacency_sha256": EXPECTED_ADJACENCY_SHA256,
            "seed": list(SEED),
            "seed_mask": SEED_MASK,
            "outside": list(audit["outside"]),
        },
        "quantifier": {
            "raw_eligible_covers": len(audit["raw_covers"]),
            "raw_eligible_covers_sha256": stable_hash(list(audit["raw_covers"])),
            "current_covers": list(audit["current_covers"]),
            "prior_layer_eliminated_raw_covers": audit[
                "prior_layer_eliminated_raw_covers"
            ],
            "current_infeasible_covers": list(audit["current_infeasible_covers"]),
            "sole_new_branch": audit["sole_new_branch"],
            "sole_new_branch_passing_families": 1,
            "graph_rejection_quantifier": (
                "Every realization must realize this required K7 seed. Its "
                "actual Z/support branch is exhaustive above; all branches fail."
            ),
        },
        "exact_pattern": {
            "local_rank7_clique": list(CLIQUE),
            "global_rank7_clique": [audit["outside"][index] for index in CLIQUE],
            "local_remainder": list(REMAINDER),
            "global_remainder": [audit["outside"][index] for index in REMAINDER],
            "propagated_support_supersets": list(audit["propagated_masks"]),
            "basis_masks": list(pattern["basis_masks"]),
            "pair_order": [list(pair) for pair in PAIR_ORDER],
            "pair_targets": list(pattern["pair_targets"]),
            "pair_target_string": "".join(map(str, pattern["pair_targets"])),
            "zero_remainder_pairs": zero_pairs,
            "zero_entry_justification": (
                "Every zero in the displayed basis/pair pattern has disjoint "
                "propagated support supersets; it is not inferred from a "
                "candidate nonedge alone."
            ),
        },
        "algebra": algebra,
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
        default=Path("/tmp/d6_k7_schur_3949382_probe.json"),
    )
    args = parser.parse_args()
    report = build_report()
    atomic_json(args.output.resolve(), report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "index": TARGET_INDEX,
                "raw_eligible_covers": report["quantifier"]["raw_eligible_covers"],
                "current_covers": report["quantifier"]["current_covers"],
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
