#!/usr/bin/env python3
"""Independent checker for the exact graph-3949382 K7 Schur increment.

This checker imports neither the production builder nor the exploratory probe.
It reconstructs the graph/seed/cover boundary and every exact Schur entry from
the pinned artifacts, then verifies a different polynomial elimination from
the one archived by the producer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
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
PAIR_ORDER = tuple((first, second) for second in range(1, 5) for first in range(second))
EXPECTED_BASIS_MASKS = (83, 47, 98, 86, 90)
EXPECTED_PAIR_TARGETS = tuple(map(int, "0111011011"))
EXPECTED_PROPAGATED = (84, 36, 127, 66, 34, 33, 81, 88, 40, 29, 82, 61)
EXPECTED_COVERS = (0, 2048, 3072)
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
FROZEN_EXPLORATORY = {
    "d6_k7_schur_3949382.md": (
        "660de8bb17947dc7848886b51e11327f0cc4645e55dc993ae5df423a8de0aaf4"
    ),
    "probe_d6_k7_schur_3949382.py": (
        "96750e6591aae18dd25e59d382ed64c316575431942780b002be12da4594cb86"
    ),
    "test_probe_d6_k7_schur_3949382.py": (
        "ab3842c4f55e1fa6216ad10e18a197b99dc4c470e9c0e00cba7c26e6eab911ab"
    ),
}
SOURCE_FILES = {
    *FROZEN_EXPLORATORY,
    "build_d6_k7_schur_3949382_increment.py",
    "verify_d6_k7_schur_3949382_increment.py",
    "test_d6_k7_schur_3949382_increment.py",
    "d6_k7_schur_3949382_increment.md",
}
EXPECTED_INPUT = (
    316173,
    2581209,
    2593240,
    3595554,
    3648882,
    3729907,
    3935560,
    3936177,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
    3949382,
)
EXPECTED_INPUT_SHA256 = (
    "b0235ce7e95918197f0f3b4a57d26c8c83632ff51fd678a989f50197591eee8e"
)
EXPECTED_OUTPUT_SHA256 = (
    "715c52011421fb6d8721d52339ed3437ef1cec344ddf3f03a3b0b618f35b00d7"
)
EXPECTED_SEMANTICS = {
    "all_actual_support_branches_for_seed_are_quantified": True,
    "candidate_nonedges_optional": True,
    "floating_point_enters_rejection": False,
    "one_infeasible_required_K7_seed_rejects_graph": True,
    "propagated_masks_are_support_supersets": True,
    "required_edges_give_unit_K_entries": True,
    "strict_positive_diagonal_parameters_encode_distinct_points": True,
    "zero_K_entries_require_disjoint_propagated_support_supersets": True,
}
EXPECTED_STAR_SEMANTICS = {
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


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def bits(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(mask.bit_length()) if mask & (1 << index))


def validate_adjacency(adjacency: Sequence[int]) -> None:
    require(len(adjacency) == 19, "graph order")
    full = (1 << len(adjacency)) - 1
    for vertex, row in enumerate(adjacency):
        require(not row & ~full, "adjacency outside graph")
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


def induced_graph(adjacency: Sequence[int], vertices: Sequence[int]) -> tuple[int, ...]:
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
    require(is_clique(adjacency, SEED), "required K7 seed")
    seed_set = set(SEED)
    seed_position = {vertex: position for position, vertex in enumerate(SEED)}
    outside = tuple(vertex for vertex in range(19) if vertex not in seed_set)
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
        1 << index for index, defect in enumerate(defects) if defect.bit_count() >= 3
    )
    return outside, defects, tuple(ladj), eligible


def is_vertex_cover(ladj: Sequence[int], zmask: int) -> bool:
    remaining = ((1 << len(ladj)) - 1) & ~zmask
    return all(not ladj[vertex] & remaining for vertex in bits(remaining))


def raw_covers(ladj: Sequence[int], eligible: int) -> tuple[int, ...]:
    return tuple(
        zmask
        for zmask in range(1 << len(ladj))
        if not zmask & ~eligible
        and zmask.bit_count() <= 7
        and is_vertex_cover(ladj, zmask)
    )


def load_upstream() -> tuple[dict, dict, dict]:
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
        "manifest boundary",
    )
    require(
        manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v6.json"],
        "manifest checker boundary",
    )
    require(
        star.get("kind")
        == "d6_k7_one_free_neighbour_two_free_center_increment"
        and star.get("status") == "COMPLETE"
        and star.get("semantics") == EXPECTED_STAR_SEMANTICS,
        "star report boundary",
    )
    require(
        star_check.get("kind") == "d6_k7_one_two_star_increment_verification"
        and star_check.get("status") == "PASS"
        and star_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_one_two_star_increment_report.json"]
        and star_check.get("semantics") == EXPECTED_STAR_SEMANTICS
        and star_check.get("checked", {}).get("eligible_covers") == 19_932
        and star_check.get("checked", {}).get("current_passing_families") == 88,
        "star checker boundary",
    )
    return manifest, star, star_check


def validate_report_schema(report: dict) -> None:
    output = [index for index in EXPECTED_INPUT if index != TARGET_INDEX]
    require(report.get("schema") == 1, "report schema")
    require(report.get("kind") == "d6_k7_schur_3949382_exact_increment", "report kind")
    require(report.get("status") == "COMPLETE_EXACT_REJECTION", "report status")
    require(report.get("semantics") == EXPECTED_SEMANTICS, "report semantics")
    require(report.get("upstream_sha256") == dict(sorted(UPSTREAM.items())), "upstream roots")
    require(
        report.get("frozen_exploratory_sha256")
        == dict(sorted(FROZEN_EXPLORATORY.items())),
        "exploratory roots",
    )
    require(
        report.get("input", {}).get("ordered_indices") == list(EXPECTED_INPUT)
        and report.get("input", {}).get("ordered_indices_sha256")
        == EXPECTED_INPUT_SHA256,
        "input residue",
    )
    summary = report.get("summary", {})
    require(
        summary.get("graphs_rejected") == 1
        and summary.get("rejected_indices") == [TARGET_INDEX]
        and summary.get("rejected_indices_sha256") == stable_hash([TARGET_INDEX])
        and summary.get("graphs_surviving") == 15
        and summary.get("ordered_survivor_indices") == output
        and summary.get("ordered_survivor_indices_sha256") == EXPECTED_OUTPUT_SHA256,
        "union-ready summary",
    )
    require(
        report.get("certificate_sha256") == stable_hash(report.get("certificate")),
        "certificate hash",
    )


def validate_source_boundary(report: dict) -> dict:
    sources = report.get("source_sha256")
    require(isinstance(sources, dict) and set(sources) == SOURCE_FILES, "source set")
    for name, expected in FROZEN_EXPLORATORY.items():
        require(sources.get(name) == expected, f"frozen source in report: {name}")
    launch = report.get("execution", {}).get("git", {})
    require(
        launch.get("branch") == "codex/dimension6"
        and launch.get("proof_and_checker_sources_equal_committed_blobs") is True,
        "launch branch/source claim",
    )
    commit = launch.get("commit")
    require(isinstance(commit, str) and len(commit) == 40, "launch commit")
    for name, expected in sources.items():
        require(sha256(ROOT / name) == expected, f"working source hash: {name}")
        blob = subprocess.run(
            ["git", "show", f"{commit}:{name}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        require(hashlib.sha256(blob).hexdigest() == expected, f"committed source: {name}")
    porcelain = "\n".join(launch.get("porcelain_lines", ()))
    require(
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        == launch.get("porcelain_sha256"),
        "launch porcelain binding",
    )
    return {
        "commit": commit,
        "branch": launch["branch"],
        "source_files": len(sources),
        "proof_and_checker_sources_equal_committed_blobs": True,
    }


def independent_quantifier_and_pattern(manifest: dict, star: dict) -> dict:
    graph = next(
        record
        for record in manifest["classes"]["K7"]["graphs"]
        if int(record["index"]) == TARGET_INDEX
    )
    adjacency = tuple(map(int, graph["adjacency"]))
    validate_adjacency(adjacency)
    require(stable_hash(list(adjacency)) == EXPECTED_ADJACENCY_SHA256, "adjacency hash")
    outside, defects, ladj, eligible = reconstruct_seed(adjacency)
    covers = raw_covers(ladj, eligible)
    require(len(covers) == 128, "raw cover count")
    require(stable_hash(list(covers)) == EXPECTED_RAW_COVERS_SHA256, "raw cover hash")

    graph_record = next(
        record for record in star["records"] if int(record["index"]) == TARGET_INDEX
    )
    require(graph_record.get("decision") == "SURVIVOR", "upstream target status")
    seed_record = next(
        record for record in graph_record["seeds"] if int(record["seed_mask"]) == SEED_MASK
    )
    require(tuple(seed_record["seed"]) == SEED, "archived seed")
    archived = seed_record["current_covers"]
    require(tuple(int(row["zmask"]) for row in archived) == EXPECTED_COVERS, "current covers")
    by_z = {int(row["zmask"]): row for row in archived}
    for zmask in (2048, 3072):
        require(
            by_z[zmask]["status"] == "INFEASIBLE"
            and int(by_z[zmask]["current_passing_families"]) == 0,
            f"infeasible cover {zmask}",
        )
    zero = by_z[0]
    require(
        zero["status"] == "PASSING"
        and int(zero["current_passing_families"]) == 1
        and int(zero["star_passing_families"]) == 1,
        "sole current family",
    )
    witness = zero["first_passing_witness"]
    require(witness.get("z_supports") == [], "z=0 witness")
    propagated = tuple(map(int, witness["propagated_masks"]))
    require(propagated == EXPECTED_PROPAGATED, "propagated supports")

    graph_n = induced_graph(adjacency, outside)
    cliques = tuple(
        vertices
        for vertices in combinations(range(len(outside)), 7)
        if is_clique(graph_n, vertices)
    )
    require(cliques == (CLIQUE,), "unique required rank-seven clique")

    reasons: list[str] = []

    def exact_entry(first: int, second: int) -> int:
        if graph_n[first] & (1 << second):
            require(propagated[first] & propagated[second], "edge with disjoint supports")
            reasons.append("required_unit_edge")
            return 1
        require(
            not propagated[first] & propagated[second],
            f"optional unresolved entry {first},{second}",
        )
        reasons.append("disjoint_propagated_support_supersets")
        return 0

    basis_masks = []
    for vertex in REMAINDER:
        mask = 0
        for coordinate, basis_vertex in enumerate(CLIQUE):
            mask |= exact_entry(vertex, basis_vertex) << coordinate
        basis_masks.append(mask)
    pair_targets = tuple(
        exact_entry(REMAINDER[first], REMAINDER[second])
        for first, second in PAIR_ORDER
    )
    require(tuple(basis_masks) == EXPECTED_BASIS_MASKS, "basis masks")
    require(pair_targets == EXPECTED_PAIR_TARGETS, "pair targets")

    certificate = {
        "target": {
            "index": TARGET_INDEX,
            "adjacency_sha256": EXPECTED_ADJACENCY_SHA256,
            "seed": list(SEED),
            "seed_mask": SEED_MASK,
            "outside": list(outside),
        },
        "quantifier": {
            "raw_eligible_covers": len(covers),
            "raw_eligible_covers_sha256": stable_hash(list(covers)),
            "current_covers": list(EXPECTED_COVERS),
            "prior_layer_eliminated_raw_covers": 125,
            "current_infeasible_covers": [2048, 3072],
            "sole_new_branch": 0,
            "sole_new_branch_passing_families": 1,
        },
        "pattern": {
            "local_rank7_clique": list(CLIQUE),
            "global_rank7_clique": [outside[index] for index in CLIQUE],
            "local_remainder": list(REMAINDER),
            "global_remainder": [outside[index] for index in REMAINDER],
            "propagated_support_supersets": list(propagated),
            "basis_masks": basis_masks,
            "pair_order": [list(pair) for pair in PAIR_ORDER],
            "pair_targets": list(pair_targets),
            "pair_target_string": "".join(map(str, pair_targets)),
        },
        "entry_reasons": {
            "required_unit_edges": reasons.count("required_unit_edge"),
            "disjoint_support_zeros": reasons.count(
                "disjoint_propagated_support_supersets"
            ),
            "unresolved_optional_entries": 0,
        },
    }
    return certificate


def mask_sum(mask: int, variables: Sequence[sp.Symbol]) -> sp.Expr:
    return sum((variables[index] for index in bits(mask)), sp.Integer(0))


def polynomial_zero(expression: sp.Expr, variables: Iterable[sp.Symbol], label: str) -> None:
    require(
        sp.Poly(sp.together(expression), *variables, domain=sp.QQ).is_zero,
        label,
    )


def independent_algebra() -> dict:
    variables = sp.symbols("a b c d e f g")
    a, b, c, d, e, f, g = variables
    total7 = 1 + sum(variables)
    sums = [mask_sum(mask, variables) for mask in EXPECTED_BASIS_MASKS]
    equations: dict[tuple[int, int], sp.Expr] = {}
    for pair, target in zip(PAIR_ORDER, EXPECTED_PAIR_TARGETS, strict=True):
        first, second = pair
        intersection = mask_sum(
            EXPECTED_BASIS_MASKS[first] & EXPECTED_BASIS_MASKS[second], variables
        )
        equations[pair] = sp.expand(
            total7 * intersection - sums[first] * sums[second] - target * total7
        )

    positive_sum = b + f + g
    polynomial_zero(
        equations[(0, 2)] - equations[(2, 3)] + (a - c) * positive_sum,
        variables,
        "a=c identity",
    )
    polynomial_zero(
        equations[(2, 3)] - equations[(2, 4)] + (c - d) * positive_sum,
        variables,
        "c=d identity",
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
    for pair, selected in zip(PAIR_ORDER, mapping, strict=True):
        polynomial_zero(
            equations[pair].subs(substitution) - reduced[selected],
            (t, b, e, f, g),
            f"ten-to-four reduction {pair}",
        )
    e0, e1, e2, e3 = reduced

    # Independent ideal-membership route, distinct from the producer's
    # sequential P/u elimination.
    a0 = b * (f + t) - f - t**2 - 2 * t
    b0 = b * f + b * t + b - f * t - 2 * f - g * t - 4 * t - 2
    polynomial_zero(
        (f + 2 * t + 1) * a0 - (f + t) * e0 - (f + 2 * t) * e1,
        (t, b, e, f, g),
        "A0 ideal identity",
    )
    polynomial_zero(
        (b + t - 1) * b0 - (b - 2) * (e0 + e1) - t * (e2 + e3),
        (t, b, e, f, g),
        "B0 ideal identity",
    )

    b_solution = sp.cancel((f + t**2 + 2 * t) / (f + t))
    g_solution = sp.cancel(
        (-f**2 * t - f**2 - 3 * f * t - f + t**3 - t**2)
        / (t * (f + t))
    )
    e_solution = sp.cancel((f * t + f - t**2 + 2 * t + 1) / t)
    solved = {b: b_solution, e: e_solution, g: g_solution}
    polynomial_zero(e0.subs(solved), (t, f), "solved E0")
    polynomial_zero(e1.subs(solved), (t, f), "solved E1")
    final = t**2 - t - f * (t + 1)
    multiplier = (f + 2 * t + 1) ** 2 / (t * (f + t))
    polynomial_zero(e2.subs(solved) - multiplier * final, (t, f), "E2 final")
    polynomial_zero(e3.subs(solved) + multiplier * final, (t, f), "E3 final")

    f_solution = sp.cancel(t * (t - 1) / (t + 1))
    final_g = sp.factor(g_solution.subs(f, f_solution))
    polynomial_zero(
        final_g + (t - 1) * (t + 1) / (2 * t**2),
        (t,),
        "negative g formula",
    )

    boundary = dict(zip(variables, (1, 3, 1, 1, 2, 0, 0), strict=True))
    require(
        all(sp.expand(equation.subs(boundary)) == 0 for equation in equations.values()),
        "nonnegative boundary control",
    )
    return {
        "method": "independent_ideal_membership_A0_B0_elimination",
        "variables": [str(variable) for variable in variables],
        "pair_equations": [str(equations[pair]) for pair in PAIR_ORDER],
        "reduced_equations": [str(expression) for expression in reduced],
        "identities_checked_coefficientwise": 2 + len(PAIR_ORDER) + 4 + 1,
        "strict_sign_chain": [
            "a=c=d=t because b+f+g>0",
            "A0=0 gives b=1+t*(t+1)/(f+t)>1",
            "therefore b+t-1>0 and the B0 identity gives B0=0",
            "solving A0=B0=E0 gives the displayed b,g,e rational forms",
            "E2=0 gives f=t*(t-1)/(t+1)",
            "f>0 and t>0 force t>1",
            "g=-(t-1)*(t+1)/(2*t^2)<0, contradicting g>0",
        ],
        "solved_forms": {
            "b": str(b_solution),
            "e": str(e_solution),
            "g": str(g_solution),
            "f_after_final_factor": str(f_solution),
            "g_after_final_factor": str(final_g),
        },
        "boundary_control": {
            "point": [1, 3, 1, 1, 2, 0, 0],
            "all_ten_equations_zero": True,
            "strict_domain": False,
        },
        "status": "EXACT_STRICT_POSITIVITY_CONTRADICTION",
    }


def validate_report_certificate(report: dict, rebuilt: dict, algebra: dict) -> None:
    certificate = report.get("certificate", {})
    target = certificate.get("target", {})
    require(target == rebuilt["target"], "target certificate")
    archived_quantifier = certificate.get("quantifier", {})
    for key, expected in rebuilt["quantifier"].items():
        require(archived_quantifier.get(key) == expected, f"quantifier field {key}")
    archived_pattern = certificate.get("exact_pattern", {})
    for key, expected in rebuilt["pattern"].items():
        require(archived_pattern.get(key) == expected, f"pattern field {key}")
    require(
        certificate.get("algebra", {}).get("pair_equations")
        == algebra["pair_equations"],
        "archived pair equations",
    )
    require(
        certificate.get("algebra", {}).get("reduced_equations")
        == algebra["reduced_equations"],
        "archived reduced equations",
    )
    require(
        certificate.get("rank_argument")
        == {
            "normalized_gram_rank_upper_bound": 7,
            "required_clique_order": 7,
            "clique_principal_block_positive_definite": True,
            "schur_complement_for_five_remainder_vertices": "identically zero",
            "clique_inverse_method": "Sherman--Morrison over exact rationals",
            "strict_domain": "a,b,c,d,e,f,g > 0",
        },
        "rank argument",
    )


def verify_report(
    report_path: Path,
    expected_report_sha256: str,
    *,
    enforce_source_boundary: bool = True,
) -> dict:
    report_path = report_path.resolve()
    report_hash = sha256(report_path)
    require(report_hash == expected_report_sha256, "report hash")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    validate_report_schema(report)
    provenance = (
        validate_source_boundary(report)
        if enforce_source_boundary
        else {"status": "SKIPPED_FOR_TEST"}
    )
    manifest, star, star_check = load_upstream()
    require(
        tuple(map(int, star["summary"]["ordered_survivor_indices"]))
        == EXPECTED_INPUT,
        "independently loaded input list",
    )
    rebuilt = independent_quantifier_and_pattern(manifest, star)
    algebra = independent_algebra()
    validate_report_certificate(report, rebuilt, algebra)
    return {
        "schema": 1,
        "kind": "d6_k7_schur_3949382_increment_verification",
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": report_hash},
        "semantics": EXPECTED_SEMANTICS,
        "source_boundary": provenance,
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "checked": {
            "target_index": TARGET_INDEX,
            "raw_eligible_covers": rebuilt["quantifier"]["raw_eligible_covers"],
            "prior_eliminated_raw_covers": rebuilt["quantifier"][
                "prior_layer_eliminated_raw_covers"
            ],
            "upstream_current_covers": len(EXPECTED_COVERS),
            "upstream_infeasible_covers": 2,
            "new_families": 1,
            "basis_entries": len(REMAINDER) * len(CLIQUE),
            "remainder_pair_entries": len(PAIR_ORDER),
            "unresolved_optional_entries": rebuilt["entry_reasons"][
                "unresolved_optional_entries"
            ],
            "upstream_independently_checked_eligible_covers": star_check["checked"][
                "eligible_covers"
            ],
            "upstream_independently_checked_current_passing_families": star_check[
                "checked"
            ]["current_passing_families"],
            "input_graphs": len(EXPECTED_INPUT),
            "output_graphs": len(EXPECTED_INPUT) - 1,
        },
        "independent_algebra": algebra,
        "conclusion": {
            "rejected_indices": [TARGET_INDEX],
            "ordered_survivor_indices_sha256": EXPECTED_OUTPUT_SHA256,
            "floating_point_used": False,
        },
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "finished_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "verifier_sha256": sha256(Path(__file__).resolve()),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k7_schur_3949382_increment_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_schur_3949382_increment_verification.json",
    )
    args = parser.parse_args()
    result = verify_report(args.report, args.report_sha256)
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "sha256": sha256(args.output.resolve()),
                "status": result["status"],
                "report_sha256": result["report"]["sha256"],
                "rejected_indices": result["conclusion"]["rejected_indices"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
