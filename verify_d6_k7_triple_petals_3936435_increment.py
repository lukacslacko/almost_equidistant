#!/usr/bin/env python3
"""Independent checker for the graph-3936435 triple-petal increment.

This verifier imports neither the producer nor any discovery probe.  It
reloads the hash-pinned current K7 and one/two-star boundaries, reconstructs
the target graph and K7 seed, enumerates the 255 eligible covers in a
different order, checks both surviving support families, and verifies the
quadratic-field contradiction with symbolic polynomial reduction.
"""

from __future__ import annotations

import argparse
import ast
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
from typing import Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 3_936_435
SEED = (1, 3, 6, 8, 11, 13, 17)
SEED_MASK = 141_642
FIVE_LOCALS = (0, 1, 3, 6, 9)
HUB_LOCALS = (0, 6)
PETAL_LOCALS = (1, 3, 9)
CENTER = 1
LEAVES = (3, 2, 5)
HUB_MASKS = (60, 45)
PETAL_MASKS = (10, 6, 34)
FIVE_MASKS = (60, 10, 6, 45, 34)
CURRENT_COVERS = (0, 2048)
EXPECTED_RAW_COVERS = 255
EXPECTED_RAW_COVERS_SHA256 = (
    "cf605db4a5fed03da735e569687d1e84d2bdde24c45d95a9fc4b61490e1a491d"
)
EXPECTED_ADJACENCY_SHA256 = (
    "15af539b50d518368164e3908f8609159d7c4ac04f1ab45be695eec27d7b9c62"
)
EXPECTED_INPUT = (
    316173,
    2581209,
    3648882,
    3729907,
    3935560,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
)
EXPECTED_INPUT_SHA256 = (
    "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c"
)
EXPECTED_OUTPUT_SHA256 = (
    "380936282d98f2e561c71680a04d04401fabc032e8a00cebf12ae4d18967a874"
)
EXPECTED_PROPAGATED = {
    0: (60, 10, 121, 6, 67, 117, 45, 82, 93, 34, 64, 127),
    2048: (60, 10, 121, 6, 67, 117, 45, 82, 93, 34, 64),
}
EXPECTED_Z_SUPPORTS = {0: (), 2048: (63,)}
UPSTREAM = {
    "d6_current_residue_manifest_v8.json": (
        "9ea10a7794f033e66152c477a130f11c6c2862e87b4bdd705c14021521c18285"
    ),
    "d6_current_residue_manifest_v8_verification.json": (
        "117ac6833a24eb69cec7a514445bb4c2913c3f3ac7c413aa81fad19fe305602e"
    ),
    "d6_k7_one_two_star_increment_report.json": (
        "12b3d18b1ea81961f58d831d4c7c322fbceb2300e7533b64161e8011fbe9d1ec"
    ),
    "d6_k7_one_two_star_increment_verification.json": (
        "c6f40578685036cb1156cb0ea8bf06d6004a70d5f376aff5ac1e16a1c0ca09eb"
    ),
}
PACKAGE_SOURCES = (
    "build_d6_k7_triple_petals_3936435_increment.py",
    "verify_d6_k7_triple_petals_3936435_increment.py",
    "test_d6_k7_triple_petals_3936435_increment.py",
    "d6_k7_triple_petals_3936435_increment.md",
)
REPORT_KIND = "d6_k7_triple_petals_3936435_exact_increment"
VERIFICATION_KIND = "d6_k7_triple_petals_3936435_increment_verification"
SEMANTICS = {
    "actual_supports_are_subsets_of_propagated_masks": True,
    "all_actual_cover_support_branches_for_seed_are_quantified": True,
    "candidate_nonedges_optional": True,
    "floating_point_enters_rejection": False,
    "one_infeasible_required_K7_seed_rejects_graph": True,
    "only_required_edges_enter_inner_product_equations": True,
    "zero_factor_vertices_excluded_from_normalization": True,
}
STAR_SEMANTICS = {
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
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("ascii")
    ).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def bits_descending(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(mask.bit_length() - 1, -1, -1) if mask & (1 << index))


def validate_adjacency(adjacency: Sequence[int], order: int) -> None:
    require(len(adjacency) == order, "graph order")
    full = (1 << order) - 1
    for first in range(order):
        require(not adjacency[first] & ~full, "adjacency outside order")
        require(not adjacency[first] & (1 << first), "loop")
        for second in range(first + 1, order):
            require(
                bool(adjacency[first] & (1 << second))
                == bool(adjacency[second] & (1 << first)),
                "asymmetry",
            )


def clique(adjacency: Sequence[int], vertices: Sequence[int]) -> bool:
    chosen = tuple(vertices)
    return all(
        adjacency[chosen[right]] & (1 << chosen[left])
        for right in range(1, len(chosen)) for left in range(right)
    )


def independent_seed_data(adjacency: Sequence[int]):
    require(clique(adjacency, SEED), "required seed")
    outside = tuple(sorted(set(range(19)).difference(SEED)))
    coordinate = {seed_vertex: index for index, seed_vertex in enumerate(SEED)}
    defects = []
    for vertex in outside:
        mask = 0
        for seed_vertex in SEED:
            if not adjacency[vertex] & (1 << seed_vertex):
                mask |= 1 << coordinate[seed_vertex]
        defects.append(mask)
    lorentz_edges = []
    for first, second in combinations(range(12), 2):
        if (
            adjacency[outside[first]] & (1 << outside[second])
            and not defects[first] & defects[second]
        ):
            lorentz_edges.append((first, second))
    eligible_vertices = tuple(
        local for local, mask in enumerate(defects) if mask.bit_count() >= 3
    )
    return outside, tuple(defects), tuple(lorentz_edges), eligible_vertices


def independently_enumerate_covers(
    lorentz_edges: Sequence[tuple[int, int]], eligible_vertices: Sequence[int]
) -> tuple[int, ...]:
    """Combination-ordered replay, then numerically sort the cover masks."""

    found = []
    for size in range(min(7, len(eligible_vertices)) + 1):
        for chosen in combinations(eligible_vertices, size):
            mask = sum(1 << local for local in chosen)
            if all(mask & (1 << first) or mask & (1 << second) for first, second in lorentz_edges):
                found.append(mask)
    return tuple(sorted(found))


def independent_nonzero_mask_map(zmask: int, witness: dict) -> dict[int, int]:
    surviving = [local for local in range(12) if not (zmask >> local) & 1]
    values = list(map(int, witness["propagated_masks"]))
    require(len(surviving) == len(values), "nonzero support list length")
    return {surviving[position]: values[position] for position in range(len(values))}


def assert_import_independence() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "build_d6_k7_triple_petals_3936435_increment",
        "probe_d6_k7_triple_petals_3936435",
    }
    require(not imports & forbidden, f"checker imports producer/probe: {imports & forbidden}")


def load_upstream() -> tuple[dict, dict, dict]:
    for name, expected in UPSTREAM.items():
        require(sha256(ROOT / name) == expected, f"upstream hash: {name}")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest_v8.json").read_text(encoding="utf-8")
    )
    manifest_check = json.loads(
        (ROOT / "d6_current_residue_manifest_v8_verification.json").read_text(
            encoding="utf-8"
        )
    )
    parent = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_report.json").read_text(encoding="utf-8")
    )
    parent_check = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v8"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "manifest boundary",
    )
    require(
        manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v8.json"],
        "manifest independent verification",
    )
    require(
        parent.get("kind")
        == "d6_k7_one_free_neighbour_two_free_center_increment"
        and parent.get("status") == "COMPLETE"
        and parent.get("semantics") == STAR_SEMANTICS,
        "parent cover boundary",
    )
    require(
        parent_check.get("kind") == "d6_k7_one_two_star_increment_verification"
        and parent_check.get("status") == "PASS"
        and parent_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_one_two_star_increment_report.json"]
        and parent_check.get("semantics") == STAR_SEMANTICS
        and parent_check.get("checked", {}).get("eligible_covers") == 19_932
        and parent_check.get("checked", {}).get("current_passing_families") == 88,
        "parent independent verification",
    )
    return manifest, parent, parent_check


def reduce_sqrt7(expression: sp.Expr, radical: sp.Symbol) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    require(not denominator.has(radical), "unexpected radical denominator")
    remainder = sp.rem(sp.Poly(numerator, radical), sp.Poly(radical**2 - 7, radical))
    return sp.cancel(remainder.as_expr() / denominator)


def q7_json(expression: sp.Expr, radical: sp.Symbol) -> list[list[int]]:
    """Canonical ``[rational, sqrt(7)-coefficient]`` serialization."""

    reduced = sp.Poly(reduce_sqrt7(expression, radical), radical)
    require(reduced.degree() <= 1, "Q(sqrt(7)) serialization degree")
    rational = sp.Rational(reduced.coeff_monomial(1))
    radical_coefficient = sp.Rational(reduced.coeff_monomial(radical))
    return [
        [int(rational.p), int(rational.q)],
        [int(radical_coefficient.p), int(radical_coefficient.q)],
    ]


def independent_algebra() -> dict:
    """Polynomial-identity route, separate from the producer's Q7 class."""

    x, y, z, radical, leaf = sp.symbols("x y z r lambda")
    xy = x * y - 1
    xz = x * z - 1
    yz = y * z - 1
    require(sp.expand(x * yz - z * xy - (z - x)) == 0, "central z=x identity")
    require(sp.expand(x * yz - y * xz - (y - x)) == 0, "central y=x identity")
    z_minus_x = x * yz - z * xy
    require(
        sp.expand(xz - x * z_minus_x - (x**2 - 1)) == 0,
        "central x-square identity",
    )
    hub_a = sp.symbols("a2 a3 a5")
    hub_d = sp.symbols("d2 d3 d5")
    hub_a_errors = tuple(leaf * coordinate - 1 for coordinate in hub_a)
    hub_d_errors = tuple(leaf * coordinate - 1 for coordinate in hub_d)
    hub_edge = sum(
        first * second for first, second in zip(hub_a, hub_d, strict=True)
    ) - 1
    # This ideal identity mechanizes the full hub implication without division:
    # six equations lambda*a_i=lambda*d_i=1 plus hub_edge=0 force
    # 3-lambda^2=0.
    hub_ideal_identity = sp.expand(
        leaf**2 * hub_edge
        - (3 - leaf**2)
        - sum(
            first * second + first + second
            for first, second in zip(
                hub_a_errors, hub_d_errors, strict=True
            )
        )
    )
    require(hub_ideal_identity == 0, "six-edge hub ideal identity")
    cases = []
    for epsilon in (-1, 1):
        leaf_value = (radical - epsilon) / 2
        leaf_conjugate = (-radical - epsilon) / 2
        leaf_norm = reduce_sqrt7(leaf_value * leaf_conjugate, radical)
        require(leaf_norm == -sp.Rational(3, 2), "nonzero leaf norm")
        inverse_value = (radical + epsilon) / 3
        require(
            reduce_sqrt7(leaf_value * inverse_value, radical) == 1,
            "independent leaf inverse",
        )
        diagonal = (
            epsilon**2 + leaf**2
            - 1
            - (epsilon + leaf - radical) ** 2
        )
        require(
            reduce_sqrt7(diagonal.subs(leaf, leaf_value), radical) == 0,
            "petal diagonal substitution",
        )
        gap = reduce_sqrt7(leaf_value**2 - 3, radical)
        expected_gap = -1 - sp.Rational(epsilon, 2) * radical
        require(sp.expand(gap - expected_gap) == 0, "leaf-square gap")
        # If gap vanished, epsilon*r=-2. Multiplying by its conjugate would
        # force r^2-4=0, whose remainder modulo r^2-7 is the nonzero integer 3.
        norm_contradiction = reduce_sqrt7(
            (epsilon * radical + 2) * (epsilon * radical - 2), radical
        )
        require(norm_contradiction == 3, "radical norm contradiction")
        substituted_hub_errors = tuple(
            reduce_sqrt7(error.subs(leaf, leaf_value).subs(coordinate, inverse_value), radical)
            for error, coordinate in zip(
                (*hub_a_errors, *hub_d_errors), (*hub_a, *hub_d), strict=True
            )
        )
        require(
            all(error == 0 for error in substituted_hub_errors),
            "six independently forced hub coordinates",
        )
        forced_hub_dot = reduce_sqrt7(3 * inverse_value**2, radical)
        cleared_hub_edge = reduce_sqrt7(
            leaf_value**2 * (forced_hub_dot - 1), radical
        )
        three_minus_leaf_square = reduce_sqrt7(3 - leaf_value**2, radical)
        require(
            cleared_hub_edge == three_minus_leaf_square,
            "cleared hub-edge substitution identity",
        )
        producer_case = {
            "epsilon": epsilon,
            "leaf_value": q7_json(leaf_value, radical),
            "leaf_inverse": q7_json(inverse_value, radical),
            "leaf_square": q7_json(leaf_value**2, radical),
            "leaf_square_minus_three": q7_json(gap, radical),
            "two_support_diagonal_factor": q7_json(3, radical),
            "six_hub_petal_products": [q7_json(1, radical) for _ in range(6)],
            "hub_hub_dot_forced_by_six_petal_edges": q7_json(
                forced_hub_dot, radical
            ),
            "leaf_square_times_hub_dot_minus_one": q7_json(
                cleared_hub_edge, radical
            ),
            "three_minus_leaf_square": q7_json(
                three_minus_leaf_square, radical
            ),
        }
        cases.append({
            "epsilon": epsilon,
            "lambda": str(leaf_value),
            "lambda_conjugate_product": str(leaf_norm),
            "lambda_inverse": str(inverse_value),
            "six_hub_petal_substitutions": [str(value) for value in substituted_hub_errors],
            "forced_hub_dot": str(forced_hub_dot),
            "cleared_hub_edge": str(cleared_hub_edge),
            "lambda_square_minus_three": str(gap),
            "conjugate_product_mod_r_squared_minus_7": str(norm_contradiction),
            "producer_case_expected": producer_case,
        })
    return {
        "method": "symbolic_polynomial_reduction_modulo_r_squared_minus_7",
        "central_equalities_checked_by_ideal_combinations": 3,
        "hub_equations_checked": 7,
        "hub_ideal_identity": str(hub_ideal_identity),
        "hub_ideal_consequence": "six_hub_petal_edges_and_hub_edge_imply_lambda_squared_equals_three",
        "petal_diagonal_cases": cases,
        "hub_intersection_leaves": list(LEAVES),
        "status": "EXACT_CONTRADICTION_FOR_BOTH_SIGNS",
    }


def independent_positive_controls() -> dict:
    codes = [value for value in range(1 << 5) if value.bit_count() & 1]
    adjacency = [0] * 18
    for first, second in combinations(range(16), 2):
        if (codes[first] ^ codes[second]).bit_count() == 2:
            adjacency[first] |= 1 << second
            adjacency[second] |= 1 << first
    for vertex in range(16):
        adjacency[vertex] |= (1 << 16) | (1 << 17)
        adjacency[16] |= 1 << vertex
        adjacency[17] |= 1 << vertex
    validate_adjacency(adjacency, 18)
    seeds = sum(clique(adjacency, chosen) for chosen in combinations(range(18), 7))
    require(seeds == 0, "realizable 18-point K7 seeds")

    # The exact combinatorial template must disappear when the hub pair is an
    # optional nonedge; the checker never assigns that pair target zero.
    complete = [0] * 10
    for first, second in combinations(FIVE_LOCALS, 2):
        complete[first] |= 1 << second
        complete[second] |= 1 << first
    first, second = HUB_LOCALS
    optional = list(complete)
    optional[first] ^= 1 << second
    optional[second] ^= 1 << first
    require(clique(complete, FIVE_LOCALS), "synthetic required K5")
    require(not clique(optional, FIVE_LOCALS), "optional hub edge control")
    return {
        "known_realizable_18": {"passed": True, "K7_seeds": seeds},
        "optional_hub_nonedge_is_not_an_equation": True,
    }


def independently_rebuild(manifest: dict, parent: dict) -> dict:
    indices = tuple(int(row["index"]) for row in manifest["classes"]["K7"]["graphs"])
    require(indices == EXPECTED_INPUT, "input indices")
    require(stable_hash(list(indices)) == EXPECTED_INPUT_SHA256, "input hash")
    graph = next(
        row for row in manifest["classes"]["K7"]["graphs"]
        if int(row["index"]) == TARGET_INDEX
    )
    adjacency = tuple(map(int, graph["adjacency"]))
    validate_adjacency(adjacency, 19)
    require(stable_hash(list(adjacency)) == EXPECTED_ADJACENCY_SHA256, "adjacency hash")
    outside, defects, lorentz_edges, eligible = independent_seed_data(adjacency)
    covers = independently_enumerate_covers(lorentz_edges, eligible)
    require(len(covers) == EXPECTED_RAW_COVERS, "eligible cover count")
    require(stable_hash(list(covers)) == EXPECTED_RAW_COVERS_SHA256, "eligible cover hash")

    graph_record = next(row for row in parent["records"] if int(row["index"]) == TARGET_INDEX)
    require(graph_record.get("decision") == "SURVIVOR", "parent graph status")
    seed_record = next(
        row for row in graph_record["seeds"] if int(row["seed_mask"]) == SEED_MASK
    )
    require(tuple(seed_record["seed"]) == SEED, "parent seed identity")
    archived = seed_record["current_covers"]
    require(tuple(int(row["zmask"]) for row in archived) == CURRENT_COVERS, "parent covers")
    cover_rows = []
    for row in archived:
        zmask = int(row["zmask"])
        require(
            row.get("status") == "PASSING"
            and int(row["current_passing_families"]) == 1
            and int(row["star_passing_families"]) == 1,
            f"parent family count {zmask}",
        )
        witness = row["first_passing_witness"]
        require(tuple(map(int, witness["propagated_masks"])) == EXPECTED_PROPAGATED[zmask], "propagated")
        require(tuple(map(int, witness["z_supports"])) == EXPECTED_Z_SUPPORTS[zmask], "Z supports")
        masks = independent_nonzero_mask_map(zmask, witness)
        require(not any(zmask & (1 << local) for local in FIVE_LOCALS), "pattern enters Z")
        relevant = tuple(masks[local] for local in FIVE_LOCALS)
        require(relevant == FIVE_MASKS, "five support supersets")
        globals5 = tuple(outside[local] for local in FIVE_LOCALS)
        require(clique(adjacency, globals5), "required K5")

        petal_masks = tuple(masks[local] for local in PETAL_LOCALS)
        hub_masks = tuple(masks[local] for local in HUB_LOCALS)
        require(petal_masks == PETAL_MASKS and hub_masks == HUB_MASKS, "labeled masks")
        center_bit = 1 << CENTER
        leaf_bits = tuple(1 << leaf for leaf in LEAVES)
        require(
            all(
                petal_masks[first] & petal_masks[second] == center_bit
                for first, second in combinations(range(3), 2)
            ),
            "petal center intersections",
        )
        require(
            all(hub & petal_masks[number] == leaf_bits[number] for hub in hub_masks for number in range(3)),
            "hub-petal leaf intersections",
        )
        require(hub_masks[0] & hub_masks[1] == sum(leaf_bits), "hub intersection")
        cover_rows.append({
            "zmask": zmask,
            "five_global_vertices": list(globals5),
            "five_masks": list(relevant),
            "all_ten_required_edges": True,
        })

    return {
        "indices": indices,
        "adjacency": adjacency,
        "outside": outside,
        "defects": defects,
        "covers": covers,
        "current_cover_rows": cover_rows,
    }


def validate_report_schema(report: dict) -> None:
    output = [index for index in EXPECTED_INPUT if index != TARGET_INDEX]
    require(report.get("schema") == 1, "report schema")
    require(report.get("kind") == REPORT_KIND, "report kind")
    require(report.get("status") == "COMPLETE_EXACT_REJECTION", "report status")
    require(
        report.get("claim")
        == "Graph 3936435 has no realization by distinct points in R^6 with all candidate edges at unit distance.",
        "report claim",
    )
    require(report.get("semantics") == SEMANTICS, "semantics")
    require(report.get("upstream_sha256") == dict(sorted(UPSTREAM.items())), "upstream roots")
    require(
        report.get("input", {}).get("ordered_indices") == list(EXPECTED_INPUT)
        and report.get("input", {}).get("ordered_indices_sha256") == EXPECTED_INPUT_SHA256,
        "input boundary",
    )
    summary = report.get("summary", {})
    require(
        summary.get("rejected_indices") == [TARGET_INDEX]
        and summary.get("graphs_rejected") == 1
        and summary.get("ordered_survivor_indices") == output
        and summary.get("ordered_survivor_indices_sha256") == EXPECTED_OUTPUT_SHA256
        and summary.get("graphs_surviving") == len(output),
        "summary",
    )
    require(report.get("certificate_sha256") == stable_hash(report.get("certificate")), "certificate hash")


def validate_source_boundary(report: dict) -> dict:
    sources = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    require(report.get("source_sha256") == sources, "working source hashes")
    launch = report.get("execution", {}).get("git", {})
    require(
        launch.get("branch") == "codex/dimension6"
        and launch.get("source_sha256") == sources
        and launch.get("proof_and_checker_sources_equal_committed_blobs") is True
        and launch.get("tracked_clean") is True,
        "launch source boundary",
    )
    commit = launch.get("commit")
    require(isinstance(commit, str) and len(commit) == 40, "launch commit")
    for name, expected in sources.items():
        blob = subprocess.run(
            ["git", "show", f"{commit}:{name}"], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        ).stdout
        require(hashlib.sha256(blob).hexdigest() == expected, f"committed source: {name}")
    porcelain = "\n".join(launch.get("porcelain_lines", ()))
    require(
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        == launch.get("porcelain_sha256"),
        "porcelain binding",
    )
    return {"commit": commit, "branch": launch["branch"], "sources": len(sources)}


def validate_certificate(report: dict, rebuilt: dict, algebra: dict) -> None:
    certificate = report["certificate"]
    require(certificate["target"]["index"] == TARGET_INDEX, "certificate target")
    require(certificate["target"]["seed"] == list(SEED), "certificate seed")
    require(certificate["target"]["outside"] == list(rebuilt["outside"]), "certificate outside")
    quantifier = certificate["quantifier"]
    require(
        quantifier["raw_eligible_covers"] == len(rebuilt["covers"])
        and quantifier["raw_eligible_covers_sha256"] == stable_hash(list(rebuilt["covers"]))
        and quantifier["prior_layer_eliminated_raw_covers"] == 253
        and quantifier["upstream_current_covers"] == list(CURRENT_COVERS)
        and quantifier["new_families_checked"] == 2
        and quantifier["new_families_rejected"] == 2,
        "certificate quantifier",
    )
    require(len(certificate["covers"]) == 2, "certificate cover count")
    for archived, rebuilt_row in zip(certificate["covers"], rebuilt["current_cover_rows"], strict=True):
        require(archived["zmask"] == rebuilt_row["zmask"], "cover zmask")
        require(archived["five_global_vertices"] == rebuilt_row["five_global_vertices"], "cover vertices")
        require(archived["five_propagated_masks"] == rebuilt_row["five_masks"], "cover masks")
        require(archived["all_ten_pairs_required_unit_edges"] is True, "cover edges")
        require(archived["contradiction"] == "triple_petal_Q_sqrt7", "cover contradiction")
    expected_algebra_text = {
        "field": "Q(sqrt(7)) with positive sqrt(7)",
        "central_product_deduction": (
            "The three required petal edges give xy=xz=yz=1, hence "
            "x=y=z=epsilon with epsilon in {-1,+1}."
        ),
        "two_support_diagonal_identity": (
            "For normalized w, ||w||^2=1+(sum(w)-sqrt(7))^2. "
            "A petal (epsilon,lambda) therefore has "
            "lambda=(sqrt(7)-epsilon)/2."
        ),
        "hub_edge_deduction": (
            "Each hub-petal edge forces the corresponding hub leaf to "
            "1/lambda (lambda is nonzero). All six products are checked. "
            "The two hubs share exactly three leaves, so their required edge "
            "forces 3/lambda^2=1; clearing the nonzero lambda^2 gives "
            "lambda^2=3."
        ),
        "conclusion": (
            "lambda^2=2-epsilon*sqrt(7)/2 is never 3; equality would give "
            "epsilon*sqrt(7)=-2 and hence 7=4."
        ),
    }
    for key, expected in expected_algebra_text.items():
        require(certificate["algebra"].get(key) == expected, f"algebra text: {key}")
    require(
        certificate["algebra"]["status"] == algebra["status"]
        == "EXACT_CONTRADICTION_FOR_BOTH_SIGNS",
        "algebra conclusion",
    )
    archived_cases = certificate["algebra"]["cases"]
    require([row["epsilon"] for row in archived_cases] == [-1, 1], "archived sign cases")
    require([row["epsilon"] for row in algebra["petal_diagonal_cases"]] == [-1, 1], "independent sign cases")
    independently_expected_cases = [
        row["producer_case_expected"] for row in algebra["petal_diagonal_cases"]
    ]
    require(
        archived_cases == independently_expected_cases,
        "all archived Q(sqrt(7)) algebra values",
    )
    require(
        algebra["hub_ideal_identity"] == "0"
        and algebra["hub_equations_checked"] == 7,
        "independent full hub ideal implication",
    )
    require(
        certificate.get("graph_rejection_quantifier")
        == (
            "Every realization must realize this required K7 seed. The "
            "hash-pinned parent exhausts all 255 eligible covers, eliminates "
            "253, and leaves exactly the two one-family covers checked here; "
            "both exact support/value systems are contradictory."
        ),
        "graph rejection quantifier",
    )


def verify_report(
    report_path: Path, expected_report_sha256: str, *, enforce_source_boundary: bool = True
) -> dict:
    assert_import_independence()
    report_path = report_path.resolve()
    report_hash = sha256(report_path)
    require(report_hash == expected_report_sha256, "report hash")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    validate_report_schema(report)
    source_boundary = (
        validate_source_boundary(report)
        if enforce_source_boundary else {"status": "SKIPPED_FOR_TEST"}
    )
    manifest, parent, parent_check = load_upstream()
    rebuilt = independently_rebuild(manifest, parent)
    algebra = independent_algebra()
    controls = independent_positive_controls()
    validate_certificate(report, rebuilt, algebra)
    require(report.get("controls", {}).get("known_realizable_18") == controls["known_realizable_18"], "positive control")
    return {
        "schema": 1,
        "kind": VERIFICATION_KIND,
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": report_hash},
        "semantics": SEMANTICS,
        "source_boundary": source_boundary,
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "checks": {
            "import_independence": True,
            "current_K7_parent_boundary": True,
            "independent_raw_cover_replay": True,
            "all_255_raw_covers_accounted": True,
            "both_surviving_covers_checked": True,
            "both_unique_support_families_checked": True,
            "propagated_masks_used_only_as_supersets": True,
            "only_required_edges_used": True,
            "independent_symbolic_Q_sqrt7_arithmetic": True,
            "known_realizable_18_control": True,
        },
        "checked": {
            "input_graphs": len(EXPECTED_INPUT),
            "output_graphs": len(EXPECTED_INPUT) - 1,
            "target_index": TARGET_INDEX,
            "raw_eligible_covers": len(rebuilt["covers"]),
            "prior_eliminated_raw_covers": 253,
            "new_covers": 2,
            "new_families": 2,
            "parent_independently_checked_eligible_covers": parent_check["checked"]["eligible_covers"],
            "parent_independently_checked_current_passing_families": parent_check["checked"]["current_passing_families"],
        },
        "independent_algebra": algebra,
        "controls": controls,
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
        "--report", type=Path,
        default=ROOT / "d6_k7_triple_petals_3936435_increment_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_triple_petals_3936435_increment_verification.json",
    )
    args = parser.parse_args()
    result = verify_report(args.report, args.report_sha256)
    atomic_json(args.output.resolve(), result)
    print(json.dumps({
        "status": result["status"],
        "rejected": result["conclusion"]["rejected_indices"],
        "report_sha256": result["report"]["sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
