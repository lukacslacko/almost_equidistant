#!/usr/bin/env python3
"""Build the exact triple-petal increment for K7 graph 3936435.

The package starts from the hash-pinned, independently verified current K7
residue and the independently verified one/two-free-star cover boundary.  For
one required K7 seed, that boundary leaves two covers and one propagated
support family on each.  The same five nonzero-factor vertices in both
families form a required K5 whose support pattern is contradictory in
Q(sqrt(7)).

Candidate nonedges are never assigned a distance.  Propagated masks are used
only as supersets of actual supports.  All rejecting arithmetic is exact.
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
from dataclasses import dataclass
from datetime import UTC, datetime
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 3_936_435
SEED = (1, 3, 6, 8, 11, 13, 17)
SEED_MASK = 141_642
PETAL_LOCALS = (1, 3, 9)
HUB_LOCALS = (0, 6)
FIVE_LOCALS = (0, 1, 3, 6, 9)
CENTER_COORDINATE = 1
LEAF_COORDINATES = (3, 2, 5)
EXPECTED_PETAL_MASKS = (10, 6, 34)
EXPECTED_HUB_MASKS = (60, 45)
EXPECTED_RELEVANT_MASKS = (60, 10, 6, 45, 34)
EXPECTED_CURRENT_COVERS = (0, 2048)
EXPECTED_Z_SUPPORTS = {0: (), 2048: (63,)}
EXPECTED_PROPAGATED_BY_COVER = {
    0: (60, 10, 121, 6, 67, 117, 45, 82, 93, 34, 64, 127),
    2048: (60, 10, 121, 6, 67, 117, 45, 82, 93, 34, 64),
}
EXPECTED_RAW_COVERS = 255
EXPECTED_RAW_COVERS_SHA256 = (
    "cf605db4a5fed03da735e569687d1e84d2bdde24c45d95a9fc4b61490e1a491d"
)
EXPECTED_ADJACENCY_SHA256 = (
    "15af539b50d518368164e3908f8609159d7c4ac04f1ab45be695eec27d7b9c62"
)
EXPECTED_INPUT_INDICES = (
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
REPORT = ROOT / "d6_k7_triple_petals_3936435_increment_report.json"
REPORT_KIND = "d6_k7_triple_petals_3936435_exact_increment"
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
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def bits(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(mask.bit_length()) if mask & (1 << index))


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


@dataclass(frozen=True)
class Q7:
    """An exact element ``a+b*sqrt(7)``."""

    rational: Fraction = Fraction(0)
    radical: Fraction = Fraction(0)

    @classmethod
    def of(cls, value: int | Fraction | "Q7") -> "Q7":
        return value if isinstance(value, cls) else cls(Fraction(value), Fraction(0))

    def __add__(self, other: int | Fraction | "Q7") -> "Q7":
        right = self.of(other)
        return Q7(self.rational + right.rational, self.radical + right.radical)

    __radd__ = __add__

    def __neg__(self) -> "Q7":
        return Q7(-self.rational, -self.radical)

    def __sub__(self, other: int | Fraction | "Q7") -> "Q7":
        return self + (-self.of(other))

    def __rsub__(self, other: int | Fraction | "Q7") -> "Q7":
        return self.of(other) - self

    def __mul__(self, other: int | Fraction | "Q7") -> "Q7":
        right = self.of(other)
        return Q7(
            self.rational * right.rational + 7 * self.radical * right.radical,
            self.rational * right.radical + self.radical * right.rational,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | Fraction) -> "Q7":
        divisor = Fraction(scalar)
        require(divisor != 0, "division by zero")
        return Q7(self.rational / divisor, self.radical / divisor)

    def inverse(self) -> "Q7":
        norm = self.rational**2 - 7 * self.radical**2
        require(norm != 0, "inverse of zero Q(sqrt(7)) element")
        return Q7(self.rational / norm, -self.radical / norm)

    def json(self) -> list[list[int]]:
        return [
            [self.rational.numerator, self.rational.denominator],
            [self.radical.numerator, self.radical.denominator],
        ]


SQRT7 = Q7(Fraction(0), Fraction(1))


def validate_adjacency(adjacency: Sequence[int], order: int = 19) -> None:
    require(len(adjacency) == order, "graph order")
    full = (1 << order) - 1
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


def reconstruct_seed(
    adjacency: Sequence[int],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...], int]:
    require(is_clique(adjacency, SEED), "frozen seed is not a required K7")
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
        1 << local for local, defect in enumerate(defects)
        if defect.bit_count() >= 3
    )
    return outside, defects, tuple(ladj), eligible


def is_vertex_cover(ladj: Sequence[int], zmask: int) -> bool:
    complement = ((1 << len(ladj)) - 1) & ~zmask
    return all(not ladj[local] & complement for local in bits(complement))


def raw_covers(ladj: Sequence[int], eligible: int) -> tuple[int, ...]:
    return tuple(
        zmask
        for zmask in range(1 << len(ladj))
        if not zmask & ~eligible
        and zmask.bit_count() <= 7
        and is_vertex_cover(ladj, zmask)
    )


def propagated_by_local(zmask: int, witness: dict, outside_size: int) -> dict[int, int]:
    nonzero = tuple(local for local in range(outside_size) if not zmask & (1 << local))
    propagated = tuple(map(int, witness["propagated_masks"]))
    require(len(propagated) == len(nonzero), "propagated-mask/nonzero mismatch")
    return dict(zip(nonzero, propagated, strict=True))


def triple_petal_pattern_holds(
    adjacency: Sequence[int], outside: Sequence[int], masks: dict[int, int]
) -> bool:
    """Recognize exactly the local support/required-edge hypothesis."""

    if any(local not in masks for local in FIVE_LOCALS):
        return False
    if not is_clique(adjacency, tuple(outside[local] for local in FIVE_LOCALS)):
        return False
    petals = tuple(masks[local] for local in PETAL_LOCALS)
    hubs = tuple(masks[local] for local in HUB_LOCALS)
    if petals != EXPECTED_PETAL_MASKS or hubs != EXPECTED_HUB_MASKS:
        return False
    center = 1 << CENTER_COORDINATE
    leaves = tuple(1 << coordinate for coordinate in LEAF_COORDINATES)
    if any(petals[i] & petals[j] != center for i, j in combinations(range(3), 2)):
        return False
    if any(hub & petals[i] != leaves[i] for hub in hubs for i in range(3)):
        return False
    return hubs[0] & hubs[1] == sum(leaves)


def exact_algebra_certificate() -> dict:
    """Check both possible shared petal signs in exact Q(sqrt(7))."""

    cases = []
    for epsilon_integer in (-1, 1):
        epsilon = Q7.of(epsilon_integer)
        leaf = (SQRT7 - epsilon) / 2
        diagonal_factor = (epsilon - SQRT7) * (leaf - SQRT7)
        require(diagonal_factor == Q7.of(3), "two-support diagonal identity")
        leaf_square = leaf * leaf
        require(leaf_square != Q7(), "leaf value used by a required edge vanished")
        leaf_inverse = leaf.inverse()
        require(leaf * leaf_inverse == Q7.of(1), "leaf inverse identity")
        hub_petal_products = tuple(
            leaf * leaf_inverse for _hub in range(2) for _petal in range(3)
        )
        require(
            all(product == Q7.of(1) for product in hub_petal_products),
            "six hub-petal edge equations",
        )
        hub_dot = 3 * leaf_inverse * leaf_inverse
        cleared_hub_edge = leaf_square * (hub_dot - 1)
        require(
            cleared_hub_edge == 3 - leaf_square,
            "cleared hub-hub edge identity",
        )
        contradiction = leaf_square - 3
        require(contradiction != Q7(), "hub edge would be feasible")
        cases.append({
            "epsilon": epsilon_integer,
            "leaf_value": leaf.json(),
            "leaf_inverse": leaf_inverse.json(),
            "leaf_square": leaf_square.json(),
            "leaf_square_minus_three": contradiction.json(),
            "two_support_diagonal_factor": diagonal_factor.json(),
            "six_hub_petal_products": [
                product.json() for product in hub_petal_products
            ],
            "hub_hub_dot_forced_by_six_petal_edges": hub_dot.json(),
            "leaf_square_times_hub_dot_minus_one": cleared_hub_edge.json(),
            "three_minus_leaf_square": (3 - leaf_square).json(),
        })
    return {
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
        "cases": cases,
        "conclusion": (
            "lambda^2=2-epsilon*sqrt(7)/2 is never 3; equality would give "
            "epsilon*sqrt(7)=-2 and hence 7=4."
        ),
        "status": "EXACT_CONTRADICTION_FOR_BOTH_SIGNS",
    }


def lower_bound_18_graph() -> tuple[int, ...]:
    odd = tuple(value for value in range(32) if value.bit_count() % 2 == 1)
    adjacency = [0] * 18
    for first, code in enumerate(odd):
        for second in range(first):
            if (code ^ odd[second]).bit_count() == 2:
                adjacency[first] |= 1 << second
                adjacency[second] |= 1 << first
    for apex in (16, 17):
        for vertex in range(16):
            adjacency[apex] |= 1 << vertex
            adjacency[vertex] |= 1 << apex
    return tuple(adjacency)


def positive_and_synthetic_controls() -> dict:
    positive = lower_bound_18_graph()
    validate_adjacency(positive, order=18)
    k7_seeds = sum(
        is_clique(positive, chosen)
        for chosen in combinations(range(18), 7)
    )
    require(k7_seeds == 0, "known realizable 18-point K7 seed count")

    synthetic_masks = {
        0: EXPECTED_HUB_MASKS[0],
        1: EXPECTED_PETAL_MASKS[0],
        3: EXPECTED_PETAL_MASKS[1],
        6: EXPECTED_HUB_MASKS[1],
        9: EXPECTED_PETAL_MASKS[2],
    }
    # Relabel the five graph vertices back to the frozen local labels.
    expanded_outside = list(range(10))
    expanded_adjacency = [0] * 10
    for first, second in combinations(FIVE_LOCALS, 2):
        expanded_adjacency[first] |= 1 << second
        expanded_adjacency[second] |= 1 << first
    require(
        triple_petal_pattern_holds(
            expanded_adjacency, expanded_outside, synthetic_masks
        ),
        "synthetic exact pattern was not recognized",
    )
    missing_edge = list(expanded_adjacency)
    first, second = HUB_LOCALS
    missing_edge[first] &= ~(1 << second)
    missing_edge[second] &= ~(1 << first)
    require(
        not triple_petal_pattern_holds(missing_edge, expanded_outside, synthetic_masks),
        "optional hub nonedge was treated as a required edge",
    )
    return {
        "known_realizable_18": {"passed": True, "K7_seeds": k7_seeds},
        "synthetic_exact_pattern_detected": True,
        "synthetic_missing_hub_edge_not_rejected": True,
    }


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
    star = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_report.json").read_text(encoding="utf-8")
    )
    star_check = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v8"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "current manifest boundary",
    )
    require(
        manifest_check.get("schema")
        == "d6-current-certified-residue-v8-verification-v1"
        and manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v8.json"],
        "current manifest checker boundary",
    )
    require(
        star.get("kind")
        == "d6_k7_one_free_neighbour_two_free_center_increment"
        and star.get("status") == "COMPLETE"
        and star.get("semantics") == STAR_SEMANTICS,
        "star report boundary",
    )
    require(
        star_check.get("kind") == "d6_k7_one_two_star_increment_verification"
        and star_check.get("status") == "PASS"
        and star_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_one_two_star_increment_report.json"]
        and star_check.get("semantics") == STAR_SEMANTICS
        and star_check.get("checked", {}).get("eligible_covers") == 19_932
        and star_check.get("checked", {}).get("current_passing_families") == 88,
        "star independent checker boundary",
    )
    return manifest, star, star_check


def build_certificate(manifest: dict, star: dict) -> dict:
    indices = tuple(int(row["index"]) for row in manifest["classes"]["K7"]["graphs"])
    require(indices == EXPECTED_INPUT_INDICES, "current K7 residue list")
    require(stable_hash(list(indices)) == EXPECTED_INPUT_SHA256, "current K7 hash")
    graph = next(
        row for row in manifest["classes"]["K7"]["graphs"]
        if int(row["index"]) == TARGET_INDEX
    )
    adjacency = tuple(map(int, graph["adjacency"]))
    validate_adjacency(adjacency)
    require(stable_hash(list(adjacency)) == EXPECTED_ADJACENCY_SHA256, "adjacency hash")
    outside, defects, ladj, eligible = reconstruct_seed(adjacency)
    covers = raw_covers(ladj, eligible)
    require(len(covers) == EXPECTED_RAW_COVERS, "raw eligible cover count")
    require(stable_hash(list(covers)) == EXPECTED_RAW_COVERS_SHA256, "raw cover hash")

    graph_record = next(
        row for row in star["records"] if int(row["index"]) == TARGET_INDEX
    )
    require(graph_record.get("decision") == "SURVIVOR", "star target status")
    seed_record = next(
        row for row in graph_record["seeds"] if int(row["seed_mask"]) == SEED_MASK
    )
    require(tuple(seed_record["seed"]) == SEED, "star seed identity")
    require(seed_record.get("status") == "PASSING", "star seed status")
    archived = seed_record["current_covers"]
    require(
        tuple(int(row["zmask"]) for row in archived) == EXPECTED_CURRENT_COVERS,
        "current cover list",
    )
    cover_certificates = []
    for row in archived:
        zmask = int(row["zmask"])
        require(
            row.get("status") == "PASSING"
            and int(row.get("current_passing_families", -1)) == 1
            and int(row.get("star_passing_families", -1)) == 1
            and int(row.get("star_failed_families", -1)) == 0,
            f"cover multiplicity/status {zmask}",
        )
        witness = row["first_passing_witness"]
        require(
            tuple(map(int, witness["propagated_masks"]))
            == EXPECTED_PROPAGATED_BY_COVER[zmask],
            f"propagated masks {zmask}",
        )
        require(
            tuple(map(int, witness["z_supports"])) == EXPECTED_Z_SUPPORTS[zmask],
            f"Z supports {zmask}",
        )
        by_local = propagated_by_local(zmask, witness, len(outside))
        require(not any(zmask & (1 << local) for local in FIVE_LOCALS), "pattern vertex in Z")
        relevant = tuple(by_local[local] for local in FIVE_LOCALS)
        require(relevant == EXPECTED_RELEVANT_MASKS, f"relevant masks {zmask}")
        require(
            triple_petal_pattern_holds(adjacency, outside, by_local),
            f"triple-petal pattern {zmask}",
        )
        cover_certificates.append({
            "zmask": zmask,
            "Z_local_vertices": list(bits(zmask)),
            "Z_global_vertices": [outside[local] for local in bits(zmask)],
            "Z_actual_supports": list(map(int, witness["z_supports"])),
            "nonzero_local_vertices": [
                local for local in range(len(outside)) if not zmask & (1 << local)
            ],
            "propagated_support_supersets_in_nonzero_order": list(
                map(int, witness["propagated_masks"])
            ),
            "five_local_vertices": list(FIVE_LOCALS),
            "five_global_vertices": [outside[local] for local in FIVE_LOCALS],
            "five_propagated_masks": list(relevant),
            "all_ten_pairs_required_unit_edges": True,
            "contradiction": "triple_petal_Q_sqrt7",
        })

    return {
        "target": {
            "index": TARGET_INDEX,
            "adjacency_sha256": EXPECTED_ADJACENCY_SHA256,
            "seed": list(SEED),
            "seed_mask": SEED_MASK,
            "outside": list(outside),
            "defect_masks": list(defects),
        },
        "quantifier": {
            "raw_eligible_covers": len(covers),
            "raw_eligible_covers_sha256": stable_hash(list(covers)),
            "prior_layer_eliminated_raw_covers": len(covers) - len(archived),
            "upstream_current_covers": list(EXPECTED_CURRENT_COVERS),
            "current_passing_families_per_cover": 1,
            "new_families_checked": len(archived),
            "new_families_rejected": len(archived),
        },
        "pattern": {
            "hub_local_vertices": list(HUB_LOCALS),
            "petal_local_vertices": list(PETAL_LOCALS),
            "center_coordinate": CENTER_COORDINATE,
            "leaf_coordinates": list(LEAF_COORDINATES),
            "hub_masks": list(EXPECTED_HUB_MASKS),
            "petal_masks": list(EXPECTED_PETAL_MASKS),
            "five_pair_target": "required_K5_all_inner_products_one",
        },
        "covers": cover_certificates,
        "algebra": exact_algebra_certificate(),
        "graph_rejection_quantifier": (
            "Every realization must realize this required K7 seed. The "
            "hash-pinned parent exhausts all 255 eligible covers, eliminates "
            "253, and leaves exactly the two one-family covers checked here; "
            "both exact support/value systems are contradictory."
        ),
    }


def git(*arguments: str, binary: bool = False):
    return subprocess.run(
        ["git", *arguments], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=not binary,
    ).stdout


def committed_source_boundary() -> dict:
    commit = git("rev-parse", "HEAD").strip()
    branch = git("branch", "--show-current").strip()
    require(branch == "codex/dimension6", "wrong production branch")
    sources = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    for name, expected in sources.items():
        blob = git("show", f"{commit}:{name}", binary=True)
        require(hashlib.sha256(blob).hexdigest() == expected, f"uncommitted source: {name}")
    porcelain = git("status", "--porcelain=v1", "--untracked-files=all").strip()
    lines = porcelain.splitlines() if porcelain else []
    require(all(line.startswith("?? ") for line in lines), "tracked/staged dirt at launch")
    return {
        "commit": commit,
        "branch": branch,
        "source_sha256": sources,
        "porcelain_lines": lines,
        "porcelain_sha256": hashlib.sha256(porcelain.encode("utf-8")).hexdigest(),
        "proof_and_checker_sources_equal_committed_blobs": True,
        "tracked_clean": True,
    }


def build_report(source_boundary: dict | None = None) -> dict:
    manifest, star, star_check = load_upstream()
    certificate = build_certificate(manifest, star)
    controls = positive_and_synthetic_controls()
    if source_boundary is None:
        source_boundary = committed_source_boundary()
    source_hashes = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    require(source_boundary.get("source_sha256") == source_hashes, "source boundary hashes")
    output = [index for index in EXPECTED_INPUT_INDICES if index != TARGET_INDEX]
    require(stable_hash(output) == EXPECTED_OUTPUT_SHA256, "output residue hash")
    return {
        "schema": 1,
        "kind": REPORT_KIND,
        "status": "COMPLETE_EXACT_REJECTION",
        "claim": (
            "Graph 3936435 has no realization by distinct points in R^6 with "
            "all candidate edges at unit distance."
        ),
        "semantics": SEMANTICS,
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "source_sha256": source_hashes,
        "input": {
            "class": "current_K7_residue_v8",
            "ordered_indices": list(EXPECTED_INPUT_INDICES),
            "ordered_indices_sha256": EXPECTED_INPUT_SHA256,
            "graphs": len(EXPECTED_INPUT_INDICES),
        },
        "summary": {
            "graphs_rejected": 1,
            "rejected_indices": [TARGET_INDEX],
            "rejected_indices_sha256": stable_hash([TARGET_INDEX]),
            "graphs_surviving": len(output),
            "ordered_survivor_indices": output,
            "ordered_survivor_indices_sha256": EXPECTED_OUTPUT_SHA256,
        },
        "certificate": certificate,
        "certificate_sha256": stable_hash(certificate),
        "controls": controls,
        "parent_independent_check_counts": {
            "eligible_covers": star_check["checked"]["eligible_covers"],
            "current_passing_families": star_check["checked"][
                "current_passing_families"
            ],
            "positive_18_control": star_check["checked"]["positive_18_control"],
        },
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "finished_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "git": source_boundary,
        },
        "nonclaims": [
            "surviving graphs are not asserted realizable",
            "candidate nonedges are not asserted nonunit",
            "this increment does not settle the remaining K7 or K6-only residue",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPORT)
    args = parser.parse_args()
    report = build_report()
    atomic_json(args.output.resolve(), report)
    print(json.dumps({
        "status": report["status"],
        "rejected": report["summary"]["rejected_indices"],
        "surviving": report["summary"]["graphs_surviving"],
        "source_commit": report["execution"]["git"]["commit"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
