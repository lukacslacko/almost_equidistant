#!/usr/bin/env python3
"""Build the exact five-petal/basis increment for K7 graph 3936310.

The package begins at the independently verified graph-3936435 increment.
For one required K7 seed of graph 3936310, the pinned one/two-star boundary
leaves two covers and one propagated support family on each.  Five petal
vectors, five complementary basis vectors, and one singleton vector give an
exact contradiction in Q(sqrt(7)).

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
TARGET_INDEX = 3_936_310
SEED = (0, 3, 4, 9, 10, 14, 17)
SEED_MASK = 149_017
PETAL_LOCALS = (1, 2, 4, 7, 9)
PETAL_MASKS = (10, 6, 3, 18, 34)
CENTER_COORDINATE = 1
LEAF_COORDINATES = (3, 2, 0, 4, 5)
# Basis i omits the leaf of petal i and contains all four other leaves plus 6.
BASIS_LOCALS = (5, 3, 0, 6, 8)
BASIS_MASKS = (117, 121, 124, 109, 93)
SINGLETON_LOCAL = 11
SINGLETON_MASK = 64
SINGLETON_COORDINATE = 6
USED_LOCALS = PETAL_LOCALS + BASIS_LOCALS + (SINGLETON_LOCAL,)
EXPECTED_CURRENT_COVERS = (0, 1024)
EXPECTED_Z_SUPPORTS = {0: (), 1024: (63,)}
EXPECTED_PROPAGATED_BY_COVER = {
    0: (124, 10, 6, 121, 3, 117, 109, 18, 93, 34, 127, 64),
    1024: (124, 10, 6, 121, 3, 117, 109, 18, 93, 34, 64),
}
EXPECTED_RAW_COVERS = 64
EXPECTED_RAW_COVERS_SHA256 = (
    "b023ca3ca660d6b64c26c6242230717e4c71c1fadf415ff6ff21103fa78475eb"
)
EXPECTED_ADJACENCY_SHA256 = (
    "a7e155482e1d8a05365a61ecb08f751d491f196314ba14c9bd76c44ad0a3f17c"
)
EXPECTED_INPUT_INDICES = (
    316173,
    2581209,
    3648882,
    3729907,
    3935560,
    3936310,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
)
EXPECTED_INPUT_SHA256 = (
    "380936282d98f2e561c71680a04d04401fabc032e8a00cebf12ae4d18967a874"
)
EXPECTED_OUTPUT_SHA256 = (
    "609155c31c237e6e59c2abc0f4ddeb635fcda621adffa578932f9dde76b1936b"
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
    "d6_k7_triple_petals_3936435_increment_report.json": (
        "b13240a303f3e304dd2d095cc6dd28f3caa1f6571e0ba210637cb7e0ab0f168a"
    ),
    "d6_k7_triple_petals_3936435_increment_verification.json": (
        "d07916b38ec2e640d06eaf44fe351e92d32b32f522455d5b5ff008d233ad7c0e"
    ),
}
PACKAGE_SOURCES = (
    "build_d6_k7_five_petals_3936310_increment.py",
    "verify_d6_k7_five_petals_3936310_increment.py",
    "test_d6_k7_five_petals_3936310_increment.py",
    "d6_k7_five_petals_3936310_increment.md",
)
REPORT = ROOT / "d6_k7_five_petals_3936310_increment_report.json"
REPORT_KIND = "d6_k7_five_petals_basis_3936310_exact_increment"
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


def bits(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(mask.bit_length()) if mask & (1 << index))


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

    def conjugate(self) -> "Q7":
        return Q7(self.rational, -self.radical)

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
    lorentz_adjacency = [0] * len(outside)
    for first, second in combinations(range(len(outside)), 2):
        if (
            adjacency[outside[first]] & (1 << outside[second])
            and not defects[first] & defects[second]
        ):
            lorentz_adjacency[first] |= 1 << second
            lorentz_adjacency[second] |= 1 << first
    eligible = sum(
        1 << local
        for local, defect in enumerate(defects)
        if defect.bit_count() >= 3
    )
    return outside, defects, tuple(lorentz_adjacency), eligible


def raw_covers(lorentz_adjacency: Sequence[int], eligible: int) -> tuple[int, ...]:
    result = []
    full = (1 << len(lorentz_adjacency)) - 1
    for zmask in range(1 << len(lorentz_adjacency)):
        if zmask & ~eligible or zmask.bit_count() > 7:
            continue
        complement = full & ~zmask
        if all(not lorentz_adjacency[local] & complement for local in bits(complement)):
            result.append(zmask)
    return tuple(result)


def propagated_by_local(zmask: int, witness: dict, outside_size: int) -> dict[int, int]:
    nonzero = tuple(local for local in range(outside_size) if not zmask & (1 << local))
    propagated = tuple(map(int, witness["propagated_masks"]))
    require(len(propagated) == len(nonzero), "propagated-mask/nonzero mismatch")
    return dict(zip(nonzero, propagated, strict=True))


def five_petal_basis_pattern_holds(
    adjacency: Sequence[int], outside: Sequence[int], masks: dict[int, int]
) -> bool:
    """Recognize exactly the support and required-edge hypotheses."""

    if any(local not in masks for local in USED_LOCALS):
        return False
    petals = tuple(masks[local] for local in PETAL_LOCALS)
    bases = tuple(masks[local] for local in BASIS_LOCALS)
    if petals != PETAL_MASKS or bases != BASIS_MASKS:
        return False
    if masks[SINGLETON_LOCAL] != SINGLETON_MASK:
        return False
    petal_globals = tuple(outside[local] for local in PETAL_LOCALS)
    basis_globals = tuple(outside[local] for local in BASIS_LOCALS)
    singleton_global = outside[SINGLETON_LOCAL]
    if not is_clique(adjacency, petal_globals):
        return False
    if not is_clique(adjacency, basis_globals):
        return False
    if any(
        not adjacency[singleton_global] & (1 << basis)
        for basis in basis_globals
    ):
        return False
    if any(
        petal_number != basis_number
        and not adjacency[basis_globals[basis_number]]
        & (1 << petal_globals[petal_number])
        for basis_number in range(5)
        for petal_number in range(5)
    ):
        return False
    center_bit = 1 << CENTER_COORDINATE
    leaf_bits = tuple(1 << coordinate for coordinate in LEAF_COORDINATES)
    singleton_bit = 1 << SINGLETON_COORDINATE
    if any(
        petals[first] & petals[second] != center_bit
        for first, second in combinations(range(5), 2)
    ):
        return False
    for basis_number, basis_mask in enumerate(bases):
        if basis_mask & petals[basis_number]:
            return False
        for petal_number in range(5):
            if petal_number != basis_number:
                if basis_mask & petals[petal_number] != leaf_bits[petal_number]:
                    return False
    if any(
        bases[first] & bases[second]
        != singleton_bit
        + sum(
            leaf_bits[number]
            for number in range(5)
            if number not in (first, second)
        )
        for first, second in combinations(range(5), 2)
    ):
        return False
    return all(mask & SINGLETON_MASK == singleton_bit for mask in bases)


def exact_algebra_certificate() -> dict:
    """Check both common petal-center signs in exact Q(sqrt(7))."""

    singleton = 4 * SQRT7 / 7
    singleton_diagonal = singleton * singleton - 1 - (singleton - SQRT7) * (
        singleton - SQRT7
    )
    require(singleton_diagonal == Q7(), "singleton diagonal identity")
    singleton_inverse = singleton.inverse()
    require(singleton * singleton_inverse == Q7.of(1), "singleton inverse")
    require(singleton_inverse == SQRT7 / 4, "singleton forced basis coordinate")

    cases = []
    for epsilon_integer in (-1, 1):
        epsilon = Q7.of(epsilon_integer)
        leaf = (SQRT7 - epsilon) / 2
        diagonal_factor = (epsilon - SQRT7) * (leaf - SQRT7)
        require(diagonal_factor == Q7.of(3), "two-support diagonal identity")
        require(leaf != Q7(), "petal leaf value vanished")
        leaf_inverse = leaf.inverse()
        require(leaf * leaf_inverse == Q7.of(1), "petal leaf inverse")
        petal_products = tuple(leaf * leaf_inverse for _ in range(20))
        singleton_products = tuple(
            singleton * singleton_inverse for _ in range(5)
        )
        require(all(value == Q7.of(1) for value in petal_products), "20 cross edges")
        require(
            all(value == Q7.of(1) for value in singleton_products),
            "five singleton-basis edges",
        )
        leaf_square = leaf * leaf
        forced_basis_dot = 3 * leaf_inverse * leaf_inverse + (
            singleton_inverse * singleton_inverse
        )
        cleared_basis_edge = leaf_square * (forced_basis_dot - 1)
        expected_cleared = 3 - Fraction(9, 16) * leaf_square
        require(cleared_basis_edge == expected_cleared, "cleared basis-edge identity")
        contradiction = 3 * leaf_square - 16
        require(
            cleared_basis_edge == Fraction(-3, 16) * contradiction,
            "basis edge equivalent to 3 lambda^2=16",
        )
        require(contradiction != Q7(), "basis edge would be feasible")
        norm = contradiction * contradiction.conjugate()
        require(norm == Q7(Fraction(337, 4)), "contradiction norm")
        cases.append({
            "epsilon": epsilon_integer,
            "petal_leaf": leaf.json(),
            "petal_leaf_inverse": leaf_inverse.json(),
            "petal_leaf_square": leaf_square.json(),
            "petal_two_support_diagonal_factor": diagonal_factor.json(),
            "singleton_value": singleton.json(),
            "singleton_inverse": singleton_inverse.json(),
            "singleton_diagonal_residual": singleton_diagonal.json(),
            "twenty_basis_petal_products": [value.json() for value in petal_products],
            "five_singleton_basis_products": [
                value.json() for value in singleton_products
            ],
            "basis_pair_dot_forced": forced_basis_dot.json(),
            "leaf_square_times_basis_dot_minus_one": cleared_basis_edge.json(),
            "three_minus_nine_sixteenths_leaf_square": expected_cleared.json(),
            "three_leaf_square_minus_sixteen": contradiction.json(),
            "contradiction_conjugate_product": norm.json(),
        })
    return {
        "field": "Q(sqrt(7)) with positive sqrt(7)",
        "petal_center_deduction": (
            "The required petal K5 gives x_i*x_j=1 for every pair, so all "
            "five center values equal epsilon in {-1,+1}."
        ),
        "petal_leaf_deduction": (
            "The normalized diagonal identity uniquely forces every distinct "
            "petal leaf to lambda=(sqrt(7)-epsilon)/2."
        ),
        "basis_coordinate_deduction": (
            "Twenty required basis-petal edges force 1/lambda on each of the "
            "four allowed leaves of every basis.  The singleton diagonal "
            "forces sigma=4/sqrt(7), and five singleton-basis edges force "
            "sqrt(7)/4 on basis coordinate 6."
        ),
        "basis_edge_deduction": (
            "Each required basis-basis edge has three shared leaves and shared "
            "coordinate 6, so 3/lambda^2+7/16=1 and lambda^2=16/3."
        ),
        "cases": cases,
        "conclusion": (
            "lambda^2=2-epsilon*sqrt(7)/2 cannot equal 16/3; equality would "
            "give epsilon*sqrt(7)=-20/3 and hence 63=400."
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


def positive_controls() -> dict:
    positive = lower_bound_18_graph()
    validate_adjacency(positive, order=18)
    k7_seeds = sum(
        is_clique(positive, chosen) for chosen in combinations(range(18), 7)
    )
    require(k7_seeds == 0, "known realizable 18-point K7 seed count")
    return {"known_realizable_18": {"passed": True, "K7_seeds": k7_seeds}}


def load_upstream() -> tuple[dict, dict, dict, dict, dict]:
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
    prior = json.loads(
        (ROOT / "d6_k7_triple_petals_3936435_increment_report.json").read_text(
            encoding="utf-8"
        )
    )
    prior_check = json.loads(
        (ROOT / "d6_k7_triple_petals_3936435_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v8"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "current manifest boundary",
    )
    require(
        manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v8.json"],
        "manifest checker boundary",
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
    require(
        prior.get("kind") == "d6_k7_triple_petals_3936435_exact_increment"
        and prior.get("status") == "COMPLETE_EXACT_REJECTION"
        and prior.get("summary", {}).get("ordered_survivor_indices")
        == list(EXPECTED_INPUT_INDICES)
        and prior.get("summary", {}).get("ordered_survivor_indices_sha256")
        == EXPECTED_INPUT_SHA256,
        "prior exact increment boundary",
    )
    require(
        prior_check.get("kind")
        == "d6_k7_triple_petals_3936435_increment_verification"
        and prior_check.get("status") == "PASS"
        and prior_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_triple_petals_3936435_increment_report.json"]
        and prior_check.get("conclusion", {}).get("rejected_indices") == [3936435],
        "prior independent verification boundary",
    )
    return manifest, star, star_check, prior, prior_check


def build_certificate(manifest: dict, star: dict) -> dict:
    graph = next(
        row
        for row in manifest["classes"]["K7"]["graphs"]
        if int(row["index"]) == TARGET_INDEX
    )
    adjacency = tuple(map(int, graph["adjacency"]))
    validate_adjacency(adjacency)
    require(stable_hash(list(adjacency)) == EXPECTED_ADJACENCY_SHA256, "adjacency hash")
    outside, defects, lorentz_adjacency, eligible = reconstruct_seed(adjacency)
    covers = raw_covers(lorentz_adjacency, eligible)
    require(len(covers) == EXPECTED_RAW_COVERS, "raw eligible cover count")
    require(stable_hash(list(covers)) == EXPECTED_RAW_COVERS_SHA256, "raw cover hash")

    graph_record = next(row for row in star["records"] if int(row["index"]) == TARGET_INDEX)
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
        require(not any(zmask & (1 << local) for local in USED_LOCALS), "pattern enters Z")
        require(
            five_petal_basis_pattern_holds(adjacency, outside, by_local),
            f"five-petal/basis pattern {zmask}",
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
            "petal_local_vertices": list(PETAL_LOCALS),
            "petal_global_vertices": [outside[local] for local in PETAL_LOCALS],
            "petal_masks": list(PETAL_MASKS),
            "basis_local_vertices_aligned_to_petals": list(BASIS_LOCALS),
            "basis_global_vertices_aligned_to_petals": [
                outside[local] for local in BASIS_LOCALS
            ],
            "basis_masks_aligned_to_petals": list(BASIS_MASKS),
            "singleton_local_vertex": SINGLETON_LOCAL,
            "singleton_global_vertex": outside[SINGLETON_LOCAL],
            "singleton_mask": SINGLETON_MASK,
            "required_edges_used": {
                "petal_K5": 10,
                "basis_K5": 10,
                "off_diagonal_basis_petal": 20,
                "singleton_basis": 5,
            },
            "contradiction": "five_petal_complementary_basis_Q_sqrt7",
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
            "petal_local_vertices": list(PETAL_LOCALS),
            "petal_masks": list(PETAL_MASKS),
            "center_coordinate": CENTER_COORDINATE,
            "leaf_coordinates_aligned_to_petals": list(LEAF_COORDINATES),
            "basis_local_vertices_aligned_to_petals": list(BASIS_LOCALS),
            "basis_masks_aligned_to_petals": list(BASIS_MASKS),
            "singleton_local_vertex": SINGLETON_LOCAL,
            "singleton_coordinate": SINGLETON_COORDINATE,
        },
        "covers": cover_certificates,
        "algebra": exact_algebra_certificate(),
        "graph_rejection_quantifier": (
            "Every realization must realize this required K7 seed. The "
            "hash-pinned parent exhausts all 64 eligible covers, eliminates "
            "62, and leaves exactly the two one-family covers checked here; "
            "both exact support/value systems are contradictory."
        ),
    }


def git(*arguments: str, binary: bool = False):
    return subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=not binary,
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
    manifest, star, star_check, prior, prior_check = load_upstream()
    certificate = build_certificate(manifest, star)
    controls = positive_controls()
    if source_boundary is None:
        source_boundary = committed_source_boundary()
    source_hashes = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    require(source_boundary.get("source_sha256") == source_hashes, "source boundary hashes")
    require(
        prior["summary"]["ordered_survivor_indices"] == list(EXPECTED_INPUT_INDICES),
        "input inherited from prior increment",
    )
    output = [index for index in EXPECTED_INPUT_INDICES if index != TARGET_INDEX]
    require(stable_hash(output) == EXPECTED_OUTPUT_SHA256, "output residue hash")
    return {
        "schema": 1,
        "kind": REPORT_KIND,
        "status": "COMPLETE_EXACT_REJECTION",
        "claim": (
            "Graph 3936310 has no realization by distinct points in R^6 with "
            "all candidate edges at unit distance."
        ),
        "semantics": SEMANTICS,
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "source_sha256": source_hashes,
        "input": {
            "class": "K7_residue_after_graph_3936435_increment",
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
            "prior_exact_rejected_indices": prior_check["conclusion"][
                "rejected_indices"
            ],
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
