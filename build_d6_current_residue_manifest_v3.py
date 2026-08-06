#!/usr/bin/env python3
"""Build the self-contained exact dimension-six residue manifest v3.

Starting from the 258 K7 and 977 K6-only adjacencies embedded in the committed
v2 manifest, this package applies the independently verified K7
support-or-pentad result and the complete three-stage K6 PSD-Z/Hall chain.  It
embeds the ordered 155 + 822 = 977 surviving adjacency records.  Membership is
exact set arithmetic; survival is not a realizability claim.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import time
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
K7_RESULT_COMMIT = "01c5b18413b6bab5225558fa950971efeda58e0d"
K6_RESULT_COMMIT = "a21db74bacdf4c9c5c841ac137af326b331d8137"
K7_SOURCE_COMMIT = "4982e564eecf552f83160644a61fb78662b13e96"
K6_SOURCE_COMMIT = "4e1d50492ad76867b67d748792f21ca5c0a6a94f"

EXPECTED_FILES = {
    "build_d6_current_residue_manifest_v2.py": (
        "f3db5ce7bbca61f432a73ca56184528504184ddce2d2c2cacf9e70db4d7ad5f1"
    ),
    "verify_d6_current_residue_manifest_v2.py": (
        "f015d4ee199ec928655a21a2fa7bd09aa85c70ad02fd14d755d2678c53ecd64c"
    ),
    "d6_current_residue_manifest_v2.json": (
        "961dac1f9b44bb541e2c5bd1626027826eeea7bc628bf5ff0a85ab26afa949f4"
    ),
    "d6_current_residue_manifest_v2_verification.json": (
        "35cb7a37c759d6626e34b4fa1174e8cc2d78e5667cee37be3f61af35227f80b3"
    ),
    "build_d6_k7_pentad_joint_conjunction.py": (
        "383cc8d55298e5d82f974e2b17d79bb4c830d4316151b1d34bbb2d67c5c81768"
    ),
    "verify_d6_k7_pentad_joint_conjunction.py": (
        "b5ee8c3c3b75d59cabbb43da79770d4dbc10196d8aa55f0ed5f210c5b077f217"
    ),
    "test_d6_k7_pentad_joint_conjunction.py": (
        "87016da50967397a46c88e90ee433fe8da553e8f928d147d4615ef6de3d13cc1"
    ),
    "d6_k7_pentad_joint_conjunction_report.json": (
        "4528b03a4c8c7b644051e6349dbd10dc3a634f81b35c56d46b09c471eb0436fa"
    ),
    "d6_k7_pentad_joint_conjunction_verification.json": (
        "5c47730967dc7a9a41fb7810642ed354d50ddd30ded59bca56025cc58b92acb7"
    ),
    "d6_k7_pentad_joint_conjunction_result_manifest.json": (
        "bea52786aa6100d15a37c44354d873e74a5ad7cd2eea49506d2417d367d44866"
    ),
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "verify_d6_k6_psd_zmatrix.py": (
        "ca4ca62f75d02846f12e59e16eb9eab4473fe5438fc496e72a25a39db77f9042"
    ),
    "test_d6_k6_psd_zmatrix.py": (
        "1d8662efcc8e091c166aac3383fdbbf154eeeac6d5a7f6ceeb90d02af6a072f2"
    ),
    "d6_k6_psd_zmatrix_report.json": (
        "abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30"
    ),
    "d6_k6_psd_zmatrix_verification.json": (
        "8261b58e788f120e3d1eee6dc481a4989e3eb75f9f03babbd552474b4a5302eb"
    ),
    "d6_k6_psd_zmatrix_certificates.json": (
        "62fc7170d97056941bcaf553c45cfa1cabdf20117ffc69fbd8429b3764ec999c"
    ),
    "d6_k6_psd_z_hereditary.py": (
        "7b869895cb9b4f7f3bbf1a1d1b11ec71507b6355c8ec0b7d7f8eaec2a51b180d"
    ),
    "verify_d6_k6_psd_z_hereditary.py": (
        "297e749923ae60b84cda12170d19509da237b665965d01af6214465f4e4df17f"
    ),
    "test_d6_k6_psd_z_hereditary.py": (
        "8206708620b323c418f3dd0b1cf4d35854bcf15cfaaddc545954a4233079294d"
    ),
    "d6_k6_psd_z_hereditary_report.json": (
        "202a844d6505d3c68c9a0bea8a5d82a983a7e44b96cb9f3c711d80c47f6c00c1"
    ),
    "d6_k6_psd_z_hereditary_verification.json": (
        "e9f159863e7605d0ef75485dad77a045e418c87687e0c4b7138d710b071924be"
    ),
    "d6_k6_psd_z_hereditary_certificates.json": (
        "799602d9516bf44f29594c5ded0a32836a32d64bba82129754d865fa323810ca"
    ),
    "d6_k6_psd_z_hereditary_checkpoint.json": (
        "1c838e0fff3cedfc2deac7bc6ed37c801a98fc48fe43cc651dba2350ba34ff42"
    ),
    "d6_k6_arbitrary_subset_hall.py": (
        "cadc0048fcc7a77299bc2950e2eca57bfed87d54ed6ffd59886d155b2c9a88ff"
    ),
    "verify_d6_k6_arbitrary_subset_hall.py": (
        "4ad4f886d67b3d1ebf149277050a639cda243f657ef02b553e6a3fa241f2853e"
    ),
    "test_d6_k6_arbitrary_subset_hall.py": (
        "e92896a2d426d9f98a03e210990d76ca314ddfa0f15a56a765f516bfb719cde1"
    ),
    "d6_k6_arbitrary_subset_hall_report.json": (
        "bcc4baecf28a95443a5b13e87a796cbba0a97cd42fa76573119bfd81ef3735fc"
    ),
    "d6_k6_arbitrary_subset_hall_verification.json": (
        "e6042edfff3683043903f9be544322379e85673207e35ee87349637ff4c73058"
    ),
    "d6_k6_arbitrary_subset_hall_certificates.json": (
        "b7fb04f769ddaf5d36630372f5134f2e4a32e572f620b26059667e9cf65f4a23"
    ),
    "d6_k6_arbitrary_subset_hall_checkpoint.json": (
        "43c1b65da336bc0db006313b97630867a660b944b06c814433ee3a8cdc8297e6"
    ),
    "d6_k6_arbitrary_subset_hall_result_manifest.json": (
        "a245fb4ae68d5c93b75e707799346e360358e78f17cc97c7b8f4ecb4705674bd"
    ),
}


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


def load(name: str) -> dict:
    value = json.loads((ROOT / name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{name} is not a JSON object")
    return value


def integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    if len(value) != len(set(value)):
        raise ValueError(f"{label} repeats an index")
    return list(value)


def verify_hashes() -> dict[str, str]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED_FILES}
    if observed != EXPECTED_FILES:
        raise ValueError("an upstream source/result hash changed")
    return observed


def validate_adjacency(rows: object, vertices: int, label: str) -> list[int]:
    if (
        not isinstance(rows, list) or len(rows) != vertices
        or any(type(row) is not int for row in rows)
    ):
        raise ValueError(f"{label} is not a {vertices}-row integer adjacency")
    full = (1 << vertices) - 1
    for vertex, row in enumerate(rows):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError(f"{label} has an invalid row {vertex}")
        neighbours_back = sum(
            1 << other for other in range(vertices)
            if rows[other] & (1 << vertex)
        )
        if row != neighbours_back:
            raise ValueError(f"{label} is asymmetric at vertex {vertex}")
    return list(rows)


def alpha_at_most_two(rows: Sequence[int]) -> bool:
    full = (1 << len(rows)) - 1
    for first, row in enumerate(rows):
        nonneighbours = full & ~row & ~(1 << first)
        remaining = nonneighbours & ~((1 << (first + 1)) - 1)
        while remaining:
            bit = remaining & -remaining
            second = bit.bit_length() - 1
            third = nonneighbours & ~rows[second] & ~(1 << second)
            if third:
                return False
            remaining ^= bit
    return True


def contains_clique(rows: Sequence[int], target: int) -> bool:
    full = (1 << len(rows)) - 1

    def search(candidates: int, need: int) -> bool:
        if need == 0:
            return True
        if candidates.bit_count() < need:
            return False
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            if search(candidates & rows[vertex], need - 1):
                return True
            if candidates.bit_count() < need:
                return False
        return False

    return search(full, target)


def count_cliques(rows: Sequence[int], size: int) -> int:
    return sum(
        all(rows[first] & (1 << second) for first, second in itertools.combinations(group, 2))
        for group in itertools.combinations(range(len(rows)), size)
    )


def positive_18_graph() -> list[int]:
    base = [word for word in range(32) if word.bit_count() % 2 == 1]
    rows = [0] * 18
    for first, word in enumerate(base):
        for second, other in enumerate(base[:first]):
            if (word ^ other).bit_count() == 2:
                rows[first] |= 1 << second
                rows[second] |= 1 << first
    for apex in (16, 17):
        for vertex in range(16):
            rows[apex] |= 1 << vertex
            rows[vertex] |= 1 << apex
    return rows


def graph_lookup(records: object, label: str) -> tuple[list[int], dict[int, list[int]]]:
    if not isinstance(records, list):
        raise ValueError(f"{label} records are not a list")
    order = []
    lookup = {}
    for position, record in enumerate(records):
        if not isinstance(record, dict) or type(record.get("index")) is not int:
            raise ValueError(f"{label} record {position} has no integer index")
        index = int(record["index"])
        if index in lookup:
            raise ValueError(f"{label} repeats graph {index}")
        order.append(index)
        lookup[index] = validate_adjacency(
            record.get("adjacency"), 19, f"{label} graph {index}"
        )
    return order, lookup


def ordered_records(indices: Sequence[int], lookup: dict[int, list[int]]) -> list[dict]:
    if any(index not in lookup for index in indices):
        raise ValueError("a residue index has no v2 adjacency")
    return [{"index": index, "adjacency": lookup[index]} for index in indices]


def check_report_stage(
    report: dict,
    verification: dict,
    input_indices: Sequence[int],
    report_schema: str,
    verification_schema: str,
    rejected_count: int,
    survivor_count: int,
) -> tuple[list[int], list[int]]:
    order = [
        int(record["index"]) for record in report.get("graph_results", [])
        if isinstance(record, dict) and type(record.get("index")) is int
    ]
    rejected = integer_list(report.get("rejected_indices"), "stage rejections")
    decisions = [
        int(record["index"]) for record in report.get("graph_results", [])
        if record.get("decision", {}).get("rejected") is True
    ]
    checked = integer_list(
        verification.get("rejected_indices"), "verified stage rejections"
    )
    if (
        report.get("schema") != report_schema
        or report.get("status") != "COMPLETE"
        or verification.get("schema") != verification_schema
        or verification.get("status") != "PASS"
        or order != list(input_indices)
        or report.get("input_graphs") != len(input_indices)
        or report.get("input_indices_sha256") != stable_hash(list(input_indices))
        or report.get("graphs_rejected") != rejected_count
        or report.get("graphs_surviving") != survivor_count
        or verification.get("graphs_recomputed") != len(input_indices)
        or verification.get("graphs_rejected") != rejected_count
        or verification.get("graphs_surviving") != survivor_count
        or rejected != decisions or rejected != checked
        or len(rejected) != rejected_count
        or not set(rejected) <= set(input_indices)
        or report.get("positive_control", {}).get("passed") is not True
        or verification.get("positive_control", {}).get("passed") is not True
    ):
        raise ValueError(f"stage gate failed for {report_schema}")
    survivors = [index for index in input_indices if index not in set(rejected)]
    if len(survivors) != survivor_count:
        raise ValueError("stage ordered complement count failed")
    return rejected, survivors


def build_manifest() -> dict:
    observed = verify_hashes()
    v2 = load("d6_current_residue_manifest_v2.json")
    v2_verification = load("d6_current_residue_manifest_v2_verification.json")
    if (
        v2.get("schema") != "d6-current-exact-residue-v2"
        or v2.get("status") != "COMPLETE"
        or v2_verification.get("schema")
        != "d6-current-exact-residue-v2-verification-v1"
        or v2_verification.get("status") != "PASS"
        or v2_verification.get("manifest", {}).get("sha256")
        != observed["d6_current_residue_manifest_v2.json"]
        or v2.get("class_order") != ["K7", "K6_only"]
    ):
        raise ValueError("v2 exact boundary gate failed")
    k7_v2 = v2["classes"]["K7"]
    k6_v2 = v2["classes"]["K6_only"]
    k7_base = integer_list(k7_v2.get("residue_indices"), "v2 K7 residue")
    k6_base = integer_list(k6_v2.get("residue_indices"), "v2 K6-only residue")
    k7_order, k7_lookup = graph_lookup(k7_v2.get("graphs"), "v2 K7")
    k6_order, k6_lookup = graph_lookup(k6_v2.get("graphs"), "v2 K6-only")
    if (
        k7_base != k7_order or len(k7_base) != 258
        or stable_hash(k7_base) != k7_v2.get("residue_indices_sha256")
        or stable_hash(k7_v2["graphs"]) != k7_v2.get("graphs_sha256")
        or k6_base != k6_order or len(k6_base) != 977
        or stable_hash(k6_base) != k6_v2.get("residue_indices_sha256")
        or stable_hash(k6_v2["graphs"]) != k6_v2.get("graphs_sha256")
    ):
        raise ValueError("v2 embedded class order/hash failed")

    # K7 258 -> 155.
    k7_report = load("d6_k7_pentad_joint_conjunction_report.json")
    k7_check = load("d6_k7_pentad_joint_conjunction_verification.json")
    k7_result = load("d6_k7_pentad_joint_conjunction_result_manifest.json")
    joint = integer_list(
        k7_report.get("joint_support_rejections", {}).get("ordered_indices"),
        "K7 joint-support rejections",
    )
    incremental = integer_list(
        k7_report.get("new_seed_conjunction_rejections", {}).get("ordered_indices"),
        "K7 pentad increment",
    )
    k7_rejected = integer_list(
        k7_report.get("combined_rejections", {}).get("ordered_indices"),
        "K7 combined rejections",
    )
    k7_residue = integer_list(
        k7_report.get("exact_residue", {}).get("ordered_indices"),
        "K7 exact residue",
    )
    joint_residue = [index for index in k7_base if index not in set(joint)]
    if (
        k7_report.get("kind") != "d6_k7_pentad_joint_seed_conjunction"
        or k7_report.get("status") != "COMPLETE"
        or k7_check.get("kind")
        != "d6_k7_pentad_joint_seed_conjunction_verification"
        or k7_check.get("status") != "PASS"
        or k7_check.get("report", {}).get("sha256")
        != observed["d6_k7_pentad_joint_conjunction_report.json"]
        or k7_report.get("selection", {}).get("graphs") != len(k7_base)
        or k7_report.get("selection", {}).get("ordered_indices_sha256")
        != stable_hash(k7_base)
        or len(joint) != 69 or len(incremental) != 34
        or set(joint) & set(incremental)
        or not set(incremental) <= set(joint_residue)
        or k7_rejected != [index for index in k7_base if index in set(joint) | set(incremental)]
        or k7_residue != [index for index in k7_base if index not in set(k7_rejected)]
        or k7_residue != [
            index for index in joint_residue if index not in set(incremental)
        ]
        or len(k7_rejected) != 103 or len(k7_residue) != 155
        or k7_report["combined_rejections"]["ordered_indices_sha256"]
        != stable_hash(k7_rejected)
        or k7_report["exact_residue"]["ordered_indices_sha256"]
        != stable_hash(k7_residue)
        or k7_check.get("ordered_rejections_sha256") != stable_hash(k7_rejected)
        or k7_check.get("ordered_residue_sha256") != stable_hash(k7_residue)
    ):
        raise ValueError("K7 exact 258-to-155 gate failed")
    if (
        k7_result.get("kind")
        != "d6-k7-support-pentad-conjunction-result-manifest"
        or k7_result.get("status") != "COMPLETE"
        or k7_result.get("source_boundary", {}).get("commit") != K7_SOURCE_COMMIT
        or k7_result.get("source_boundary", {}).get("builder_sha256")
        != observed["build_d6_k7_pentad_joint_conjunction.py"]
        or k7_result.get("source_boundary", {}).get("verifier_sha256")
        != observed["verify_d6_k7_pentad_joint_conjunction.py"]
        or k7_result.get("artifacts", {}).get("report", {}).get("sha256")
        != observed["d6_k7_pentad_joint_conjunction_report.json"]
        or k7_result.get("artifacts", {}).get("report", {}).get("bytes")
        != (ROOT / "d6_k7_pentad_joint_conjunction_report.json").stat().st_size
        or k7_result.get("artifacts", {}).get("verification", {}).get("sha256")
        != observed["d6_k7_pentad_joint_conjunction_verification.json"]
        or k7_result.get("artifacts", {}).get("verification", {}).get("bytes")
        != (ROOT / "d6_k7_pentad_joint_conjunction_verification.json").stat().st_size
        or k7_result.get("exact_result", {}).get("input_graphs") != 258
        or k7_result.get("exact_result", {}).get("combined_rejected_graphs") != 103
        or k7_result.get("exact_result", {}).get("combined_rejected_indices_sha256")
        != stable_hash(k7_rejected)
        or k7_result.get("exact_result", {}).get("residue_graphs") != 155
        or k7_result.get("exact_result", {}).get("residue_indices_sha256")
        != stable_hash(k7_residue)
    ):
        raise ValueError("K7 result-manifest gate failed")

    # K6-only 977 -> 861 -> 831 -> 822.
    psd_report = load("d6_k6_psd_zmatrix_report.json")
    psd_check = load("d6_k6_psd_zmatrix_verification.json")
    hereditary_report = load("d6_k6_psd_z_hereditary_report.json")
    hereditary_check = load("d6_k6_psd_z_hereditary_verification.json")
    arbitrary_report = load("d6_k6_arbitrary_subset_hall_report.json")
    arbitrary_check = load("d6_k6_arbitrary_subset_hall_verification.json")
    arbitrary_result = load("d6_k6_arbitrary_subset_hall_result_manifest.json")
    psd_rejected, psd_residue = check_report_stage(
        psd_report, psd_check, k6_base,
        "d6-k6-psd-zmatrix-v1", "d6-k6-psd-zmatrix-verification-v1",
        116, 861,
    )
    hereditary_rejected, hereditary_residue = check_report_stage(
        hereditary_report, hereditary_check, psd_residue,
        "d6-k6-psd-z-hereditary-v1",
        "d6-k6-psd-z-hereditary-verification-v1", 30, 831,
    )
    arbitrary_rejected, k6_residue = check_report_stage(
        arbitrary_report, arbitrary_check, hereditary_residue,
        "d6-k6-arbitrary-subset-hall-v1",
        "d6-k6-arbitrary-subset-hall-verification-v1", 9, 822,
    )
    for report, certificate_name, checkpoint_name in (
        (psd_report, "d6_k6_psd_zmatrix_certificates.json", None),
        (
            hereditary_report,
            "d6_k6_psd_z_hereditary_certificates.json",
            "d6_k6_psd_z_hereditary_checkpoint.json",
        ),
        (
            arbitrary_report,
            "d6_k6_arbitrary_subset_hall_certificates.json",
            "d6_k6_arbitrary_subset_hall_checkpoint.json",
        ),
    ):
        certificate = report.get("certificate_archive", {})
        if (
            certificate.get("path") != certificate_name
            or certificate.get("sha256") != observed[certificate_name]
            or certificate.get("bytes") != (ROOT / certificate_name).stat().st_size
        ):
            raise ValueError(f"certificate archive binding failed: {certificate_name}")
        if checkpoint_name is not None:
            checkpoint = report.get("checkpoint", {})
            if (
                checkpoint.get("path") != checkpoint_name
                or checkpoint.get("sha256") != observed[checkpoint_name]
            ):
                raise ValueError(f"checkpoint binding failed: {checkpoint_name}")
    if (
        arbitrary_result.get("kind")
        != "d6-k6-arbitrary-subset-hall-result-manifest"
        or arbitrary_result.get("status") != "COMPLETE"
        or arbitrary_result.get("source_boundary", {}).get("commit")
        != K6_SOURCE_COMMIT
        or arbitrary_result.get("source_boundary", {}).get("production_sha256")
        != observed["d6_k6_arbitrary_subset_hall.py"]
        or arbitrary_result.get("source_boundary", {}).get("verifier_sha256")
        != observed["verify_d6_k6_arbitrary_subset_hall.py"]
        or arbitrary_result.get("parent", {}).get("sha256")
        != observed["d6_k6_psd_z_hereditary_report.json"]
        or arbitrary_result.get("artifacts", {}).get("report", {}).get("sha256")
        != observed["d6_k6_arbitrary_subset_hall_report.json"]
        or arbitrary_result.get("artifacts", {}).get("verification", {}).get("sha256")
        != observed["d6_k6_arbitrary_subset_hall_verification.json"]
        or arbitrary_result.get("artifacts", {}).get("certificates", {}).get("sha256")
        != observed["d6_k6_arbitrary_subset_hall_certificates.json"]
        or arbitrary_result.get("artifacts", {}).get("checkpoint", {}).get("sha256")
        != observed["d6_k6_arbitrary_subset_hall_checkpoint.json"]
        or arbitrary_result.get("exact_result", {}).get("input_graphs") != 831
        or arbitrary_result.get("exact_result", {}).get("rejected_graphs") != 9
        or arbitrary_result.get("exact_result", {}).get("rejected_indices")
        != arbitrary_rejected
        or arbitrary_result.get("exact_result", {}).get("rejected_indices_sha256")
        != stable_hash(arbitrary_rejected)
        or arbitrary_result.get("exact_result", {}).get("surviving_graphs") != 822
    ):
        raise ValueError("K6 arbitrary-subset result-manifest gate failed")
    k6_stage_sets = [set(psd_rejected), set(hereditary_rejected), set(arbitrary_rejected)]
    if any(k6_stage_sets[i] & k6_stage_sets[j] for i in range(3) for j in range(i)):
        raise ValueError("K6 stage rejection sets overlap")
    k6_rejected_set = set().union(*k6_stage_sets)
    k6_rejected = [index for index in k6_base if index in k6_rejected_set]
    if (
        len(k6_rejected) != 155
        or k6_residue != [index for index in k6_base if index not in k6_rejected_set]
        or len(k6_residue) != 822
        or arbitrary_result.get("exact_result", {}).get("input_indices_sha256")
        != stable_hash(hereditary_residue)
    ):
        raise ValueError("K6 exact 977-to-822 partition failed")

    k7_graphs = ordered_records(k7_residue, k7_lookup)
    k6_graphs = ordered_records(k6_residue, k6_lookup)
    for record in k7_graphs:
        rows = record["adjacency"]
        if not alpha_at_most_two(rows) or not contains_clique(rows, 7):
            raise ValueError(f"K7 graph {record['index']} failed class validation")
    for record in k6_graphs:
        rows = record["adjacency"]
        if (
            not alpha_at_most_two(rows)
            or not contains_clique(rows, 6)
            or contains_clique(rows, 7)
        ):
            raise ValueError(f"K6-only graph {record['index']} failed class validation")
    if set(k7_residue) & set(k6_residue):
        raise ValueError("final graph classes overlap")

    positive = validate_adjacency(positive_18_graph(), 18, "positive18")
    if (
        not alpha_at_most_two(positive)
        or count_cliques(positive, 6) != 32
        or count_cliques(positive, 7) != 0
    ):
        raise ValueError("positive18 structural control failed")

    def stage_payload(
        name: str, inputs: Sequence[int], rejected: Sequence[int], residue: Sequence[int]
    ) -> dict:
        return {
            "name": name,
            "input_count": len(inputs),
            "input_indices_sha256": stable_hash(list(inputs)),
            "rejection_count": len(rejected),
            "rejected_indices": list(rejected),
            "rejected_indices_sha256": stable_hash(list(rejected)),
            "residue_count": len(residue),
            "residue_indices_sha256": stable_hash(list(residue)),
        }

    combined_indices = k7_residue + k6_residue
    tagged_indices = (
        [{"class": "K7", "index": index} for index in k7_residue]
        + [{"class": "K6_only", "index": index} for index in k6_residue]
    )
    tagged_graphs = (
        [{"class": "K7", **record} for record in k7_graphs]
        + [{"class": "K6_only", **record} for record in k6_graphs]
    )
    source_hashes = dict(observed)
    source_hashes[Path(__file__).name] = sha256(Path(__file__).resolve())
    return {
        "schema": "d6-current-exact-residue-v3",
        "status": "COMPLETE",
        "description": (
            "Self-contained exact unresolved level-19 boundary after the K7 "
            "support-or-pentad conjunction and the K6-only PSD-Z, hereditary, "
            "and arbitrary-subset Hall reductions."
        ),
        "upstream_result_commits": {
            "K7": K7_RESULT_COMMIT,
            "K6_only": K6_RESULT_COMMIT,
        },
        "upstream_source_commits": {
            "K7": K7_SOURCE_COMMIT,
            "K6_only": K6_SOURCE_COMMIT,
        },
        "class_order": ["K7", "K6_only"],
        "classes": {
            "K7": {
                "v2_boundary_count": len(k7_base),
                "v2_boundary_indices_sha256": stable_hash(k7_base),
                "stage_rejections": {
                    "joint_support": stage_payload(
                        "joint_support", k7_base, joint, joint_residue,
                    ),
                    "pentad_conjunction_increment": stage_payload(
                        "pentad_conjunction_increment", joint_residue,
                        incremental, k7_residue,
                    ),
                    "combined": stage_payload(
                        "combined_support_or_pentad", k7_base,
                        k7_rejected, k7_residue,
                    ),
                },
                "exact_rejections_since_v2": len(k7_rejected),
                "exact_rejected_indices": k7_rejected,
                "exact_rejected_indices_sha256": stable_hash(k7_rejected),
                "residue_count": len(k7_residue),
                "residue_indices": k7_residue,
                "residue_indices_sha256": stable_hash(k7_residue),
                "graphs": k7_graphs,
                "graphs_sha256": stable_hash(k7_graphs),
                "graph_validation": {
                    "symmetric_loopless_19_vertex": len(k7_graphs),
                    "alpha_at_most_two": len(k7_graphs),
                    "contains_K7": len(k7_graphs),
                },
                "adjacency_origin": {
                    "path": "d6_current_residue_manifest_v2.json",
                    "sha256": observed["d6_current_residue_manifest_v2.json"],
                    "embedded": True,
                },
            },
            "K6_only": {
                "v2_boundary_count": len(k6_base),
                "v2_boundary_indices_sha256": stable_hash(k6_base),
                "stage_rejections": [
                    stage_payload("psd_zmatrix", k6_base, psd_rejected, psd_residue),
                    stage_payload(
                        "psd_z_hereditary", psd_residue,
                        hereditary_rejected, hereditary_residue,
                    ),
                    stage_payload(
                        "arbitrary_subset_hall", hereditary_residue,
                        arbitrary_rejected, k6_residue,
                    ),
                ],
                "exact_rejections_since_v2": len(k6_rejected),
                "exact_rejected_indices": k6_rejected,
                "exact_rejected_indices_sha256": stable_hash(k6_rejected),
                "residue_count": len(k6_residue),
                "residue_indices": k6_residue,
                "residue_indices_sha256": stable_hash(k6_residue),
                "graphs": k6_graphs,
                "graphs_sha256": stable_hash(k6_graphs),
                "graph_validation": {
                    "symmetric_loopless_19_vertex": len(k6_graphs),
                    "alpha_at_most_two": len(k6_graphs),
                    "contains_K6": len(k6_graphs),
                    "contains_K7": 0,
                },
                "adjacency_origin": {
                    "path": "d6_current_residue_manifest_v2.json",
                    "sha256": observed["d6_current_residue_manifest_v2.json"],
                    "embedded": True,
                },
            },
        },
        "combined": {
            "count": len(combined_indices),
            "class_counts": {"K7": len(k7_residue), "K6_only": len(k6_residue)},
            "ordered_indices_sha256": stable_hash(combined_indices),
            "sorted_indices_sha256": stable_hash(sorted(combined_indices)),
            "ordered_class_index_records_sha256": stable_hash(tagged_indices),
            "ordered_class_graph_records_sha256": stable_hash(tagged_graphs),
            "cross_class_overlap": 0,
            "cross_class_overlap_indices_sha256": stable_hash([]),
        },
        "positive_18_control": {
            "name": "known realizable 18-point construction",
            "in_level_19_corpus": False,
            "vertices": 18,
            "adjacency": positive,
            "adjacency_sha256": stable_hash(positive),
            "symmetric_loopless": True,
            "alpha_at_most_two": True,
            "K6_seeds": 32,
            "K7_seeds": 0,
            "upstream_exact_checks_passed": True,
        },
        "verification_gates": {
            "v2": {
                "manifest_sha256": observed["d6_current_residue_manifest_v2.json"],
                "verification_sha256": observed[
                    "d6_current_residue_manifest_v2_verification.json"
                ],
                "status": "PASS",
            },
            "K7": {
                "result_commit": K7_RESULT_COMMIT,
                "source_commit": K7_SOURCE_COMMIT,
                "report_sha256": observed[
                    "d6_k7_pentad_joint_conjunction_report.json"
                ],
                "verification_sha256": observed[
                    "d6_k7_pentad_joint_conjunction_verification.json"
                ],
                "result_manifest_sha256": observed[
                    "d6_k7_pentad_joint_conjunction_result_manifest.json"
                ],
                "status": "PASS",
            },
            "K6_only": {
                "result_commit": K6_RESULT_COMMIT,
                "source_commit": K6_SOURCE_COMMIT,
                "psd_zmatrix_report_sha256": observed[
                    "d6_k6_psd_zmatrix_report.json"
                ],
                "psd_zmatrix_verification_sha256": observed[
                    "d6_k6_psd_zmatrix_verification.json"
                ],
                "hereditary_report_sha256": observed[
                    "d6_k6_psd_z_hereditary_report.json"
                ],
                "hereditary_verification_sha256": observed[
                    "d6_k6_psd_z_hereditary_verification.json"
                ],
                "arbitrary_report_sha256": observed[
                    "d6_k6_arbitrary_subset_hall_report.json"
                ],
                "arbitrary_verification_sha256": observed[
                    "d6_k6_arbitrary_subset_hall_verification.json"
                ],
                "result_manifest_sha256": observed[
                    "d6_k6_arbitrary_subset_hall_result_manifest.json"
                ],
                "status": "PASS",
            },
        },
        "source_hashes": source_hashes,
        "semantics": {
            "unit_graph_edges": "required to have Euclidean distance exactly 1",
            "candidate_nonedges": "unconstrained and may also have distance 1",
            "points": "must be distinct",
            "allowed_defect_coordinates": "upper bounds only and may be zero",
            "K7_rule": "full seed-local joint-support-or-pentad disjunction",
            "K6_arbitrary_subset_rule": "one selected subset per original span",
            "K6_one_vertex_extension": "not used",
            "membership_arithmetic": "exact set filtering and integer bitmasks",
            "unresolved": "never counted as rejected",
            "residue": "filter non-rejection only; not a realizability claim",
        },
        "nonclaims": [
            "No graph in the residue is claimed realizable in R6.",
            "No candidate nonedge is constrained to be non-unit.",
            "The combined count 977 is not yet a proof that f(6)=18.",
            "The positive 18-point control is separate from the level-19 corpus.",
            "This manifest does not replace upstream exact certificates.",
        ],
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
        "--output", type=Path,
        default=ROOT / "d6_current_residue_manifest_v3.json",
    )
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output, manifest)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "K7": manifest["classes"]["K7"]["residue_count"],
        "K6_only": manifest["classes"]["K6_only"]["residue_count"],
        "combined": manifest["combined"]["count"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
