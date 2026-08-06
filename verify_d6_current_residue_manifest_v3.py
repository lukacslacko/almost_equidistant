#!/usr/bin/env python3
"""Structurally independent checker for the exact v3 residue manifest.

The checker imports neither the v3 builder nor a K7/K6 production engine.  It
reconstructs the ordered 155-graph K7 and 822-graph K6-only subsets directly
from the committed v2 adjacencies and upstream result reports, checks the full
K6 977 -> 861 -> 831 -> 822 chain, validates every embedded graph, and checks
the separate positive18 control.
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
BUILDER = "build_d6_current_residue_manifest_v3.py"
EXPECTED_BUILDER_SHA256 = (
    "c0c94c461f5c62576ac79a3390777d6f43bf50b650f4e421a6ed951c2aeaf0d9"
)
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

POSITIVE_18 = [
    202622, 207805, 216539, 233191, 225767, 242395, 251069, 256126,
    228887, 245035, 252749, 255886, 255857, 252850, 245204, 229096,
    65535, 65535,
]


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
        raise ValueError(f"{name} is not an object")
    return value


def as_indices(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    result = list(value)
    if len(result) != len(set(result)):
        raise ValueError(f"{label} has duplicates")
    return result


def require_equal(observed: object, expected: object, label: str) -> None:
    if observed != expected:
        raise ValueError(f"{label} mismatch")


def validate_rows(rows: object, vertices: int, label: str) -> tuple[int, ...]:
    if (
        not isinstance(rows, list) or len(rows) != vertices
        or any(type(row) is not int for row in rows)
    ):
        raise ValueError(f"{label} has a bad adjacency shape")
    bound = 1 << vertices
    result = tuple(rows)
    for vertex, row in enumerate(result):
        if row < 0 or row >= bound or row & (1 << vertex):
            raise ValueError(f"{label} has an invalid row")
    for first, second in itertools.combinations(range(vertices), 2):
        if bool(result[first] & (1 << second)) != bool(
            result[second] & (1 << first)
        ):
            raise ValueError(f"{label} is asymmetric")
    return result


def independent_alpha_two(rows: Sequence[int]) -> bool:
    return not any(
        not (rows[a] & (1 << b))
        and not (rows[a] & (1 << c))
        and not (rows[b] & (1 << c))
        for a, b, c in itertools.combinations(range(len(rows)), 3)
    )


def independent_has_clique(rows: Sequence[int], size: int) -> bool:
    def extend(candidates: tuple[int, ...], need: int) -> bool:
        if need == 0:
            return True
        if len(candidates) < need:
            return False
        for position, vertex in enumerate(candidates):
            tail = tuple(
                other for other in candidates[position + 1:]
                if rows[vertex] & (1 << other)
            )
            if extend(tail, need - 1):
                return True
        return False

    return extend(tuple(range(len(rows))), size)


def count_cliques(rows: Sequence[int], size: int) -> int:
    count = 0
    for group in itertools.combinations(range(len(rows)), size):
        if all(rows[a] & (1 << b) for a, b in itertools.combinations(group, 2)):
            count += 1
    return count


def source_graphs(class_payload: dict, label: str) -> tuple[list[int], dict[int, list[int]]]:
    records = class_payload.get("graphs")
    if not isinstance(records, list):
        raise ValueError(f"{label} graph records missing")
    order = []
    lookup = {}
    for record in records:
        if not isinstance(record, dict) or type(record.get("index")) is not int:
            raise ValueError(f"{label} has a malformed graph record")
        index = int(record["index"])
        if index in lookup:
            raise ValueError(f"{label} repeats an index")
        order.append(index)
        lookup[index] = list(validate_rows(record.get("adjacency"), 19, label))
    return order, lookup


def independently_check_stage(
    report: dict,
    verification: dict,
    input_order: Sequence[int],
    report_schema: str,
    verifier_schema: str,
    rejection_count: int,
) -> tuple[list[int], list[int]]:
    rows = report.get("graph_results")
    if not isinstance(rows, list):
        raise ValueError("stage graph results missing")
    order = [int(row["index"]) for row in rows]
    rejected_by_rows = [
        int(row["index"]) for row in rows
        if row.get("decision", {}).get("rejected") is True
    ]
    rejected = as_indices(report.get("rejected_indices"), "stage rejected")
    independently_rejected = as_indices(
        verification.get("rejected_indices"), "verified stage rejected"
    )
    survivors = [index for index in input_order if index not in set(rejected)]
    if (
        report.get("schema") != report_schema or report.get("status") != "COMPLETE"
        or verification.get("schema") != verifier_schema
        or verification.get("status") != "PASS"
        or order != list(input_order)
        or report.get("input_graphs") != len(input_order)
        or report.get("input_indices_sha256") != stable_hash(list(input_order))
        or rejected != rejected_by_rows or rejected != independently_rejected
        or len(rejected) != rejection_count
        or report.get("graphs_rejected") != rejection_count
        or report.get("graphs_surviving") != len(survivors)
        or verification.get("graphs_recomputed") != len(input_order)
        or verification.get("graphs_rejected") != rejection_count
        or verification.get("graphs_surviving") != len(survivors)
        or report.get("positive_control", {}).get("passed") is not True
        or verification.get("positive_control", {}).get("passed") is not True
    ):
        raise ValueError(f"independent stage check failed: {report_schema}")
    return rejected, survivors


def stage_record(
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


def verify(manifest_path: Path) -> dict:
    observed_files = {name: sha256(ROOT / name) for name in EXPECTED_FILES}
    if observed_files != EXPECTED_FILES:
        raise ValueError("upstream root hash mismatch")
    builder_hash = sha256(ROOT / BUILDER)
    if builder_hash != EXPECTED_BUILDER_SHA256:
        raise ValueError("v3 builder source differs from checker pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("v3 manifest is not an object")

    v2 = load("d6_current_residue_manifest_v2.json")
    v2_check = load("d6_current_residue_manifest_v2_verification.json")
    if (
        v2.get("schema") != "d6-current-exact-residue-v2"
        or v2.get("status") != "COMPLETE"
        or v2_check.get("status") != "PASS"
        or v2_check.get("manifest", {}).get("sha256")
        != EXPECTED_FILES["d6_current_residue_manifest_v2.json"]
    ):
        raise ValueError("v2 independent gate failed")
    k7_source = v2["classes"]["K7"]
    k6_source = v2["classes"]["K6_only"]
    k7_base = as_indices(k7_source.get("residue_indices"), "v2 K7")
    k6_base = as_indices(k6_source.get("residue_indices"), "v2 K6")
    k7_order, k7_lookup = source_graphs(k7_source, "v2 K7")
    k6_order, k6_lookup = source_graphs(k6_source, "v2 K6")
    if (
        k7_base != k7_order or len(k7_base) != 258
        or stable_hash(k7_base) != k7_source.get("residue_indices_sha256")
        or stable_hash(k7_source["graphs"]) != k7_source.get("graphs_sha256")
        or k6_base != k6_order or len(k6_base) != 977
        or stable_hash(k6_base) != k6_source.get("residue_indices_sha256")
        or stable_hash(k6_source["graphs"]) != k6_source.get("graphs_sha256")
    ):
        raise ValueError("v2 embedded adjacency/order gate failed")

    k7_report = load("d6_k7_pentad_joint_conjunction_report.json")
    k7_verification = load("d6_k7_pentad_joint_conjunction_verification.json")
    k7_result = load("d6_k7_pentad_joint_conjunction_result_manifest.json")
    joint = as_indices(
        k7_report["joint_support_rejections"].get("ordered_indices"), "K7 joint"
    )
    incremental = as_indices(
        k7_report["new_seed_conjunction_rejections"].get("ordered_indices"),
        "K7 increment",
    )
    joint_residue = [index for index in k7_base if index not in set(joint)]
    rejection_set = set(joint) | set(incremental)
    k7_rejected = [index for index in k7_base if index in rejection_set]
    k7_residue = [index for index in k7_base if index not in rejection_set]
    if (
        len(joint) != 69 or len(incremental) != 34 or set(joint) & set(incremental)
        or not set(incremental) <= set(joint_residue)
        or len(k7_rejected) != 103 or len(k7_residue) != 155
        or k7_report.get("status") != "COMPLETE"
        or k7_report["combined_rejections"].get("ordered_indices") != k7_rejected
        or k7_report["exact_residue"].get("ordered_indices") != k7_residue
        or k7_verification.get("status") != "PASS"
        or k7_verification.get("report", {}).get("sha256")
        != EXPECTED_FILES["d6_k7_pentad_joint_conjunction_report.json"]
        or k7_verification.get("ordered_rejections_sha256") != stable_hash(k7_rejected)
        or k7_verification.get("ordered_residue_sha256") != stable_hash(k7_residue)
        or k7_result.get("status") != "COMPLETE"
        or k7_result.get("source_boundary", {}).get("commit") != K7_SOURCE_COMMIT
        or k7_result.get("source_boundary", {}).get("builder_sha256")
        != EXPECTED_FILES["build_d6_k7_pentad_joint_conjunction.py"]
        or k7_result.get("source_boundary", {}).get("verifier_sha256")
        != EXPECTED_FILES["verify_d6_k7_pentad_joint_conjunction.py"]
        or k7_result.get("artifacts", {}).get("report", {}).get("path")
        != "d6_k7_pentad_joint_conjunction_report.json"
        or k7_result.get("artifacts", {}).get("report", {}).get("sha256")
        != EXPECTED_FILES["d6_k7_pentad_joint_conjunction_report.json"]
        or k7_result.get("artifacts", {}).get("report", {}).get("bytes")
        != (ROOT / "d6_k7_pentad_joint_conjunction_report.json").stat().st_size
        or k7_result.get("artifacts", {}).get("verification", {}).get("path")
        != "d6_k7_pentad_joint_conjunction_verification.json"
        or k7_result.get("artifacts", {}).get("verification", {}).get("sha256")
        != EXPECTED_FILES["d6_k7_pentad_joint_conjunction_verification.json"]
        or k7_result.get("artifacts", {}).get("verification", {}).get("bytes")
        != (ROOT / "d6_k7_pentad_joint_conjunction_verification.json").stat().st_size
        or k7_result.get("exact_result", {}).get("input_graphs") != 258
        or k7_result.get("exact_result", {}).get("input_indices_sha256")
        != stable_hash(k7_base)
        or k7_result.get("exact_result", {}).get("prior_joint_support_rejections")
        != 69
        or k7_result.get("exact_result", {}).get("new_pentad_conjunction_rejections")
        != 34
        or k7_result.get("exact_result", {}).get("combined_rejected_graphs")
        != 103
        or k7_result.get("exact_result", {}).get(
            "combined_rejected_indices_sha256"
        ) != stable_hash(k7_rejected)
        or k7_result.get("exact_result", {}).get("residue_graphs") != 155
        or k7_result.get("exact_result", {}).get("residue_indices_sha256")
        != stable_hash(k7_residue)
    ):
        raise ValueError("independent K7 258-to-155 reconstruction failed")

    psd = load("d6_k6_psd_zmatrix_report.json")
    psd_v = load("d6_k6_psd_zmatrix_verification.json")
    hereditary = load("d6_k6_psd_z_hereditary_report.json")
    hereditary_v = load("d6_k6_psd_z_hereditary_verification.json")
    arbitrary = load("d6_k6_arbitrary_subset_hall_report.json")
    arbitrary_v = load("d6_k6_arbitrary_subset_hall_verification.json")
    arbitrary_result = load("d6_k6_arbitrary_subset_hall_result_manifest.json")
    psd_rejected, psd_residue = independently_check_stage(
        psd, psd_v, k6_base, "d6-k6-psd-zmatrix-v1",
        "d6-k6-psd-zmatrix-verification-v1", 116,
    )
    hereditary_rejected, hereditary_residue = independently_check_stage(
        hereditary, hereditary_v, psd_residue, "d6-k6-psd-z-hereditary-v1",
        "d6-k6-psd-z-hereditary-verification-v1", 30,
    )
    arbitrary_rejected, k6_residue = independently_check_stage(
        arbitrary, arbitrary_v, hereditary_residue,
        "d6-k6-arbitrary-subset-hall-v1",
        "d6-k6-arbitrary-subset-hall-verification-v1", 9,
    )
    for report, certificate_name, checkpoint_name in (
        (psd, "d6_k6_psd_zmatrix_certificates.json", None),
        (
            hereditary,
            "d6_k6_psd_z_hereditary_certificates.json",
            "d6_k6_psd_z_hereditary_checkpoint.json",
        ),
        (
            arbitrary,
            "d6_k6_arbitrary_subset_hall_certificates.json",
            "d6_k6_arbitrary_subset_hall_checkpoint.json",
        ),
    ):
        archive = report.get("certificate_archive", {})
        if (
            archive.get("path") != certificate_name
            or archive.get("sha256") != EXPECTED_FILES[certificate_name]
            or archive.get("bytes") != (ROOT / certificate_name).stat().st_size
        ):
            raise ValueError(f"K6 certificate binding failed: {certificate_name}")
        if checkpoint_name is not None:
            checkpoint = report.get("checkpoint", {})
            if (
                checkpoint.get("path") != checkpoint_name
                or checkpoint.get("sha256") != EXPECTED_FILES[checkpoint_name]
            ):
                raise ValueError(f"K6 checkpoint binding failed: {checkpoint_name}")
    stage_sets = [set(psd_rejected), set(hereditary_rejected), set(arbitrary_rejected)]
    if any(stage_sets[a] & stage_sets[b] for a, b in itertools.combinations(range(3), 2)):
        raise ValueError("K6 independent stage sets overlap")
    k6_rejection_set = set().union(*stage_sets)
    k6_rejected = [index for index in k6_base if index in k6_rejection_set]
    if (
        len(k6_rejected) != 155 or len(k6_residue) != 822
        or k6_residue != [index for index in k6_base if index not in k6_rejection_set]
        or arbitrary_result.get("status") != "COMPLETE"
        or arbitrary_result.get("source_boundary", {}).get("commit") != K6_SOURCE_COMMIT
        or arbitrary_result.get("source_boundary", {}).get("production_sha256")
        != EXPECTED_FILES["d6_k6_arbitrary_subset_hall.py"]
        or arbitrary_result.get("source_boundary", {}).get("verifier_sha256")
        != EXPECTED_FILES["verify_d6_k6_arbitrary_subset_hall.py"]
        or arbitrary_result.get("parent", {}).get("sha256")
        != EXPECTED_FILES["d6_k6_psd_z_hereditary_report.json"]
        or arbitrary_result.get("artifacts", {}).get("report", {}).get("path")
        != "d6_k6_arbitrary_subset_hall_report.json"
        or arbitrary_result.get("artifacts", {}).get("report", {}).get("sha256")
        != EXPECTED_FILES["d6_k6_arbitrary_subset_hall_report.json"]
        or arbitrary_result.get("artifacts", {}).get("verification", {}).get("path")
        != "d6_k6_arbitrary_subset_hall_verification.json"
        or arbitrary_result.get("artifacts", {}).get("verification", {}).get("sha256")
        != EXPECTED_FILES["d6_k6_arbitrary_subset_hall_verification.json"]
        or arbitrary_result.get("artifacts", {}).get("certificates", {}).get("path")
        != "d6_k6_arbitrary_subset_hall_certificates.json"
        or arbitrary_result.get("artifacts", {}).get("certificates", {}).get("sha256")
        != EXPECTED_FILES["d6_k6_arbitrary_subset_hall_certificates.json"]
        or arbitrary_result.get("artifacts", {}).get("checkpoint", {}).get("path")
        != "d6_k6_arbitrary_subset_hall_checkpoint.json"
        or arbitrary_result.get("artifacts", {}).get("checkpoint", {}).get("sha256")
        != EXPECTED_FILES["d6_k6_arbitrary_subset_hall_checkpoint.json"]
        or any(
            arbitrary_result.get("artifacts", {}).get(key, {}).get("bytes")
            != (ROOT / path).stat().st_size
            for key, path in (
                ("report", "d6_k6_arbitrary_subset_hall_report.json"),
                ("verification", "d6_k6_arbitrary_subset_hall_verification.json"),
                ("certificates", "d6_k6_arbitrary_subset_hall_certificates.json"),
                ("checkpoint", "d6_k6_arbitrary_subset_hall_checkpoint.json"),
            )
        )
        or arbitrary_result.get("exact_result", {}).get("input_graphs") != 831
        or arbitrary_result.get("exact_result", {}).get("input_indices_sha256")
        != stable_hash(hereditary_residue)
        or arbitrary_result.get("exact_result", {}).get("rejected_graphs") != 9
        or arbitrary_result.get("exact_result", {}).get("rejected_indices")
        != arbitrary_rejected
        or arbitrary_result.get("exact_result", {}).get("rejected_indices_sha256")
        != stable_hash(arbitrary_rejected)
        or arbitrary_result.get("exact_result", {}).get("surviving_graphs") != 822
    ):
        raise ValueError("independent K6 977-to-822 reconstruction failed")

    k7_graphs = [{"index": index, "adjacency": k7_lookup[index]} for index in k7_residue]
    k6_graphs = [{"index": index, "adjacency": k6_lookup[index]} for index in k6_residue]
    for record in k7_graphs:
        rows = validate_rows(record["adjacency"], 19, "final K7")
        if not independent_alpha_two(rows) or not independent_has_clique(rows, 7):
            raise ValueError(f"final K7 graph {record['index']} is misclassified")
    for record in k6_graphs:
        rows = validate_rows(record["adjacency"], 19, "final K6-only")
        if (
            not independent_alpha_two(rows)
            or not independent_has_clique(rows, 6)
            or independent_has_clique(rows, 7)
        ):
            raise ValueError(f"final K6-only graph {record['index']} is misclassified")
    if set(k7_residue) & set(k6_residue):
        raise ValueError("independent final class disjointness failed")

    expected_k7 = {
        "v2_boundary_count": 258,
        "v2_boundary_indices_sha256": stable_hash(k7_base),
        "stage_rejections": {
            "joint_support": stage_record(
                "joint_support", k7_base, joint, joint_residue,
            ),
            "pentad_conjunction_increment": stage_record(
                "pentad_conjunction_increment",
                joint_residue,
                incremental,
                k7_residue,
            ),
            "combined": stage_record(
                "combined_support_or_pentad", k7_base, k7_rejected, k7_residue
            ),
        },
        "exact_rejections_since_v2": 103,
        "exact_rejected_indices": k7_rejected,
        "exact_rejected_indices_sha256": stable_hash(k7_rejected),
        "residue_count": 155,
        "residue_indices": k7_residue,
        "residue_indices_sha256": stable_hash(k7_residue),
        "graphs": k7_graphs,
        "graphs_sha256": stable_hash(k7_graphs),
        "graph_validation": {
            "symmetric_loopless_19_vertex": 155,
            "alpha_at_most_two": 155,
            "contains_K7": 155,
        },
        "adjacency_origin": {
            "path": "d6_current_residue_manifest_v2.json",
            "sha256": EXPECTED_FILES["d6_current_residue_manifest_v2.json"],
            "embedded": True,
        },
    }
    expected_k6 = {
        "v2_boundary_count": 977,
        "v2_boundary_indices_sha256": stable_hash(k6_base),
        "stage_rejections": [
            stage_record("psd_zmatrix", k6_base, psd_rejected, psd_residue),
            stage_record(
                "psd_z_hereditary", psd_residue,
                hereditary_rejected, hereditary_residue,
            ),
            stage_record(
                "arbitrary_subset_hall", hereditary_residue,
                arbitrary_rejected, k6_residue,
            ),
        ],
        "exact_rejections_since_v2": 155,
        "exact_rejected_indices": k6_rejected,
        "exact_rejected_indices_sha256": stable_hash(k6_rejected),
        "residue_count": 822,
        "residue_indices": k6_residue,
        "residue_indices_sha256": stable_hash(k6_residue),
        "graphs": k6_graphs,
        "graphs_sha256": stable_hash(k6_graphs),
        "graph_validation": {
            "symmetric_loopless_19_vertex": 822,
            "alpha_at_most_two": 822,
            "contains_K6": 822,
            "contains_K7": 0,
        },
        "adjacency_origin": {
            "path": "d6_current_residue_manifest_v2.json",
            "sha256": EXPECTED_FILES["d6_current_residue_manifest_v2.json"],
            "embedded": True,
        },
    }
    require_equal(manifest.get("classes", {}).get("K7"), expected_k7, "K7 class")
    require_equal(manifest.get("classes", {}).get("K6_only"), expected_k6, "K6 class")
    require_equal(set(manifest.get("classes", {})), {"K7", "K6_only"}, "class keys")

    concatenated = k7_residue + k6_residue
    tagged_indices = (
        [{"class": "K7", "index": index} for index in k7_residue]
        + [{"class": "K6_only", "index": index} for index in k6_residue]
    )
    tagged_graphs = (
        [{"class": "K7", **record} for record in k7_graphs]
        + [{"class": "K6_only", **record} for record in k6_graphs]
    )
    expected_combined = {
        "count": 977,
        "class_counts": {"K7": 155, "K6_only": 822},
        "ordered_indices_sha256": stable_hash(concatenated),
        "sorted_indices_sha256": stable_hash(sorted(concatenated)),
        "ordered_class_index_records_sha256": stable_hash(tagged_indices),
        "ordered_class_graph_records_sha256": stable_hash(tagged_graphs),
        "cross_class_overlap": 0,
        "cross_class_overlap_indices_sha256": stable_hash([]),
    }
    require_equal(manifest.get("combined"), expected_combined, "combined accounting")

    positive = validate_rows(list(POSITIVE_18), 18, "positive18")
    if (
        not independent_alpha_two(positive)
        or count_cliques(positive, 6) != 32
        or count_cliques(positive, 7) != 0
    ):
        raise ValueError("independent positive18 reconstruction failed")
    expected_positive = {
        "name": "known realizable 18-point construction",
        "in_level_19_corpus": False,
        "vertices": 18,
        "adjacency": list(POSITIVE_18),
        "adjacency_sha256": stable_hash(list(POSITIVE_18)),
        "symmetric_loopless": True,
        "alpha_at_most_two": True,
        "K6_seeds": 32,
        "K7_seeds": 0,
        "upstream_exact_checks_passed": True,
    }
    require_equal(manifest.get("positive_18_control"), expected_positive, "positive18")

    expected_gates = {
        "v2": {
            "manifest_sha256": EXPECTED_FILES["d6_current_residue_manifest_v2.json"],
            "verification_sha256": EXPECTED_FILES[
                "d6_current_residue_manifest_v2_verification.json"
            ],
            "status": "PASS",
        },
        "K7": {
            "result_commit": K7_RESULT_COMMIT,
            "source_commit": K7_SOURCE_COMMIT,
            "report_sha256": EXPECTED_FILES[
                "d6_k7_pentad_joint_conjunction_report.json"
            ],
            "verification_sha256": EXPECTED_FILES[
                "d6_k7_pentad_joint_conjunction_verification.json"
            ],
            "result_manifest_sha256": EXPECTED_FILES[
                "d6_k7_pentad_joint_conjunction_result_manifest.json"
            ],
            "status": "PASS",
        },
        "K6_only": {
            "result_commit": K6_RESULT_COMMIT,
            "source_commit": K6_SOURCE_COMMIT,
            "psd_zmatrix_report_sha256": EXPECTED_FILES[
                "d6_k6_psd_zmatrix_report.json"
            ],
            "psd_zmatrix_verification_sha256": EXPECTED_FILES[
                "d6_k6_psd_zmatrix_verification.json"
            ],
            "hereditary_report_sha256": EXPECTED_FILES[
                "d6_k6_psd_z_hereditary_report.json"
            ],
            "hereditary_verification_sha256": EXPECTED_FILES[
                "d6_k6_psd_z_hereditary_verification.json"
            ],
            "arbitrary_report_sha256": EXPECTED_FILES[
                "d6_k6_arbitrary_subset_hall_report.json"
            ],
            "arbitrary_verification_sha256": EXPECTED_FILES[
                "d6_k6_arbitrary_subset_hall_verification.json"
            ],
            "result_manifest_sha256": EXPECTED_FILES[
                "d6_k6_arbitrary_subset_hall_result_manifest.json"
            ],
            "status": "PASS",
        },
    }
    require_equal(manifest.get("verification_gates"), expected_gates, "gates")
    expected_sources = dict(EXPECTED_FILES)
    expected_sources[BUILDER] = builder_hash
    require_equal(manifest.get("source_hashes"), expected_sources, "source hashes")

    expected_semantics = {
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
    }
    expected_nonclaims = [
        "No graph in the residue is claimed realizable in R6.",
        "No candidate nonedge is constrained to be non-unit.",
        "The combined count 977 is not yet a proof that f(6)=18.",
        "The positive 18-point control is separate from the level-19 corpus.",
        "This manifest does not replace upstream exact certificates.",
    ]
    require_equal(manifest.get("semantics"), expected_semantics, "semantics")
    require_equal(manifest.get("nonclaims"), expected_nonclaims, "nonclaims")
    require_equal(manifest.get("class_order"), ["K7", "K6_only"], "class order")
    require_equal(
        manifest.get("upstream_result_commits"),
        {"K7": K7_RESULT_COMMIT, "K6_only": K6_RESULT_COMMIT},
        "result commits",
    )
    require_equal(
        manifest.get("upstream_source_commits"),
        {"K7": K7_SOURCE_COMMIT, "K6_only": K6_SOURCE_COMMIT},
        "source commits",
    )
    require_equal(
        set(manifest),
        {
            "schema", "status", "description", "upstream_result_commits",
            "upstream_source_commits", "class_order", "classes", "combined",
            "positive_18_control", "verification_gates", "source_hashes",
            "semantics", "nonclaims",
        },
        "top-level keys",
    )
    require_equal(manifest.get("schema"), "d6-current-exact-residue-v3", "schema")
    require_equal(manifest.get("status"), "COMPLETE", "status")

    return {
        "schema": "d6-current-exact-residue-v3-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": sha256(manifest_path)},
        "upstream_result_commits": {
            "K7": K7_RESULT_COMMIT,
            "K6_only": K6_RESULT_COMMIT,
        },
        "classes": {
            "K7": {
                "v2_input": 258,
                "exact_rejections": 103,
                "residue": 155,
                "residue_indices_sha256": stable_hash(k7_residue),
                "adjacency_records_checked": len(k7_graphs),
                "adjacency_records_sha256": stable_hash(k7_graphs),
            },
            "K6_only": {
                "v2_input": 977,
                "stage_rejections": [116, 30, 9],
                "exact_rejections": 155,
                "residue": 822,
                "residue_indices_sha256": stable_hash(k6_residue),
                "adjacency_records_checked": len(k6_graphs),
                "adjacency_records_sha256": stable_hash(k6_graphs),
            },
        },
        "combined": expected_combined,
        "positive_18_control": {
            "vertices": 18,
            "adjacency_sha256": stable_hash(list(POSITIVE_18)),
            "K6_seeds": 32,
            "K7_seeds": 0,
            "separate_from_level_19_corpus": True,
        },
        "checks": {
            "all_upstream_hashes": True,
            "K7_258_to_155_partition": True,
            "K6_977_to_861_to_831_to_822_partition": True,
            "class_disjointness": True,
            "all_977_adjacencies_match_v2": True,
            "all_977_symmetric_loopless": True,
            "all_977_alpha_at_most_two": True,
            "all_155_K7_graphs_contain_K7": True,
            "all_822_K6_only_graphs_contain_K6_not_K7": True,
            "positive18_separate_control": True,
        },
        "inputs": expected_sources,
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
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
        "--manifest", type=Path,
        default=ROOT / "d6_current_residue_manifest_v3.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_current_residue_manifest_v3_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.manifest)
    atomic_json(args.output, result)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "status": result["status"],
        "K7": result["classes"]["K7"]["residue"],
        "K6_only": result["classes"]["K6_only"]["residue"],
        "combined": result["combined"]["count"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
