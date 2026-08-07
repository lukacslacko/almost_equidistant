#!/usr/bin/env python3
"""Independent structural checker for the mixed-certificate d=6 v6 residue."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import time
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
BUILDER = "build_d6_current_residue_manifest_v6.py"
EXPECTED_BUILDER_SHA256 = (
    "32345fc272efec6e04cab76a158411f65b4ea7cc6970ad4d1442dc8ca8471918"
)
PRIMARY = {
    "d6_current_residue_manifest_v5.json": (
        "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5"
    ),
    "d6_current_residue_manifest_v5_verification.json": (
        "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46"
    ),
    "d6_k7_one_free_conjunction_report.json": (
        "181f63015d373fb69f37f4f390cdd8ad50bc32b5a2d97a36fd0bac876d0877a2"
    ),
    "d6_k7_one_free_conjunction_verification.json": (
        "ebaacfe663728eff427c9ba2280fc4261d26dfe10fad21e1a00982b9d31ac062"
    ),
    "d6_interval_v5_k7_cap100000_report.json": (
        "85b175f3b3348204c3ba9065b30432458528bdd2e22b471e64fc56341480bac3"
    ),
    "d6_interval_v5_k7_cap100000_verification.json": (
        "a56fb1561b778c55bffb433bf0eea90e5cca823573d0717e74f2a73528bc1b10"
    ),
    "d6_k6_repeated_two_support_arm_report.json": (
        "ea8d9438d062b92857a1057b950e76b8ae08bbe1bf97b4327130fd4058171ebe"
    ),
    "d6_k6_repeated_two_support_arm_verification.json": (
        "3abf48bab092671a9dff03c2dc47c0063a75515cd219e1d5a4bfa9a5fa797b5b"
    ),
}
EXPECTED_K7 = [
    316173,
    2581209,
    2592657,
    2593240,
    3595554,
    3648882,
    3729907,
    3785980,
    3888410,
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
]
EXPECTED_K7_SHA256 = (
    "af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290"
)
EXPECTED_K6_SHA256 = (
    "04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1"
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load(name_or_path: str | Path) -> dict:
    path = Path(name_or_path)
    if not path.is_absolute():
        path = ROOT / path
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path.name}: JSON object")
    return value


def check_named_hashes(mapping: object, label: str) -> dict[str, str]:
    require(isinstance(mapping, dict) and bool(mapping), f"{label}: hash map")
    answer = {}
    for name, expected in mapping.items():
        require(
            isinstance(name, str)
            and Path(name).name == name
            and isinstance(expected, str)
            and len(expected) == 64,
            f"{label}: malformed entry",
        )
        require(sha256(ROOT / name) == expected, f"{label}: {name}")
        answer[name] = expected
    return dict(sorted(answer.items()))


def validate_rows(value: object, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == 19
        and all(type(row) is int for row in value),
        f"{label}: adjacency",
    )
    rows = tuple(value)
    for first, row in enumerate(rows):
        require(0 <= row < (1 << 19), f"{label}: row range")
        require(not row & (1 << first), f"{label}: loop")
        for second in range(first):
            require(
                bool(row & (1 << second))
                == bool(rows[second] & (1 << first)),
                f"{label}: asymmetry",
            )
    return rows


def alpha_at_most_two(rows: Sequence[int]) -> bool:
    for first in range(19):
        for second in range(first):
            if rows[first] & (1 << second):
                continue
            for third in range(second):
                if (
                    not rows[first] & (1 << third)
                    and not rows[second] & (1 << third)
                ):
                    return False
    return True


def contains_clique(rows: Sequence[int], size: int) -> bool:
    def visit(candidates: int, need: int) -> bool:
        if not need:
            return True
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            if visit(candidates & rows[vertex], need - 1):
                return True
        return False

    return visit((1 << len(rows)) - 1, size)


def exact_partition(
    report: dict,
    verification: dict,
    *,
    schema: str,
    verification_schema: str,
    report_name: str,
    input_count: int,
    rejection_count: int,
    residue_count: int,
) -> tuple[list[int], list[int], list[int], dict[str, str]]:
    require(report.get("schema") == schema, f"{report_name}: schema")
    require(report.get("status") == "COMPLETE", f"{report_name}: status")
    require(report.get("input_graphs") == input_count, f"{report_name}: input")
    require(
        report.get("graphs_rejected") == rejection_count,
        f"{report_name}: rejected count",
    )
    require(
        report.get("graphs_surviving") == residue_count,
        f"{report_name}: residue count",
    )
    require(
        verification.get("schema") == verification_schema
        and verification.get("status") == "PASS",
        f"{report_name}: verification",
    )
    require(
        verification.get("report", {}).get("sha256")
        == sha256(ROOT / report_name),
        f"{report_name}: verification binding",
    )
    require(
        verification.get("graphs_recomputed") == input_count
        and verification.get("graphs_rejected") == rejection_count
        and verification.get("graphs_surviving") == residue_count,
        f"{report_name}: independently recomputed counts",
    )
    inputs = list(report.get("ordered_input_indices", ()))
    rejected = list(report.get("rejected_indices", ()))
    residue = list(report.get("ordered_residue_indices", ()))
    require(
        len(inputs) == input_count
        and len(set(inputs)) == input_count
        and len(rejected) == rejection_count
        and len(set(rejected)) == rejection_count
        and len(residue) == residue_count
        and len(set(residue)) == residue_count
        and all(type(index) is int for index in inputs + rejected + residue),
        f"{report_name}: index lists",
    )
    rejected_set = set(rejected)
    require(
        [index for index in inputs if index not in rejected_set] == residue,
        f"{report_name}: ordered subtraction",
    )
    require(
        stable_hash(inputs) == report.get("ordered_input_indices_sha256")
        and stable_hash(rejected) == report.get("rejected_indices_sha256")
        and stable_hash(residue) == report.get("ordered_residue_indices_sha256"),
        f"{report_name}: report hashes",
    )
    require(
        stable_hash(rejected) == verification.get("rejected_indices_sha256")
        and stable_hash(residue)
        == verification.get("ordered_residue_indices_sha256"),
        f"{report_name}: verification hashes",
    )
    dependencies = check_named_hashes(
        report.get("dependencies"), f"{report_name} dependencies"
    )
    archive = report.get("certificate_archive", {})
    require(
        isinstance(archive, dict)
        and isinstance(archive.get("path"), str)
        and Path(archive["path"]).name == archive["path"]
        and sha256(ROOT / archive["path"]) == archive.get("sha256")
        and archive.get("rejected_graphs") == rejection_count,
        f"{report_name}: archive",
    )
    source = report_name.replace("_report.json", ".py")
    require(
        sha256(ROOT / source) == report.get("production_source_sha256"),
        f"{report_name}: source",
    )
    semantics = report.get("semantics", {})
    require(
        "exact" in semantics.get("arithmetic", "")
        and "unconstrained" in semantics.get("candidate_nonedges", ""),
        f"{report_name}: semantics",
    )
    require(
        report.get("positive_18_control", {}).get("passed") is True
        and verification.get("positive_18_control", {}).get("passed") is True,
        f"{report_name}: positive control",
    )
    return inputs, rejected, residue, dependencies


def verify(manifest_path: Path, expected_manifest_sha256: str | None) -> dict:
    checks: dict[str, bool] = {}
    manifest_hash = sha256(manifest_path)
    if expected_manifest_sha256 is not None:
        require(manifest_hash == expected_manifest_sha256, "manifest hash")
    require(sha256(ROOT / BUILDER) == EXPECTED_BUILDER_SHA256, "builder hash")
    check_named_hashes(PRIMARY, "primary artifacts")
    checks["primary_hash_roots"] = True

    manifest = load(manifest_path)
    v5 = load("d6_current_residue_manifest_v5.json")
    v5_check = load("d6_current_residue_manifest_v5_verification.json")
    require(
        v5.get("schema") == "d6-current-exact-residue-v5"
        and v5.get("status") == "COMPLETE_EXACT_FILTER_UNION"
        and v5_check.get("status") == "PASS"
        and v5_check.get("manifest", {}).get("sha256")
        == PRIMARY["d6_current_residue_manifest_v5.json"],
        "v5 verified root",
    )
    v5_k7 = list(v5["classes"]["K7"]["graphs"])
    v5_k6 = list(v5["classes"]["K6_only"]["graphs"])
    v5_k7_indices = [record["index"] for record in v5_k7]
    v5_k6_indices = [record["index"] for record in v5_k6]
    require(len(v5_k7) == 155 and len(v5_k6) == 756, "v5 counts")

    exact = load("d6_k7_one_free_conjunction_report.json")
    exact_check = load("d6_k7_one_free_conjunction_verification.json")
    exact_upstream = check_named_hashes(
        exact.get("upstream_artifact_sha256"), "K7 exact upstream"
    )
    exact_sources = check_named_hashes(
        exact.get("source_sha256"), "K7 exact sources"
    )
    require(
        exact.get("kind") == "d6_k7_one_free_edge_seed_conjunction"
        and exact.get("status") == "COMPLETE"
        and exact_check.get("kind")
        == "d6_k7_one_free_edge_seed_conjunction_verification"
        and exact_check.get("status") == "PASS"
        and exact_check.get("report", {}).get("sha256")
        == PRIMARY["d6_k7_one_free_conjunction_report.json"]
        and exact_check.get("upstream_artifact_sha256") == exact_upstream,
        "K7 exact verified root",
    )
    base = load("d6_k7_double_pin_conjunction_report.json")
    base_check = load("d6_k7_double_pin_conjunction_verification.json")
    full = load("d6_k7_full_pin_increment_report.json")
    full_check = load("d6_k7_full_pin_increment_verification.json")
    require(
        base.get("input", {}).get("graphs") == 155
        and base.get("input", {}).get("indices_sha256")
        == stable_hash(v5_k7_indices)
        and base.get("input", {}).get("graphs_sha256") == stable_hash(v5_k7)
        and base_check.get("status") == "PASS"
        and base_check.get("report", {}).get("sha256")
        == exact_upstream["d6_k7_double_pin_conjunction_report.json"]
        and full.get("status") == "COMPLETE"
        and full_check.get("status") == "PASS"
        and full_check.get("report", {}).get("sha256")
        == exact_upstream["d6_k7_full_pin_increment_report.json"]
        and exact.get("input_indices") == full.get("exact_survivors"),
        "K7 exact transitive chain",
    )
    exact22 = list(exact["summary"]["survivors_by_pipeline"]["combined"])
    require(
        len(exact22) == 22
        and stable_hash(exact22)
        == "cce80d8065f5c97c15091343128f9626447b13e7bd6e7638b59c9dcef5296cf6"
        and exact_check.get("exact_survivors_sha256") == stable_hash(exact22),
        "K7 exact 22-list",
    )
    exact22_set = set(exact22)
    exact_rejected = [index for index in v5_k7_indices if index not in exact22_set]
    checks["K7_exact_chain"] = True

    interval = load("d6_interval_v5_k7_cap100000_report.json")
    interval_check = load("d6_interval_v5_k7_cap100000_verification.json")
    require(
        interval.get("schema") == 2
        and interval.get("completed_graphs") == 155
        and interval.get("total_graphs") == 155
        and interval.get("certified_killed") == 24
        and interval_check.get("schema")
        == "d6-interval-v5-verification-v1"
        and interval_check.get("status") == "PASS"
        and interval_check.get("report", {}).get("sha256")
        == PRIMARY["d6_interval_v5_k7_cap100000_report.json"]
        and interval_check.get("selection", {}).get("graphs") == 155
        and interval_check.get("selection", {}).get("indices_sha256")
        == stable_hash(v5_k7_indices)
        and interval_check.get("kernel_controls", {}).get("status") == "PASS"
        and interval_check.get("witness_replay", {}).get("status") == "PASS",
        "interval verified root",
    )
    killed = list(interval_check.get("killed_indices", ()))
    require(
        len(killed) == 24
        and len(set(killed)) == 24
        and set(killed) <= set(v5_k7_indices)
        and interval.get("status_counts", {}).get("KILLED") == 24
        and interval_check.get("status_counts", {}).get("KILLED") == 24,
        "interval KILLED list",
    )
    interval_new = [index for index in exact22 if index in set(killed)]
    interval_old = [index for index in killed if index not in exact22_set]
    k7_indices = [index for index in exact22 if index not in set(killed)]
    require(
        interval_new == [423661, 424226, 3936176]
        and len(interval_old) == 21
        and k7_indices == EXPECTED_K7
        and stable_hash(k7_indices) == EXPECTED_K7_SHA256,
        "K7 cross-method union",
    )
    trust = interval_check.get("trust_assumptions")
    require(
        isinstance(trust, dict)
        and set(trust)
        == {"basic_arithmetic", "candidate_nonedges", "transcendentals"}
        and "IEEE-754" in trust["basic_arithmetic"]
        and "8 ulps" in trust["transcendentals"],
        "interval trust assumptions",
    )
    checks["K7_interval_KILLED_union"] = True

    repeated = load("d6_k6_repeated_two_support_arm_report.json")
    repeated_check = load("d6_k6_repeated_two_support_arm_verification.json")
    fan = load("d6_k6_singleton_fan_report.json")
    fan_check = load("d6_k6_singleton_fan_verification.json")
    tight = load("d6_k6_tight_hall_support_report.json")
    tight_check = load("d6_k6_tight_hall_support_verification.json")
    tight_input, tight_bad, tight_left, tight_dependencies = exact_partition(
        tight,
        tight_check,
        schema="d6-k6-tight-hall-support-v1",
        verification_schema="d6-k6-tight-hall-support-verification-v1",
        report_name="d6_k6_tight_hall_support_report.json",
        input_count=756,
        rejection_count=107,
        residue_count=649,
    )
    fan_input, fan_bad, fan_left, fan_dependencies = exact_partition(
        fan,
        fan_check,
        schema="d6-k6-singleton-fan-v1",
        verification_schema="d6-k6-singleton-fan-verification-v1",
        report_name="d6_k6_singleton_fan_report.json",
        input_count=649,
        rejection_count=15,
        residue_count=634,
    )
    arm_input, arm_bad, k6_indices, arm_dependencies = exact_partition(
        repeated,
        repeated_check,
        schema="d6-k6-repeated-two-support-arm-v1",
        verification_schema="d6-k6-repeated-two-support-arm-verification-v1",
        report_name="d6_k6_repeated_two_support_arm_report.json",
        input_count=634,
        rejection_count=9,
        residue_count=625,
    )
    require(
        tight_input == v5_k6_indices
        and tight_left == fan_input
        and fan_left == arm_input
        and stable_hash(k6_indices) == EXPECTED_K6_SHA256,
        "K6 exact chain",
    )
    k6_rejected = tight_bad + fan_bad + arm_bad
    require(
        len(k6_rejected) == 131 and len(set(k6_rejected)) == 131,
        "K6 rejection partition",
    )
    checks["K6_exact_chain"] = True

    require(
        manifest.get("schema") == "d6-current-certified-residue-v6"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and manifest.get("class_order") == ["K7", "K6_only"],
        "manifest schema/status",
    )
    embedded_k7 = manifest["classes"]["K7"]["graphs"]
    embedded_k6 = manifest["classes"]["K6_only"]["graphs"]
    expected_k7_graphs = [
        record for record in v5_k7 if record["index"] in set(k7_indices)
    ]
    expected_k6_graphs = [
        record for record in v5_k6 if record["index"] in set(k6_indices)
    ]
    require(embedded_k7 == expected_k7_graphs, "embedded K7 adjacency records")
    require(embedded_k6 == expected_k6_graphs, "embedded K6 adjacency records")
    require(
        [record["index"] for record in embedded_k7] == k7_indices
        and [record["index"] for record in embedded_k6] == k6_indices,
        "embedded class order",
    )
    for class_name, records in (("K7", embedded_k7), ("K6_only", embedded_k6)):
        for record in records:
            rows = validate_rows(
                record.get("adjacency"), f"{class_name} {record.get('index')}"
            )
            require(alpha_at_most_two(rows), "independent triple")
            require(contains_clique(rows, 6), "missing K6")
            if class_name == "K7":
                require(contains_clique(rows, 7), "K7 record lacks K7")
            else:
                require(not contains_clique(rows, 7), "K6-only record has K7")
    checks["adjacency_and_graph_classes"] = True

    k7 = manifest["classes"]["K7"]
    k6 = manifest["classes"]["K6_only"]
    require(
        k7.get("v5_count") == 155
        and k7.get("exact_algebra_rejection_count") == 133
        and k7.get("exact_algebra_rejected_indices") == exact_rejected
        and k7.get("exact_algebra_survivor_indices") == exact22
        and k7.get("interval_certified_killed_indices") == killed
        and k7.get("interval_killed_already_exact_indices") == interval_old
        and k7.get("interval_new_rejected_indices") == interval_new
        and k7.get("union_rejection_count") == 136
        and k7.get("count") == 19
        and k7.get("indices") == k7_indices
        and k7.get("indices_sha256") == EXPECTED_K7_SHA256
        and k7.get("graphs_sha256") == stable_hash(embedded_k7),
        "K7 manifest accounting",
    )
    layers = k6.get("exact_layers")
    require(
        isinstance(layers, list)
        and [layer.get("rejection_count") for layer in layers] == [107, 15, 9]
        and [layer.get("residue_count") for layer in layers] == [649, 634, 625]
        and k6.get("exact_rejection_count") == 131
        and k6.get("exact_rejected_indices") == k6_rejected
        and k6.get("count") == 625
        and k6.get("indices") == k6_indices
        and k6.get("indices_sha256") == EXPECTED_K6_SHA256
        and k6.get("graphs_sha256") == stable_hash(embedded_k6),
        "K6 manifest accounting",
    )
    combined_indices = k7_indices + k6_indices
    combined = manifest.get("combined", {})
    require(
        combined.get("v5_count") == 911
        and combined.get("exact_only_residue_count") == 647
        and combined.get("interval_incremental_rejection_count") == 3
        and combined.get("count") == 644
        and combined.get("class_counts") == {"K7": 19, "K6_only": 625}
        and combined.get("ordered_indices_sha256") == stable_hash(combined_indices)
        and combined.get("sorted_indices_sha256")
        == stable_hash(sorted(combined_indices))
        and combined.get("cross_class_overlap") == 0,
        "combined manifest accounting",
    )
    checks["counts_and_hashes"] = True

    accounting = manifest.get("certificate_accounting", {})
    require(
        accounting.get("exact_algebra_and_graph_logic", {}).get(
            "floating_point_enters_rejection"
        )
        is False
        and accounting.get("exact_algebra_and_graph_logic", {}).get(
            "total_rejections"
        )
        == 264
        and accounting.get("interval", {}).get("incremental_rejections") == 3
        and accounting.get("interval", {}).get(
            "non_KILLED_statuses_used_for_rejection"
        )
        is False
        and accounting.get("interval", {}).get("trust_assumptions") == trust
        and accounting.get("union_rejections_from_v5") == 267,
        "trust-tier accounting",
    )
    semantics = manifest.get("semantics", {})
    require(
        "unconstrained" in semantics.get("nonedges", "")
        and "IEEE-754" in semantics.get("interval_rejection", "")
        and "non-rejection" in semantics.get("residue", ""),
        "manifest semantics",
    )
    require(
        manifest.get("positive_18_control", {}).get("K6_exact_chain", {}).get(
            "passed"
        )
        is True
        and manifest.get("positive_18_control", {}).get(
            "K6_independent_verification", {}
        ).get("passed")
        is True,
        "manifest positive controls",
    )
    checks["trust_semantics_and_controls"] = True

    transitive_expected = dict(
        sorted(
            {
                **exact_upstream,
                **exact_sources,
                **tight_dependencies,
                **fan_dependencies,
                **arm_dependencies,
                tight["certificate_archive"]["path"]: tight[
                    "certificate_archive"
                ]["sha256"],
                fan["certificate_archive"]["path"]: fan[
                    "certificate_archive"
                ]["sha256"],
                repeated["certificate_archive"]["path"]: repeated[
                    "certificate_archive"
                ]["sha256"],
            }.items()
        )
    )
    boundary = manifest.get("source_boundary", {})
    require(
        boundary.get("primary_artifacts") == dict(sorted(PRIMARY.items()))
        and boundary.get("transitive_files") == transitive_expected
        and boundary.get("K7_cross_method_survivors_sha256")
        == EXPECTED_K7_SHA256
        and boundary.get("K6_exact_survivors_sha256") == EXPECTED_K6_SHA256,
        "manifest source boundary",
    )
    check_named_hashes(boundary["transitive_files"], "manifest transitive files")
    checks["transitive_hash_roots"] = True

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require(
        "build_d6_current_residue_manifest_v6" not in imports,
        "checker imports builder",
    )
    require(
        not any(
            name.startswith(("d6_k6_", "d6_k7_", "verify_d6_k6_", "verify_d6_k7_"))
            for name in imports
        ),
        "checker imports a production kernel/verifier",
    )
    checks["import_independence"] = True

    require(all(checks.values()), "a v6 check is false")
    return {
        "schema": "d6-current-certified-residue-v6-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "builder_sha256": EXPECTED_BUILDER_SHA256,
        "counts": {
            "K7": 19,
            "K6_only": 625,
            "combined": 644,
            "exact_rejections_from_v5": 264,
            "interval_incremental_rejections": 3,
        },
        "class_indices_sha256": {
            "K7": EXPECTED_K7_SHA256,
            "K6_only": EXPECTED_K6_SHA256,
        },
        "ordered_residue_indices_sha256": stable_hash(combined_indices),
        "trust_tiers": {
            "exact_floating_point": False,
            "interval_incremental_rejections": 3,
            "interval_trust_assumptions": trust,
        },
        "checks": checks,
    }


def atomic_json(path: Path, value: object) -> None:
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
        "--manifest",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v6.json",
    )
    parser.add_argument("--expected-manifest-sha256")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v6_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.manifest.resolve(), args.expected_manifest_sha256)
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {"status": "PASS", **result["counts"], "output": str(args.output)},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
