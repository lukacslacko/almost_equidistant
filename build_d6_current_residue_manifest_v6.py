#!/usr/bin/env python3
"""Build the cross-method d=6 residue v6 from independently checked roots."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PRIMARY_ARTIFACTS = {
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

EXPECTED_K7_ALGEBRA_SURVIVORS_SHA256 = (
    "cce80d8065f5c97c15091343128f9626447b13e7bd6e7638b59c9dcef5296cf6"
)
EXPECTED_K7_UNION_SURVIVORS = [
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
EXPECTED_K7_UNION_SURVIVORS_SHA256 = (
    "af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290"
)
EXPECTED_K6_SURVIVORS_SHA256 = (
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


def load(name: str) -> dict:
    value = json.loads((ROOT / name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{name} is not a JSON object")
    return value


def verify_named_hashes(mapping: dict, label: str) -> dict[str, str]:
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{label} has no named hash roots")
    normalized = {}
    for name, expected in mapping.items():
        if (
            not isinstance(name, str)
            or Path(name).name != name
            or not isinstance(expected, str)
            or len(expected) != 64
        ):
            raise ValueError(f"{label} contains a malformed hash root")
        observed = sha256(ROOT / name)
        if observed != expected:
            raise ValueError(f"{label} changed: {name} {observed}")
        normalized[name] = expected
    return dict(sorted(normalized.items()))


def verify_primary_artifacts() -> None:
    verify_named_hashes(PRIMARY_ARTIFACTS, "primary artifact")


def validate_partition_report(
    report: dict,
    verification: dict,
    *,
    schema: str,
    verification_schema: str,
    report_name: str,
    input_count: int,
    rejected_count: int,
    survivor_count: int,
) -> tuple[list[int], list[int], list[int]]:
    if (
        report.get("schema") != schema
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != input_count
        or report.get("graphs_rejected") != rejected_count
        or report.get("graphs_surviving") != survivor_count
        or verification.get("schema") != verification_schema
        or verification.get("status") != "PASS"
        or verification.get("report", {}).get("sha256")
        != sha256(ROOT / report_name)
        or verification.get("graphs_recomputed") != input_count
        or verification.get("graphs_rejected") != rejected_count
        or verification.get("graphs_surviving") != survivor_count
    ):
        raise ValueError(f"incomplete exact partition: {report_name}")
    inputs = list(report.get("ordered_input_indices", ()))
    rejected = list(report.get("rejected_indices", ()))
    survivors = list(report.get("ordered_residue_indices", ()))
    if (
        len(inputs) != input_count
        or len(set(inputs)) != input_count
        or len(rejected) != rejected_count
        or len(set(rejected)) != rejected_count
        or len(survivors) != survivor_count
        or len(set(survivors)) != survivor_count
        or any(type(index) is not int for index in inputs + rejected + survivors)
        or [index for index in inputs if index not in set(rejected)] != survivors
        or stable_hash(inputs) != report.get("ordered_input_indices_sha256")
        or stable_hash(rejected) != report.get("rejected_indices_sha256")
        or stable_hash(survivors) != report.get("ordered_residue_indices_sha256")
        or verification.get("rejected_indices_sha256") != stable_hash(rejected)
        or verification.get("ordered_residue_indices_sha256")
        != stable_hash(survivors)
    ):
        raise ValueError(f"malformed exact partition: {report_name}")
    dependencies = verify_named_hashes(
        report.get("dependencies", {}), f"{report_name} dependencies"
    )
    archive = report.get("certificate_archive", {})
    archive_name = archive.get("path")
    archive_hash = archive.get("sha256")
    if (
        not isinstance(archive_name, str)
        or Path(archive_name).name != archive_name
        or sha256(ROOT / archive_name) != archive_hash
        or archive.get("rejected_graphs") != rejected_count
    ):
        raise ValueError(f"bad certificate archive: {report_name}")
    source_name = report_name.replace("_report.json", ".py")
    if (
        source_name not in dependencies
        and sha256(ROOT / source_name) != report.get("production_source_sha256")
    ):
        raise ValueError(f"production source changed: {source_name}")
    semantics = report.get("semantics", {})
    if (
        "exact" not in semantics.get("arithmetic", "")
        or "unconstrained" not in semantics.get("candidate_nonedges", "")
        or not report.get("positive_18_control", {}).get("passed")
        or not verification.get("positive_18_control", {}).get("passed")
    ):
        raise ValueError(f"unsafe semantics or control: {report_name}")
    return inputs, rejected, survivors


def build_manifest() -> dict:
    verify_primary_artifacts()
    v5 = load("d6_current_residue_manifest_v5.json")
    v5_verification = load("d6_current_residue_manifest_v5_verification.json")
    if (
        v5.get("schema") != "d6-current-exact-residue-v5"
        or v5.get("status") != "COMPLETE_EXACT_FILTER_UNION"
        or v5_verification.get("status") != "PASS"
        or v5_verification.get("manifest", {}).get("sha256")
        != PRIMARY_ARTIFACTS["d6_current_residue_manifest_v5.json"]
        or v5.get("combined", {}).get("count") != 911
    ):
        raise ValueError("v5 residue boundary is not independently complete")
    v5_k7_graphs = list(v5["classes"]["K7"]["graphs"])
    v5_k6_graphs = list(v5["classes"]["K6_only"]["graphs"])
    v5_k7_indices = [int(record["index"]) for record in v5_k7_graphs]
    v5_k6_indices = [int(record["index"]) for record in v5_k6_graphs]
    if len(v5_k7_indices) != 155 or len(v5_k6_indices) != 756:
        raise ValueError("v5 class populations changed")

    # Exact K7 algebraic chain through the independently replayed one-free
    # conjunction.  Its upstream maps bind the double-pin and full-pin roots.
    k7_exact = load("d6_k7_one_free_conjunction_report.json")
    k7_exact_verification = load(
        "d6_k7_one_free_conjunction_verification.json"
    )
    exact_upstream = verify_named_hashes(
        k7_exact.get("upstream_artifact_sha256", {}),
        "K7 exact upstream",
    )
    exact_sources = verify_named_hashes(
        k7_exact.get("source_sha256", {}), "K7 exact production source"
    )
    if (
        k7_exact.get("kind") != "d6_k7_one_free_edge_seed_conjunction"
        or k7_exact.get("status") != "COMPLETE"
        or k7_exact_verification.get("kind")
        != "d6_k7_one_free_edge_seed_conjunction_verification"
        or k7_exact_verification.get("status") != "PASS"
        or k7_exact_verification.get("report", {}).get("sha256")
        != PRIMARY_ARTIFACTS["d6_k7_one_free_conjunction_report.json"]
        or k7_exact_verification.get("upstream_artifact_sha256")
        != exact_upstream
    ):
        raise ValueError("K7 exact one-free boundary is incomplete")
    base_report = load("d6_k7_double_pin_conjunction_report.json")
    base_verification = load("d6_k7_double_pin_conjunction_verification.json")
    full_report = load("d6_k7_full_pin_increment_report.json")
    full_verification = load("d6_k7_full_pin_increment_verification.json")
    if (
        base_report.get("input", {}).get("graphs") != 155
        or base_report.get("input", {}).get("indices_sha256")
        != stable_hash(v5_k7_indices)
        or base_report.get("input", {}).get("graphs_sha256")
        != stable_hash(v5_k7_graphs)
        or base_verification.get("status") != "PASS"
        or base_verification.get("report", {}).get("sha256")
        != exact_upstream["d6_k7_double_pin_conjunction_report.json"]
        or full_report.get("status") != "COMPLETE"
        or full_verification.get("status") != "PASS"
        or full_verification.get("report", {}).get("sha256")
        != exact_upstream["d6_k7_full_pin_increment_report.json"]
        or k7_exact.get("input_indices") != full_report.get("exact_survivors")
    ):
        raise ValueError("K7 exact transitive roots do not begin at v5")
    k7_exact_survivors = list(
        k7_exact["summary"]["survivors_by_pipeline"]["combined"]
    )
    if (
        len(k7_exact_survivors) != 22
        or stable_hash(k7_exact_survivors)
        != EXPECTED_K7_ALGEBRA_SURVIVORS_SHA256
        or k7_exact_verification.get("exact_survivors_sha256")
        != EXPECTED_K7_ALGEBRA_SURVIVORS_SHA256
    ):
        raise ValueError("K7 exact algebraic residue differs from 22-list")
    exact_survivor_set = set(k7_exact_survivors)
    k7_exact_rejected = [
        index for index in v5_k7_indices if index not in exact_survivor_set
    ]

    # Separately trusted interval certificates.  Only KILLED enters the union.
    interval_report = load("d6_interval_v5_k7_cap100000_report.json")
    interval_verification = load(
        "d6_interval_v5_k7_cap100000_verification.json"
    )
    if (
        interval_report.get("schema") != 2
        or interval_report.get("completed_graphs") != 155
        or interval_report.get("total_graphs") != 155
        or interval_report.get("certified_killed") != 24
        or interval_verification.get("schema")
        != "d6-interval-v5-verification-v1"
        or interval_verification.get("status") != "PASS"
        or interval_verification.get("report", {}).get("sha256")
        != PRIMARY_ARTIFACTS["d6_interval_v5_k7_cap100000_report.json"]
        or interval_verification.get("selection", {}).get("graphs") != 155
        or interval_verification.get("selection", {}).get("indices_sha256")
        != stable_hash(v5_k7_indices)
        or interval_verification.get("kernel_controls", {}).get("status")
        != "PASS"
        or interval_verification.get("witness_replay", {}).get("status")
        != "PASS"
    ):
        raise ValueError("interval K7 certificate boundary is incomplete")
    interval_killed = list(interval_verification.get("killed_indices", ()))
    if (
        len(interval_killed) != 24
        or len(set(interval_killed)) != 24
        or not set(interval_killed) <= set(v5_k7_indices)
        or interval_report.get("status_counts", {}).get("KILLED") != 24
        or interval_verification.get("status_counts", {}).get("KILLED") != 24
    ):
        raise ValueError("interval KILLED list is malformed")
    interval_new = [
        index for index in k7_exact_survivors if index in set(interval_killed)
    ]
    interval_already_exact = [
        index for index in interval_killed if index not in exact_survivor_set
    ]
    k7_union_survivors = [
        index for index in k7_exact_survivors if index not in set(interval_killed)
    ]
    if (
        interval_new != [423661, 424226, 3936176]
        or len(interval_already_exact) != 21
        or k7_union_survivors != EXPECTED_K7_UNION_SURVIVORS
        or stable_hash(k7_union_survivors)
        != EXPECTED_K7_UNION_SURVIVORS_SHA256
    ):
        raise ValueError("cross-method K7 union differs from frozen 19-list")
    k7_union_set = set(k7_union_survivors)
    k7_graphs = [
        record for record in v5_k7_graphs if record["index"] in k7_union_set
    ]
    if [record["index"] for record in k7_graphs] != k7_union_survivors:
        raise ValueError("K7 adjacency reconstruction changed order")

    # Exact K6 chain: v5 -> tight Hall -> singleton fan -> repeated arm.
    repeated = load("d6_k6_repeated_two_support_arm_report.json")
    repeated_verification = load(
        "d6_k6_repeated_two_support_arm_verification.json"
    )
    repeated_dependencies = verify_named_hashes(
        repeated.get("dependencies", {}), "repeated-arm dependencies"
    )
    singleton = load("d6_k6_singleton_fan_report.json")
    singleton_verification = load("d6_k6_singleton_fan_verification.json")
    singleton_dependencies = verify_named_hashes(
        singleton.get("dependencies", {}), "singleton-fan dependencies"
    )
    tight = load("d6_k6_tight_hall_support_report.json")
    tight_verification = load("d6_k6_tight_hall_support_verification.json")
    tight_dependencies = verify_named_hashes(
        tight.get("dependencies", {}), "tight-Hall dependencies"
    )
    tight_input, tight_rejected, tight_residue = validate_partition_report(
        tight,
        tight_verification,
        schema="d6-k6-tight-hall-support-v1",
        verification_schema="d6-k6-tight-hall-support-verification-v1",
        report_name="d6_k6_tight_hall_support_report.json",
        input_count=756,
        rejected_count=107,
        survivor_count=649,
    )
    fan_input, fan_rejected, fan_residue = validate_partition_report(
        singleton,
        singleton_verification,
        schema="d6-k6-singleton-fan-v1",
        verification_schema="d6-k6-singleton-fan-verification-v1",
        report_name="d6_k6_singleton_fan_report.json",
        input_count=649,
        rejected_count=15,
        survivor_count=634,
    )
    arm_input, arm_rejected, k6_survivors = validate_partition_report(
        repeated,
        repeated_verification,
        schema="d6-k6-repeated-two-support-arm-v1",
        verification_schema="d6-k6-repeated-two-support-arm-verification-v1",
        report_name="d6_k6_repeated_two_support_arm_report.json",
        input_count=634,
        rejected_count=9,
        survivor_count=625,
    )
    if (
        tight_input != v5_k6_indices
        or tight_residue != fan_input
        or fan_residue != arm_input
        or stable_hash(k6_survivors) != EXPECTED_K6_SURVIVORS_SHA256
        or repeated_verification.get("ordered_residue_indices_sha256")
        != EXPECTED_K6_SURVIVORS_SHA256
    ):
        raise ValueError("exact K6 partition chain does not span v5 to 625")
    k6_rejected = tight_rejected + fan_rejected + arm_rejected
    if len(k6_rejected) != 131 or len(set(k6_rejected)) != 131:
        raise ValueError("K6 exact rejection layers overlap unexpectedly")
    k6_survivor_set = set(k6_survivors)
    k6_graphs = [
        record for record in v5_k6_graphs if record["index"] in k6_survivor_set
    ]
    if [record["index"] for record in k6_graphs] != k6_survivors:
        raise ValueError("K6 adjacency reconstruction changed order")

    combined_indices = k7_union_survivors + k6_survivors
    tagged_indices = [
        {"class": class_name, "index": record["index"]}
        for class_name, records in (("K7", k7_graphs), ("K6_only", k6_graphs))
        for record in records
    ]
    tagged_graphs = [
        {"class": class_name, **record}
        for class_name, records in (("K7", k7_graphs), ("K6_only", k6_graphs))
        for record in records
    ]
    transitive_files = dict(
        sorted(
            {
                **exact_upstream,
                **exact_sources,
                **tight_dependencies,
                **singleton_dependencies,
                **repeated_dependencies,
                tight["certificate_archive"]["path"]: tight[
                    "certificate_archive"
                ]["sha256"],
                singleton["certificate_archive"]["path"]: singleton[
                    "certificate_archive"
                ]["sha256"],
                repeated["certificate_archive"]["path"]: repeated[
                    "certificate_archive"
                ]["sha256"],
            }.items()
        )
    )
    trust_assumptions = interval_verification.get("trust_assumptions", {})
    if set(trust_assumptions) != {
        "basic_arithmetic",
        "candidate_nonedges",
        "transcendentals",
    }:
        raise ValueError("interval trust assumptions are incomplete")

    return {
        "schema": "d6-current-certified-residue-v6",
        "status": "COMPLETE_MIXED_CERTIFICATE_UNION",
        "description": (
            "The verified v5 residue after exact K7/K6 algebraic filters "
            "and the separately trusted K7 interval KILLED certificates. "
            "Survival is not realizability."
        ),
        "class_order": ["K7", "K6_only"],
        "classes": {
            "K7": {
                "v5_count": 155,
                "exact_algebra_rejection_count": len(k7_exact_rejected),
                "exact_algebra_rejected_indices": k7_exact_rejected,
                "exact_algebra_rejected_indices_sha256": stable_hash(
                    k7_exact_rejected
                ),
                "exact_algebra_survivor_count": 22,
                "exact_algebra_survivor_indices": k7_exact_survivors,
                "exact_algebra_survivor_indices_sha256": (
                    EXPECTED_K7_ALGEBRA_SURVIVORS_SHA256
                ),
                "interval_certified_killed_count": 24,
                "interval_certified_killed_indices": interval_killed,
                "interval_certified_killed_indices_sha256": stable_hash(
                    interval_killed
                ),
                "interval_killed_already_exact_count": 21,
                "interval_killed_already_exact_indices": interval_already_exact,
                "interval_new_rejection_count": 3,
                "interval_new_rejected_indices": interval_new,
                "interval_new_rejected_indices_sha256": stable_hash(interval_new),
                "union_rejection_count": 136,
                "count": 19,
                "indices": k7_union_survivors,
                "indices_sha256": EXPECTED_K7_UNION_SURVIVORS_SHA256,
                "graphs": k7_graphs,
                "graphs_sha256": stable_hash(k7_graphs),
            },
            "K6_only": {
                "v5_count": 756,
                "exact_layers": [
                    {
                        "name": "tight_hall_support",
                        "input_count": 756,
                        "rejection_count": 107,
                        "rejected_indices": tight_rejected,
                        "rejected_indices_sha256": stable_hash(tight_rejected),
                        "residue_count": 649,
                        "residue_indices_sha256": stable_hash(tight_residue),
                    },
                    {
                        "name": "singleton_fan",
                        "input_count": 649,
                        "rejection_count": 15,
                        "rejected_indices": fan_rejected,
                        "rejected_indices_sha256": stable_hash(fan_rejected),
                        "residue_count": 634,
                        "residue_indices_sha256": stable_hash(fan_residue),
                    },
                    {
                        "name": "repeated_two_support_arm",
                        "input_count": 634,
                        "rejection_count": 9,
                        "rejected_indices": arm_rejected,
                        "rejected_indices_sha256": stable_hash(arm_rejected),
                        "residue_count": 625,
                        "residue_indices_sha256": stable_hash(k6_survivors),
                    },
                ],
                "exact_rejection_count": 131,
                "exact_rejected_indices": k6_rejected,
                "exact_rejected_indices_sha256": stable_hash(k6_rejected),
                "count": 625,
                "indices": k6_survivors,
                "indices_sha256": EXPECTED_K6_SURVIVORS_SHA256,
                "graphs": k6_graphs,
                "graphs_sha256": stable_hash(k6_graphs),
            },
        },
        "combined": {
            "v5_count": 911,
            "exact_only_residue_count": 647,
            "interval_incremental_rejection_count": 3,
            "count": 644,
            "class_counts": {"K7": 19, "K6_only": 625},
            "ordered_indices_sha256": stable_hash(combined_indices),
            "sorted_indices_sha256": stable_hash(sorted(combined_indices)),
            "ordered_class_index_records_sha256": stable_hash(tagged_indices),
            "ordered_class_graph_records_sha256": stable_hash(tagged_graphs),
            "cross_class_overlap": len(
                set(k7_union_survivors) & set(k6_survivors)
            ),
        },
        "certificate_accounting": {
            "exact_algebra_and_graph_logic": {
                "K7_rejections": 133,
                "K6_only_rejections": 131,
                "total_rejections": 264,
                "floating_point_enters_rejection": False,
            },
            "interval": {
                "K7_KILLED_certificates": 24,
                "already_exactly_rejected": 21,
                "incremental_rejections": 3,
                "non_KILLED_statuses_used_for_rejection": False,
                "trust_assumptions": trust_assumptions,
            },
            "union_rejections_from_v5": 267,
        },
        "positive_18_control": {
            "v5": v5["positive_18_control"],
            "K6_exact_chain": repeated["positive_18_control"],
            "K6_independent_verification": repeated_verification[
                "positive_18_control"
            ],
            "interval_kernel_controls": interval_verification[
                "kernel_controls"
            ],
        },
        "source_boundary": {
            "primary_artifacts": dict(sorted(PRIMARY_ARTIFACTS.items())),
            "transitive_files": transitive_files,
            "v5_K7_indices_sha256": stable_hash(v5_k7_indices),
            "v5_K6_indices_sha256": stable_hash(v5_k6_indices),
            "K7_exact_algebra_survivors_sha256": (
                EXPECTED_K7_ALGEBRA_SURVIVORS_SHA256
            ),
            "K7_interval_selection_sha256": interval_verification[
                "selection"
            ]["indices_sha256"],
            "K7_cross_method_survivors_sha256": (
                EXPECTED_K7_UNION_SURVIVORS_SHA256
            ),
            "K6_exact_survivors_sha256": EXPECTED_K6_SURVIVORS_SHA256,
        },
        "semantics": {
            "edges": "required unit distances",
            "nonedges": "unconstrained and may also have distance one",
            "points": "distinct",
            "exact_rejection": (
                "sound exact algebraic, rational, integer, and graph-logic "
                "certificates with independently replayed quantifiers"
            ),
            "interval_rejection": (
                "only verified KILLED records, conditional on the stated "
                "IEEE-754 and macOS libm outward-rounding assumptions"
            ),
            "residue": "certificate-union non-rejection only",
        },
        "nonclaims": [
            "No v6 residue graph is claimed realizable.",
            "The 644-graph boundary is not yet a proof that f(6)=18.",
            "ABORT, UNRESOLVED, and INFRA_ERROR interval statuses reject nothing.",
            "No candidate nonedge is constrained to be non-unit.",
        ],
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
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v6.json",
    )
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output.resolve(), manifest)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
                "K7": 19,
                "K6_only": 625,
                "combined": 644,
                "K7_indices_sha256": EXPECTED_K7_UNION_SURVIVORS_SHA256,
                "K6_indices_sha256": EXPECTED_K6_SURVIVORS_SHA256,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
