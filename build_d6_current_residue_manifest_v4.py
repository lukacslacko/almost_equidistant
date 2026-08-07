#!/usr/bin/env python3
"""Build the exact d=6 residue after the K6 empty-support-budget layer."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EMPTY_SUPPORT_SOURCE_COMMIT = "b54b69571b5ff87c7586aff13e5c69a24ece7801"
EXPECTED_FILES = {
    "d6_current_residue_manifest_v3.json": "adbb28c93b9262eeb515e491821c562c00670499f133558ea0821e03593e10fd",
    "d6_current_residue_manifest_v3_verification.json": "e8778b550a5bbe27ab21170dfbb7c346fff83f898054b65074f14a66f9ded10a",
    "d6_k6_empty_support_budget.py": "284b8c3bda4d2a43581665ade3d52d9295b29d460418454e3a760897cabc5fac",
    "verify_d6_k6_empty_support_budget.py": "5e4929d497d7b78b7b46b399deaafa66532029b228aa280bde2f23910d04f969",
    "test_d6_k6_empty_support_budget.py": "099689f870c163ea51f5d843aa8bf21acaf6ef857910e9627a625fcc6fcbc371",
    "d6_k6_empty_support_budget_report.json": "1860dbe69ae74b55203ea283cf83afff929a4b1f5cbfcd30e09f0138688eba9f",
    "d6_k6_empty_support_budget_certificates.json": "c6fcb7ef66bccc3319fe0a979c5fd63d9f6fd9535261c5c9bfb87e8d210361d4",
    "d6_k6_empty_support_budget_checkpoint.json": "da13eecacc47559f617924a593692f42ed6d05b49d4e5deeed7b0987ebbd0bd0",
    "d6_k6_empty_support_budget_verification.json": "873c8b6babe183757d4e97b0f0c1c1a942b31f180c88eeef1c9ddf80a199276f",
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
        raise ValueError(f"{name} is not an object")
    return value


def verify_sources() -> None:
    observed = {name: sha256(ROOT / name) for name in EXPECTED_FILES}
    if observed != EXPECTED_FILES:
        raise ValueError("a pinned v4 source or result artifact changed")


def build_manifest() -> dict:
    verify_sources()
    v3 = load("d6_current_residue_manifest_v3.json")
    v3_verification = load("d6_current_residue_manifest_v3_verification.json")
    report = load("d6_k6_empty_support_budget_report.json")
    certificates = load("d6_k6_empty_support_budget_certificates.json")
    verification = load("d6_k6_empty_support_budget_verification.json")

    if v3.get("schema") != "d6-current-exact-residue-v3":
        raise ValueError("unexpected v3 schema")
    if v3_verification.get("status") != "PASS":
        raise ValueError("v3 independent verification did not pass")
    if report.get("status") != "COMPLETE" or report.get("input_graphs") != 822:
        raise ValueError("empty-support production report is incomplete")
    if certificates.get("status") != "COMPLETE":
        raise ValueError("empty-support certificate archive is incomplete")
    if verification.get("status") != "PASS":
        raise ValueError("empty-support independent verification did not pass")
    if verification.get("report_sha256") != EXPECTED_FILES[
        "d6_k6_empty_support_budget_report.json"
    ]:
        raise ValueError("verification does not bind the production report")

    rejected = report.get("rejected_indices")
    if (
        not isinstance(rejected, list)
        or len(rejected) != 17
        or any(type(index) is not int for index in rejected)
        or len(set(rejected)) != 17
    ):
        raise ValueError("invalid empty-support rejection list")
    if rejected != verification.get("rejected_indices") or rejected != certificates.get(
        "rejected_indices"
    ):
        raise ValueError("empty-support result artifacts disagree")
    rejection_set = set(rejected)

    k7_graphs = list(v3["classes"]["K7"]["graphs"])
    k6_before = list(v3["classes"]["K6_only"]["graphs"])
    if [record["index"] for record in k6_before] != report.get(
        "ordered_input_indices"
    ):
        raise ValueError("empty-support input ordering is not the v3 K6 class")
    k6_graphs = [
        record for record in k6_before if record["index"] not in rejection_set
    ]
    k6_indices = [record["index"] for record in k6_graphs]
    if k6_indices != report.get("ordered_residue_indices"):
        raise ValueError("post-empty-support residue ordering mismatch")
    if len(k7_graphs) != 155 or len(k6_graphs) != 805:
        raise ValueError("unexpected v4 class counts")

    k7_indices = [record["index"] for record in k7_graphs]
    combined_indices = k7_indices + k6_indices
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

    return {
        "schema": "d6-current-exact-residue-v4",
        "status": "COMPLETE_EXACT_FILTER_UNION",
        "description": (
            "The v3 exact residue after the independently verified K6 "
            "empty-defect support-budget layer. Survival is not realizability."
        ),
        "class_order": ["K7", "K6_only"],
        "classes": {
            "K7": {
                "count": 155,
                "indices": k7_indices,
                "indices_sha256": stable_hash(k7_indices),
                "graphs": k7_graphs,
                "graphs_sha256": stable_hash(k7_graphs),
                "change_from_v3": 0,
            },
            "K6_only": {
                "v3_count": 822,
                "empty_support_rejection_count": 17,
                "empty_support_rejected_indices": rejected,
                "empty_support_rejected_indices_sha256": stable_hash(rejected),
                "count": 805,
                "indices": k6_indices,
                "indices_sha256": stable_hash(k6_indices),
                "graphs": k6_graphs,
                "graphs_sha256": stable_hash(k6_graphs),
            },
        },
        "combined": {
            "count": 960,
            "class_counts": {"K7": 155, "K6_only": 805},
            "ordered_indices_sha256": stable_hash(combined_indices),
            "sorted_indices_sha256": stable_hash(sorted(combined_indices)),
            "ordered_class_index_records_sha256": stable_hash(tagged_indices),
            "ordered_class_graph_records_sha256": stable_hash(tagged_graphs),
            "cross_class_overlap": len(set(k7_indices) & set(k6_indices)),
        },
        "positive_18_control": v3["positive_18_control"],
        "source_boundary": {
            "empty_support_source_commit": EMPTY_SUPPORT_SOURCE_COMMIT,
            "files": EXPECTED_FILES,
            "v3_ordered_indices_sha256": v3["combined"][
                "ordered_indices_sha256"
            ],
            "empty_support_input_indices_sha256": report[
                "input_indices_sha256"
            ],
            "empty_support_residue_indices_sha256": report[
                "ordered_residue_indices_sha256"
            ],
        },
        "semantics": {
            "edges": "required unit distances",
            "nonedges": "unconstrained and may also have distance one",
            "points": "distinct",
            "rejection": (
                "one required K6 seed exhausts every exact zero-factor, "
                "lightlike, generic-orientation, support, Hall, and globally "
                "permitted empty-defect branch"
            ),
            "residue": "exact filter non-rejection only",
        },
        "nonclaims": [
            "No v4 residue graph is claimed realizable.",
            "The 960-graph boundary is not yet a proof that f(6)=18.",
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
        "--output", type=Path, default=ROOT / "d6_current_residue_manifest_v4.json"
    )
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output, manifest)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output),
                "K7": 155,
                "K6_only": 805,
                "combined": 960,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
