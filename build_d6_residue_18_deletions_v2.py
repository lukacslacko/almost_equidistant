#!/usr/bin/env python3
"""Build the rooted 18-deletion corpus of the exact d=6 v5 residue.

The v5 parent set is an exact subset of v4.  The verified v1 deletion corpus
already contains every labeled deletion occurrence of every v4 parent.  This
builder therefore performs a lossless occurrence-level subset operation,
then recomputes all group counts, hashes, rooted attachment records, and
standard18 match summaries.  The independent v2 checker reconstructs and
canonicalizes all deletions afresh from the v5 parent adjacencies.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
INPUT_NAME = "d6_current_residue_manifest_v5.json"
INPUT_VERIFICATION_NAME = "d6_current_residue_manifest_v5_verification.json"
V1_NAME = "d6_residue_18_deletions.json"
V1_VERIFICATION_NAME = "d6_residue_18_deletions_verification.json"
EXPECTED_FILES = {
    INPUT_NAME: "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5",
    INPUT_VERIFICATION_NAME: "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46",
    "build_d6_current_residue_manifest_v5.py": "2062abce9a4b9d075b794c03ad81a7f6403b43356b56691623b9fc9288a4d0a3",
    "verify_d6_current_residue_manifest_v5.py": "7f3e5c0c3b74793da50f3f7ca961e0d51c70f10998f0b3ea7e87211760e995fd",
    V1_NAME: "9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732",
    V1_VERIFICATION_NAME: "50ab5d88441b2869e60a05f3e1485db0b25008c3a816960a20a3bb21c1c4defd",
    "build_d6_residue_18_deletions.py": "f44ef42f8da3e6ba4b4c902b2db11468d781b8d9157eab3081a9df43674a08fd",
    "verify_d6_residue_18_deletions.py": "5ac2baab623a49de603fe0d5d30581402e15a1753d48a02b8ac2ce3c4ec251d4",
    "test_d6_residue_18_deletions.py": "ce96f9932eb81d9452c0feb0b45718b84e22f5c2e9cd4c8290e8987f787a0559",
}

EXPECTED_COUNTS = {
    "parents": 911,
    "deletion_occurrences": 17_309,
    "unique_deletions": 11_975,
    "unique_K7_parent_deletions": 2_288,
    "unique_K6_only_parent_deletions": 9_749,
    "cross_class_unique_overlap": 62,
    "standard_support_types": 7,
    "standard_support_labeled": 1_960,
    "standard_support_deletion_occurrences": 26,
    "standard_support_parent_graphs": 14,
    "standard_compatible_unique_deletions": 14,
    "standard_compatible_deletion_occurrences": 39,
    "standard_compatible_parent_graphs": 16,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load(name: str) -> dict:
    value = json.loads((ROOT / name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{name} is not a JSON object")
    return value


def verify_sources() -> None:
    observed = {name: sha256(ROOT / name) for name in EXPECTED_FILES}
    if observed != EXPECTED_FILES:
        raise ValueError(f"a pinned deletion-v2 source changed: {observed}")


def histogram(values: Iterable[int]) -> dict[str, int]:
    return {str(key): count for key, count in sorted(Counter(values).items())}


def occurrence_key(occurrence: dict) -> tuple[str, int, int]:
    class_name = occurrence.get("parent_class")
    parent_index = occurrence.get("parent_index")
    deleted = occurrence.get("deleted_vertex")
    if (
        class_name not in {"K7", "K6_only"}
        or type(parent_index) is not int
        or type(deleted) is not int
        or not 0 <= deleted < 19
    ):
        raise ValueError("malformed rooted deletion occurrence")
    return class_name, parent_index, deleted


def build_manifest() -> dict:
    verify_sources()
    source = load(INPUT_NAME)
    source_check = load(INPUT_VERIFICATION_NAME)
    v1 = load(V1_NAME)
    v1_check = load(V1_VERIFICATION_NAME)
    if (
        source.get("schema") != "d6-current-exact-residue-v5"
        or source.get("status") != "COMPLETE_EXACT_FILTER_UNION"
        or source.get("combined", {}).get("count") != 911
        or source_check.get("status") != "PASS"
        or source_check.get("manifest", {}).get("sha256")
        != EXPECTED_FILES[INPUT_NAME]
    ):
        raise ValueError("v5 parent boundary is incomplete")
    if (
        v1.get("schema") != "d6-residue-18-deletion-manifest-v1"
        or v1_check.get("status") != "PASS"
        or v1_check.get("manifest", {}).get("sha256")
        != EXPECTED_FILES[V1_NAME]
    ):
        raise ValueError("v1 deletion corpus boundary is incomplete")

    parents = []
    for class_name in ("K7", "K6_only"):
        for record in source["classes"][class_name]["graphs"]:
            index = record.get("index")
            if type(index) is not int:
                raise ValueError("v5 parent index is not an integer")
            parents.append((class_name, index))
    if len(parents) != 911 or len(set(parents)) != 911:
        raise ValueError("v5 parent list is not 911 distinct tagged records")

    old_by_code = {}
    old_occurrence_index: dict[tuple[str, int, int], tuple[str, dict]] = {}
    for record in v1["unique_deletions"]:
        code = record.get("canonical_graph6")
        if not isinstance(code, str) or code in old_by_code:
            raise ValueError("v1 deletion canonical code is malformed or duplicated")
        old_by_code[code] = record
        for occurrence in record.get("occurrences", []):
            key = occurrence_key(occurrence)
            if key in old_occurrence_index:
                raise ValueError("v1 rooted occurrence is duplicated")
            old_occurrence_index[key] = (code, occurrence)
    if len(old_occurrence_index) != 18_240:
        raise ValueError("v1 corpus does not contain 18,240 rooted occurrences")

    occurrences = []
    deletion_codes = []
    grouped: dict[str, list[dict]] = defaultdict(list)
    for class_name, parent_index in parents:
        for deleted in range(19):
            key = (class_name, parent_index, deleted)
            if key not in old_occurrence_index:
                raise ValueError(f"v1 corpus lacks v5 occurrence {key}")
            code, occurrence = old_occurrence_index[key]
            copied = copy.deepcopy(occurrence)
            occurrences.append(copied)
            deletion_codes.append(code)
            grouped[code].append(copied)
    if len(occurrences) != 17_309:
        raise AssertionError("911 parents did not produce 17,309 occurrences")

    class_codes: dict[str, set[str]] = {"K7": set(), "K6_only": set()}
    unique_records = []
    standard_occurrences = []
    standard_parents = set()
    compatible_occurrences = []
    compatible_parents = set()
    compatible_unique = 0
    for code in sorted(grouped):
        base = old_by_code[code]
        group = grouped[code]
        classes = sorted({occurrence["parent_class"] for occurrence in group})
        for class_name in classes:
            class_codes[class_name].add(code)
        record = copy.deepcopy(base)
        record["parent_classes"] = classes
        record["occurrence_count"] = len(group)
        record["parent_count"] = len(
            {
                (occurrence["parent_class"], occurrence["parent_index"])
                for occurrence in group
            }
        )
        record["occurrences"] = group
        if record.get("standard18_support_type") is not None:
            standard_occurrences.extend(group)
            standard_parents.update(
                (occurrence["parent_class"], occurrence["parent_index"])
                for occurrence in group
            )
        if record.get("standard18_compatible"):
            compatible_unique += 1
            compatible_occurrences.extend(group)
            compatible_parents.update(
                (occurrence["parent_class"], occurrence["parent_index"])
                for occurrence in group
            )
        unique_records.append(record)

    standard_types = copy.deepcopy(v1["standard18"]["minimal_support_types"])
    observed_counts = {
        "parents": len(parents),
        "deletion_occurrences": len(occurrences),
        "unique_deletions": len(unique_records),
        "unique_K7_parent_deletions": len(class_codes["K7"]),
        "unique_K6_only_parent_deletions": len(class_codes["K6_only"]),
        "cross_class_unique_overlap": len(
            class_codes["K7"] & class_codes["K6_only"]
        ),
        "standard_support_types": len(standard_types),
        "standard_support_labeled": sum(
            record["labeled_supports"] for record in standard_types
        ),
        "standard_support_deletion_occurrences": len(standard_occurrences),
        "standard_support_parent_graphs": len(standard_parents),
        "standard_compatible_unique_deletions": compatible_unique,
        "standard_compatible_deletion_occurrences": len(compatible_occurrences),
        "standard_compatible_parent_graphs": len(compatible_parents),
    }
    if observed_counts != EXPECTED_COUNTS:
        raise AssertionError(
            f"deletion v2 count drift: {observed_counts!r} != {EXPECTED_COUNTS!r}"
        )

    parent_indices = [index for _, index in parents]
    unique_codes = [record["canonical_graph6"] for record in unique_records]
    edge_by_code = {
        record["canonical_graph6"]: int(record["edges"])
        for record in unique_records
    }
    summary = {
        **observed_counts,
        "parent_class_counts": {
            "K7": sum(class_name == "K7" for class_name, _ in parents),
            "K6_only": sum(class_name == "K6_only" for class_name, _ in parents),
        },
        "occurrence_edge_histogram": histogram(
            edge_by_code[code] for code in deletion_codes
        ),
        "unique_edge_histogram": histogram(
            int(record["edges"]) for record in unique_records
        ),
        "unique_clique_classes": {
            "contains_K7": sum(record["contains_K7"] for record in unique_records),
            "K6_without_K7": sum(
                record["contains_K6"] and not record["contains_K7"]
                for record in unique_records
            ),
            "no_K6": sum(
                not record["contains_K6"] for record in unique_records
            ),
        },
        "parent_indices_sha256": stable_hash(parent_indices),
        "ordered_occurrences_sha256": stable_hash(occurrences),
        "ordered_deletion_canonical_graph6_sha256": stable_hash(deletion_codes),
        "unique_canonical_graph6_sha256": stable_hash(unique_codes),
        "unique_records_sha256": stable_hash(unique_records),
        "standard_support_types_sha256": stable_hash(standard_types),
        "standard_support_occurrences_sha256": stable_hash(standard_occurrences),
        "standard_compatible_occurrences_sha256": stable_hash(
            compatible_occurrences
        ),
    }
    return {
        "schema": "d6-residue-18-deletion-manifest-v2",
        "source": {
            "path": INPUT_NAME,
            "sha256": EXPECTED_FILES[INPUT_NAME],
            "schema": source["schema"],
            "combined_count": source["combined"]["count"],
            "ordered_indices_sha256": source["combined"][
                "ordered_indices_sha256"
            ],
        },
        "source_verification": {
            "path": INPUT_VERIFICATION_NAME,
            "sha256": EXPECTED_FILES[INPUT_VERIFICATION_NAME],
            "status": "PASS",
        },
        "v1_derivation": {
            "path": V1_NAME,
            "sha256": EXPECTED_FILES[V1_NAME],
            "verification_path": V1_VERIFICATION_NAME,
            "verification_sha256": EXPECTED_FILES[V1_VERIFICATION_NAME],
            "method": (
                "lossless filtering of labeled rooted occurrences by the exact "
                "v5 parent set; canonical classes and standard18 witnesses are "
                "independently rebuilt by the v2 checker"
            ),
        },
        "canonicalizer": {
            **v1["canonicalizer"],
            "role_in_v2_builder": (
                "canonical labels are inherited at the occurrence level from "
                "the verified v1 superset; no geometric decision is made"
            ),
        },
        "summary": summary,
        "standard18": {
            "full_unit_graph": copy.deepcopy(
                v1["standard18"]["full_unit_graph"]
            ),
            "full_unit_graph_sha256": v1["standard18"][
                "full_unit_graph_sha256"
            ],
            "minimal_support_types": standard_types,
            "matched_deletion_occurrences": standard_occurrences,
            "matched_parent_records": [
                {"parent_class": class_name, "parent_index": parent_index}
                for class_name, parent_index in sorted(standard_parents)
            ],
            "compatible_deletion_occurrences": compatible_occurrences,
            "compatible_parent_records": [
                {"parent_class": class_name, "parent_index": parent_index}
                for class_name, parent_index in sorted(compatible_parents)
            ],
            "semantics": v1["standard18"]["semantics"],
        },
        "unique_deletions": unique_records,
        "semantics": copy.deepcopy(v1["semantics"]),
        "nonclaims": copy.deepcopy(v1["nonclaims"]),
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
        "--output", type=Path, default=ROOT / "d6_residue_18_deletions_v2.json"
    )
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output.resolve(), manifest)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
                **manifest["summary"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
