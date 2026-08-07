#!/usr/bin/env python3
"""Build the exact d=6 residue after empty-essential virtual-K7 closure."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPECTED_FILES = {
    "d6_current_residue_manifest_v4.json": "6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4",
    "d6_current_residue_manifest_v4_verification.json": "765eceb6782d131dae8c77c9b94de78735e23a070da051c8fffbd984790e1a41",
    "d6_k6_empty_essential_profile.py": "dab7603b11946635a29c3c37daf44a895db342b2cdddb17c3699f87c182f507a",
    "d6_k6_empty_essential_profile.json": "205c0877982c8de8ad6e5a685e33eff9ee05236c9907e60d488a593792d90781",
    "verify_d6_k6_empty_essential_profile.py": "8eeea41c5ebd04125ad2a33539757229472b2709a8771f56cc6af394c8cf3c2e",
    "d6_k6_empty_essential_profile_verification.json": "c8f73e8a9bed6d7a102ab168eb3842312a916b23dbd577cceb6d089c8832f719",
    "d6_k6_empty_essential_virtual_k7.py": "dbfd44cef783be6849255504518697279d3dfc8dabf1aab72b221697f596a6c5",
    "d6_k6_empty_essential_virtual_k7.json": "63207da634a56e8fa46f18f87e1ed8b86007c62b51c49c5584eb80403f2b50ba",
    "verify_d6_k6_empty_essential_virtual_k7.py": "83c04d7b6e035f63c9e1d38867dc4be2594397709b318cb8379a4ba42b443608",
    "d6_k6_empty_essential_virtual_k7_verification.json": "d64c36856e1288c8f64ea7918bc76640bde260b9b731e881f9c277c1a1e43485",
    "test_d6_k6_empty_essential_virtual_k7.py": "c885a0150bed8bb1ab5cd9b52cf8d6c257959efeef0454a75e578e8b35d79857",
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


def verify_sources() -> None:
    observed = {name: sha256(ROOT / name) for name in EXPECTED_FILES}
    if observed != EXPECTED_FILES:
        raise ValueError(f"a pinned v5 source artifact changed: {observed}")


def build_manifest() -> dict:
    verify_sources()
    v4 = load("d6_current_residue_manifest_v4.json")
    v4_check = load("d6_current_residue_manifest_v4_verification.json")
    report = load("d6_k6_empty_essential_virtual_k7.json")
    verification = load("d6_k6_empty_essential_virtual_k7_verification.json")

    if (
        v4.get("schema") != "d6-current-exact-residue-v4"
        or v4.get("status") != "COMPLETE_EXACT_FILTER_UNION"
        or v4_check.get("status") != "PASS"
        or v4_check.get("manifest", {}).get("sha256")
        != EXPECTED_FILES["d6_current_residue_manifest_v4.json"]
    ):
        raise ValueError("v4 residue boundary is not independently complete")
    if (
        report.get("schema") != "d6-k6-empty-essential-virtual-k7-v1"
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != 805
        or report.get("target_graphs") != 49
    ):
        raise ValueError("virtual-K7 production report is incomplete")
    if (
        verification.get("schema")
        != "d6-k6-empty-essential-virtual-k7-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("report_sha256")
        != EXPECTED_FILES["d6_k6_empty_essential_virtual_k7.json"]
        or verification.get("input_graphs") != 805
        or verification.get("rejected_graphs") != 49
        or verification.get("ordered_residue_graphs") != 756
    ):
        raise ValueError("virtual-K7 independent verification is incomplete")

    rejected = report.get("aggregate", {}).get("rejected_indices")
    if (
        not isinstance(rejected, list)
        or len(rejected) != 49
        or len(set(rejected)) != 49
        or any(type(index) is not int for index in rejected)
        or rejected != verification.get("rejected_indices")
    ):
        raise ValueError("virtual-K7 rejection list is malformed or disagrees")
    rejection_set = set(rejected)

    k7_graphs = list(v4["classes"]["K7"]["graphs"])
    k6_before = list(v4["classes"]["K6_only"]["graphs"])
    k6_before_indices = [record["index"] for record in k6_before]
    if (
        len(k7_graphs) != 155
        or len(k6_before) != 805
        or stable_hash(k6_before_indices) != report.get("input_indices_sha256")
    ):
        raise ValueError("virtual-K7 input is not the ordered v4 K6 class")
    k6_graphs = [
        record for record in k6_before if record["index"] not in rejection_set
    ]
    k6_indices = [record["index"] for record in k6_graphs]
    if (
        len(k6_graphs) != 756
        or k6_indices != report["aggregate"]["ordered_residue_indices"]
        or stable_hash(k6_indices)
        != report["aggregate"]["ordered_residue_indices_sha256"]
    ):
        raise ValueError("post-virtual-K7 K6 residue ordering differs")

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
        "schema": "d6-current-exact-residue-v5",
        "status": "COMPLETE_EXACT_FILTER_UNION",
        "description": (
            "The v4 exact residue after the independently verified K6 "
            "empty-essential / virtual-K7 conjunction. Survival is not "
            "realizability."
        ),
        "class_order": ["K7", "K6_only"],
        "classes": {
            "K7": {
                "count": 155,
                "indices": k7_indices,
                "indices_sha256": stable_hash(k7_indices),
                "graphs": k7_graphs,
                "graphs_sha256": stable_hash(k7_graphs),
                "change_from_v4": 0,
            },
            "K6_only": {
                "v4_count": 805,
                "virtual_K7_rejection_count": 49,
                "virtual_K7_rejected_indices": rejected,
                "virtual_K7_rejected_indices_sha256": stable_hash(rejected),
                "count": 756,
                "indices": k6_indices,
                "indices_sha256": stable_hash(k6_indices),
                "graphs": k6_graphs,
                "graphs_sha256": stable_hash(k6_graphs),
            },
        },
        "combined": {
            "count": 911,
            "class_counts": {"K7": 155, "K6_only": 756},
            "ordered_indices_sha256": stable_hash(combined_indices),
            "sorted_indices_sha256": stable_hash(sorted(combined_indices)),
            "ordered_class_index_records_sha256": stable_hash(tagged_indices),
            "ordered_class_graph_records_sha256": stable_hash(tagged_graphs),
            "cross_class_overlap": len(set(k7_indices) & set(k6_indices)),
        },
        "positive_18_control": v4["positive_18_control"],
        "source_boundary": {
            "files": EXPECTED_FILES,
            "v4_ordered_indices_sha256": v4["combined"][
                "ordered_indices_sha256"
            ],
            "virtual_K7_input_indices_sha256": report[
                "input_indices_sha256"
            ],
            "virtual_K7_residue_indices_sha256": report["aggregate"][
                "ordered_residue_indices_sha256"
            ],
        },
        "semantics": {
            "edges": "required unit distances",
            "nonedges": "unconstrained and may also have distance one",
            "points": "distinct",
            "rejection": (
                "one required K6 seed has no all-nonempty Hall branch, and "
                "every allowed one- or two-empty state contains an apex whose "
                "branch-forced K7 exhausts every exact eligible cover"
            ),
            "residue": "exact filter non-rejection only",
        },
        "nonclaims": [
            "No v5 residue graph is claimed realizable.",
            "The 911-graph boundary is not yet a proof that f(6)=18.",
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
        "--output", type=Path, default=ROOT / "d6_current_residue_manifest_v5.json"
    )
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output.resolve(), manifest)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
                "K7": 155,
                "K6_only": 756,
                "combined": 911,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
