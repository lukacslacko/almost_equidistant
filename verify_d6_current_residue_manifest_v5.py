#!/usr/bin/env python3
"""Independent structural checker for the exact d=6 v5 residue manifest."""

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
BUILDER = "build_d6_current_residue_manifest_v5.py"
EXPECTED_BUILDER_SHA256 = (
    "2062abce9a4b9d075b794c03ad81a7f6403b43356b56691623b9fc9288a4d0a3"
)
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path.name} is not a JSON object")
    return value


def validate_rows(value: object, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == 19
        and all(type(row) is int for row in value),
        f"{label}: bad adjacency",
    )
    rows = tuple(value)
    for first in range(19):
        require(0 <= rows[first] < (1 << 19), f"{label}: row range")
        require(not rows[first] & (1 << first), f"{label}: loop")
        for second in range(first):
            require(
                bool(rows[first] & (1 << second))
                == bool(rows[second] & (1 << first)),
                f"{label}: asymmetry",
            )
    return rows


def alpha_two(rows: Sequence[int]) -> bool:
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


def clique(rows: Sequence[int], target: int) -> bool:
    def visit(candidates: list[int], chosen: int) -> bool:
        if chosen == target:
            return True
        if chosen + len(candidates) < target:
            return False
        while candidates:
            vertex = candidates.pop()
            following = [
                other for other in candidates if rows[vertex] & (1 << other)
            ]
            if visit(following, chosen + 1):
                return True
        return False

    return visit(list(range(19)), 0)


def verify(manifest_path: Path, expected_manifest_sha256: str | None) -> dict:
    checks: dict[str, bool] = {}
    manifest_hash = sha256(manifest_path)
    if expected_manifest_sha256 is not None:
        require(manifest_hash == expected_manifest_sha256, "v5 manifest hash")
    require(sha256(ROOT / BUILDER) == EXPECTED_BUILDER_SHA256, "builder hash")
    for name, expected in EXPECTED_FILES.items():
        require(sha256(ROOT / name) == expected, f"pinned file changed: {name}")
    checks["all_hash_pins"] = True

    manifest = load(manifest_path)
    v4 = load(ROOT / "d6_current_residue_manifest_v4.json")
    v4_check = load(ROOT / "d6_current_residue_manifest_v4_verification.json")
    report = load(ROOT / "d6_k6_empty_essential_virtual_k7.json")
    production_check = load(
        ROOT / "d6_k6_empty_essential_virtual_k7_verification.json"
    )
    require(manifest.get("schema") == "d6-current-exact-residue-v5", "schema")
    require(v4.get("schema") == "d6-current-exact-residue-v4", "v4 schema")
    require(v4_check.get("status") == "PASS", "v4 verification")
    require(
        v4_check.get("manifest", {}).get("sha256")
        == EXPECTED_FILES["d6_current_residue_manifest_v4.json"],
        "v4 report binding",
    )
    require(report.get("status") == "COMPLETE", "production status")
    require(production_check.get("status") == "PASS", "production verification")
    require(
        production_check.get("report_sha256")
        == EXPECTED_FILES["d6_k6_empty_essential_virtual_k7.json"],
        "production report binding",
    )
    checks["upstream_statuses"] = True

    rejected = report.get("aggregate", {}).get("rejected_indices")
    require(
        isinstance(rejected, list)
        and len(rejected) == 49
        and len(set(rejected)) == 49
        and all(type(index) is int for index in rejected),
        "rejected index list",
    )
    require(
        rejected == production_check.get("rejected_indices"),
        "verification rejection list",
    )
    rejected_set = set(rejected)

    v4_k7 = list(v4["classes"]["K7"]["graphs"])
    v4_k6 = list(v4["classes"]["K6_only"]["graphs"])
    require(
        stable_hash([record["index"] for record in v4_k6])
        == report.get("input_indices_sha256"),
        "production input order",
    )
    expected_k6 = [
        record for record in v4_k6 if record["index"] not in rejected_set
    ]
    require(
        [record["index"] for record in expected_k6]
        == report["aggregate"]["ordered_residue_indices"],
        "production residue order",
    )
    embedded_k7 = manifest["classes"]["K7"]["graphs"]
    embedded_k6 = manifest["classes"]["K6_only"]["graphs"]
    require(embedded_k7 == v4_k7, "K7 class changed")
    require(embedded_k6 == expected_k6, "K6 class exact filter mismatch")
    checks["exact_set_filter"] = True

    for class_name, records in (("K7", embedded_k7), ("K6_only", embedded_k6)):
        for record in records:
            rows = validate_rows(
                record.get("adjacency"), f"{class_name} {record.get('index')}"
            )
            require(alpha_two(rows), "independent triple")
            require(clique(rows, 6), "missing K6")
            if class_name == "K7":
                require(clique(rows, 7), "K7 class missing K7")
            else:
                require(not clique(rows, 7), "K6-only class has K7")
    checks["graph_structure"] = True

    k7_indices = [record["index"] for record in embedded_k7]
    k6_indices = [record["index"] for record in embedded_k6]
    combined = k7_indices + k6_indices
    require(len(k7_indices) == 155 and len(k6_indices) == 756, "class counts")
    require(manifest["combined"]["count"] == 911, "combined count")
    require(
        manifest["combined"]["ordered_indices_sha256"] == stable_hash(combined),
        "combined order hash",
    )
    require(not set(k7_indices) & set(k6_indices), "cross-class overlap")
    require(
        manifest["classes"]["K6_only"]["virtual_K7_rejected_indices"]
        == rejected,
        "embedded rejection list",
    )
    require(
        manifest["classes"]["K6_only"]["graphs_sha256"]
        == stable_hash(embedded_k6),
        "K6 graph hash",
    )
    checks["counts_and_hashes"] = True

    require(
        manifest.get("positive_18_control") == v4.get("positive_18_control"),
        "positive control changed",
    )
    semantics = manifest.get("semantics")
    require(isinstance(semantics, dict), "missing semantics")
    require("unconstrained" in semantics.get("nonedges", ""), "nonedge semantics")
    checks["semantics_and_positive_control"] = True

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require("build_d6_current_residue_manifest_v5" not in imports, "builder import")
    require(not any(name.startswith("d6_k6_") for name in imports), "K6 import")
    checks["import_independence"] = True

    return {
        "schema": "d6-current-exact-residue-v5-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "builder_sha256": EXPECTED_BUILDER_SHA256,
        "counts": {
            "K7": 155,
            "K6_only": 756,
            "combined": 911,
            "rejected_from_v4": 49,
        },
        "rejected_indices": rejected,
        "ordered_residue_indices_sha256": stable_hash(combined),
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
        "--manifest", type=Path, default=ROOT / "d6_current_residue_manifest_v5.json"
    )
    parser.add_argument("--expected-manifest-sha256")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v5_verification.json",
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
