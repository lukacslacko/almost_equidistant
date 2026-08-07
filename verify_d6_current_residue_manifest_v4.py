#!/usr/bin/env python3
"""Independent structural checker for the exact d=6 v4 residue manifest."""

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
BUILDER = "build_d6_current_residue_manifest_v4.py"
EXPECTED_BUILDER_SHA256 = (
    "e3c085fcf3e1e8a77e8ac43d42e368d13e103458006bff1ed28a2a9678270178"
)
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path.name} is not an object")
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
                if not rows[first] & (1 << third) and not rows[second] & (1 << third):
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
            following = [other for other in candidates if rows[vertex] & (1 << other)]
            if visit(following, chosen + 1):
                return True
        return False

    return visit(list(range(19)), 0)


def verify(manifest_path: Path, expected_manifest_sha256: str | None) -> dict:
    checks: dict[str, bool] = {}
    manifest_hash = sha256(manifest_path)
    if expected_manifest_sha256 is not None:
        require(manifest_hash == expected_manifest_sha256, "v4 manifest hash")
    require(sha256(ROOT / BUILDER) == EXPECTED_BUILDER_SHA256, "builder hash")
    for name, expected in EXPECTED_FILES.items():
        require(sha256(ROOT / name) == expected, f"pinned file changed: {name}")
    checks["all_hash_pins"] = True

    manifest = load(manifest_path)
    v3 = load(ROOT / "d6_current_residue_manifest_v3.json")
    v3_check = load(ROOT / "d6_current_residue_manifest_v3_verification.json")
    report = load(ROOT / "d6_k6_empty_support_budget_report.json")
    archive = load(ROOT / "d6_k6_empty_support_budget_certificates.json")
    production_check = load(ROOT / "d6_k6_empty_support_budget_verification.json")
    require(manifest.get("schema") == "d6-current-exact-residue-v4", "schema")
    require(v3.get("schema") == "d6-current-exact-residue-v3", "v3 schema")
    require(v3_check.get("status") == "PASS", "v3 verification")
    require(report.get("status") == "COMPLETE", "production status")
    require(archive.get("status") == "COMPLETE", "archive status")
    require(production_check.get("status") == "PASS", "production verification")
    require(
        production_check.get("report_sha256")
        == EXPECTED_FILES["d6_k6_empty_support_budget_report.json"],
        "production report binding",
    )
    checks["upstream_statuses"] = True

    rejected = report.get("rejected_indices")
    require(
        isinstance(rejected, list)
        and len(rejected) == 17
        and len(set(rejected)) == 17
        and all(type(index) is int for index in rejected),
        "rejected index list",
    )
    require(rejected == archive.get("rejected_indices"), "archive rejection set")
    require(
        rejected == production_check.get("rejected_indices"),
        "verification rejection set",
    )
    rejected_set = set(rejected)

    v3_k7 = list(v3["classes"]["K7"]["graphs"])
    v3_k6 = list(v3["classes"]["K6_only"]["graphs"])
    require(
        [record["index"] for record in v3_k6] == report.get("ordered_input_indices"),
        "production input order",
    )
    expected_k6 = [
        record for record in v3_k6 if record["index"] not in rejected_set
    ]
    require(
        [record["index"] for record in expected_k6]
        == report.get("ordered_residue_indices"),
        "production residue order",
    )
    embedded_k7 = manifest["classes"]["K7"]["graphs"]
    embedded_k6 = manifest["classes"]["K6_only"]["graphs"]
    require(embedded_k7 == v3_k7, "K7 class changed")
    require(embedded_k6 == expected_k6, "K6 class filter mismatch")
    checks["exact_set_filter"] = True

    for class_name, records in (("K7", embedded_k7), ("K6_only", embedded_k6)):
        for record in records:
            rows = validate_rows(record.get("adjacency"), f"{class_name} {record.get('index')}")
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
    require(len(k7_indices) == 155 and len(k6_indices) == 805, "class counts")
    require(manifest["combined"]["count"] == 960, "combined count")
    require(
        manifest["combined"]["ordered_indices_sha256"] == stable_hash(combined),
        "combined order hash",
    )
    require(not set(k7_indices) & set(k6_indices), "cross-class overlap")
    require(
        manifest["classes"]["K6_only"]["empty_support_rejected_indices"]
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
        manifest.get("positive_18_control") == v3.get("positive_18_control"),
        "positive control changed",
    )
    semantics = manifest.get("semantics")
    require(isinstance(semantics, dict), "missing semantics")
    require("unconstrained" in semantics.get("nonedges", ""), "nonedge semantics")
    checks["semantics_and_positive_control"] = True

    # Audit our own import boundary; this checker must not inherit a production
    # decision through a Python import.
    syntax = ast.parse((ROOT / __file__).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require("build_d6_current_residue_manifest_v4" not in imports, "builder import")
    require(not any(name.startswith("d6_k6_") for name in imports), "K6 import")
    checks["import_independence"] = True

    return {
        "schema": "d6-current-exact-residue-v4-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "builder_sha256": EXPECTED_BUILDER_SHA256,
        "counts": {"K7": 155, "K6_only": 805, "combined": 960, "rejected": 17},
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
        "--manifest", type=Path, default=ROOT / "d6_current_residue_manifest_v4.json"
    )
    parser.add_argument("--expected-manifest-sha256")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v4_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.manifest.resolve(), args.expected_manifest_sha256)
    atomic_json(args.output, result)
    print(json.dumps({"status": "PASS", **result["counts"], "output": str(args.output)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
