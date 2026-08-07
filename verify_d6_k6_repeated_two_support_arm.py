#!/usr/bin/env python3
"""Independent verifier for the K6 repeated-two-support-arm layer.

The checker imports neither the new producer nor its discovery probe.  It
extends the frozen independent singleton-fan checker with a separately
transcribed triple-loop predicate and replays every one of the 634 decisions.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import verify_d6_k6_singleton_fan as fan_verifier


ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-repeated-two-support-arm-v1"
CERTIFICATE_SCHEMA = "d6-k6-repeated-two-support-arm-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-repeated-two-support-arm-verification-v1"
EXPECTED_INPUT = 634
EXPECTED_INPUT_SHA256 = (
    "bafeef8b2f16419bbbfd1dbd8756081be78d14d0dc4187a1d954738a7edb9d88"
)
EXPECTED_PARENT_RAW_CERTIFICATE_SHA256 = (
    "b2fbf39d443679fee52d4984f3a54e4302acf416760087ed75e0953eb384d7bf"
)
EXPECTED_PARENT_ARTIFACTS = {
    "d6_k6_singleton_fan.py": (
        "b1755834995d081e633c5cf95e09333651bf47378c372b06c3377982605aee1f"
    ),
    "d6_k6_singleton_fan_report.json": (
        "bb11b82e8783108e45ec97a8e92957037bf219dd475b71fc18faee50318a4e2d"
    ),
    "d6_k6_singleton_fan_certificates.json.gz": (
        "fcce6c4e7758e4ffef107342a81d9e541d696d378592aad92701611280e5f467"
    ),
    "d6_k6_singleton_fan_verification.json": (
        "77011cf280b250aff59d43008af10525c936f853172040905c0cfc5324491820"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gunzipped_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    activate_kernel()


def assert_import_independence() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "d6_k6_repeated_two_support_arm",
        "probe_d6_k6_repeated_two_support_arm",
    }
    if imports & forbidden:
        raise AssertionError(f"independent checker imports new code: {imports & forbidden}")


FAN_SUPPORT_OBSTRUCTED = fan_verifier.strengthened_support_obstructed


def repeated_two_support_arm_obstructed(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> bool:
    """Independent ordered-triple transcription of the collision motif."""

    reduced = tuple(defects[vertex] & ~removed for vertex in absolute)
    for anchor in range(len(absolute)):
        singleton = reduced[anchor]
        if singleton.bit_count() != 1:
            continue
        for first in range(len(absolute)):
            first_mask = reduced[first]
            if (
                first == anchor
                or first_mask.bit_count() != 2
                or not (first_mask & singleton)
                or not (graph[anchor] & (1 << first))
            ):
                continue
            for second in range(first):
                if (
                    second != anchor
                    and reduced[second] == first_mask
                    and bool(graph[anchor] & (1 << second))
                ):
                    return True
    return False


def strengthened_support_obstructed(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> bool:
    return FAN_SUPPORT_OBSTRUCTED(
        graph, absolute, defects, removed
    ) or repeated_two_support_arm_obstructed(
        graph, absolute, defects, removed
    )


def activate_kernel() -> None:
    # fan_verifier.parent is the structurally independent tight-Hall checker.
    fan_verifier.parent.support_obstructed = strengthened_support_obstructed


def load_input() -> tuple[list[dict], list[int]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_PARENT_ARTIFACTS
    }
    if observed != EXPECTED_PARENT_ARTIFACTS:
        raise ValueError(f"frozen singleton-fan artifacts changed: {observed}")
    compressed = ROOT / "d6_k6_singleton_fan_certificates.json.gz"
    if gunzipped_sha256(compressed) != EXPECTED_PARENT_RAW_CERTIFICATE_SHA256:
        raise ValueError("singleton-fan compressed certificate payload changed")
    report_path = ROOT / "d6_k6_singleton_fan_report.json"
    verification_path = ROOT / "d6_k6_singleton_fan_verification.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    if (
        report.get("schema") != "d6-k6-singleton-fan-v1"
        or report.get("status") != "COMPLETE"
        or report.get("graphs_rejected") != 15
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
        or verification.get("schema")
        != "d6-k6-singleton-fan-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 15
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen singleton-fan contents changed")
    all_records, _ = fan_verifier.load_input()
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("independent 634 residue reconstruction failed")
    if any(index not in by_index for index in indices):
        raise ValueError("independent residue is outside the 649 parent input")
    return [by_index[index] for index in indices], indices


def evaluate_record(record: dict) -> dict:
    activate_kernel()
    return fan_verifier.parent.evaluate_record(record)


def positive_control() -> dict:
    activate_kernel()
    return fan_verifier.parent.positive_control()


def verify(
    report_path: Path,
    certificates_path: Path,
    output: Path,
    workers: int,
) -> dict:
    assert_import_independence()
    records, indices = load_input()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    archive = json.loads(certificates_path.read_text(encoding="utf-8"))
    source_path = ROOT / "d6_k6_repeated_two_support_arm.py"
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("production_source_sha256") != sha256(source_path)
        or report.get("ordered_input_indices") != indices
        or report.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(certificates_path)
        or archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("production_source_sha256") != sha256(source_path)
        or archive.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("repeated-arm production boundary mismatch")
    started = time.monotonic()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [row["index"] for row in results if row["rejected"]]
    residue = [row["index"] for row in results if not row["rejected"]]
    if (
        rejected != report.get("rejected_indices")
        or stable_hash(rejected) != report.get("rejected_indices_sha256")
        or residue != report.get("ordered_residue_indices")
        or stable_hash(residue) != report.get("ordered_residue_indices_sha256")
        or len(rejected) != report.get("graphs_rejected")
        or len(residue) != report.get("graphs_surviving")
    ):
        raise AssertionError("independent repeated-arm decisions differ")
    report_rows = {int(row["index"]): row for row in report["graph_results"]}
    archived = {int(row["index"]): row for row in archive["rejected_graphs"]}
    if set(archived) != set(rejected):
        raise AssertionError("certificate archive index coverage differs")
    for row in results:
        saved = report_rows[row["index"]]
        if (
            row["rejected"] != saved["rejected"]
            or row["seeds_checked"] != saved["seeds_checked"]
            or row["seed_mask"] != saved["first_impossible_seed_mask"]
            or row["seed"] != saved["first_impossible_seed"]
        ):
            raise AssertionError(f"graph detail differs at {row['index']}")
        if row["rejected"]:
            certificate = archived[row["index"]]
            if (
                certificate["seed_mask"] != row["seed_mask"]
                or certificate["seed"] != row["seed"]
                or len(certificate["rows"]) != len(row["rows"])
                or [item["Z0"] for item in certificate["rows"]]
                != [item["Z0"] for item in row["rows"]]
            ):
                raise AssertionError(f"archive coverage differs at {row['index']}")
    control = positive_control()
    if not control["passed"] or control["K6_seeds"] != 32:
        raise AssertionError("independent repeated-arm positive control failed")
    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": sha256(report_path)},
        "certificates": {
            "path": certificates_path.name,
            "sha256": sha256(certificates_path),
        },
        "checks": {
            "import_independence": True,
            "independent_ordered_triple_transcription": True,
            "independent_parent_Hall_and_inertia_kernel": True,
            "all_graph_decisions": True,
            "certificate_index_and_Z0_coverage": True,
            "positive_18_control": True,
        },
        "graphs_recomputed": len(results),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices_sha256": stable_hash(residue),
        "positive_18_control": {
            "passed": True,
            "K6_seeds": control["K6_seeds"],
        },
        "runtime": {
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k6_repeated_two_support_arm_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_repeated_two_support_arm_certificates.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_repeated_two_support_arm_verification.json",
    )
    parser.add_argument("--workers", type=int, default=11)
    args = parser.parse_args()
    result = verify(args.report, args.certificates, args.output, args.workers)
    print(json.dumps({
        "status": result["status"],
        "rejected": result["graphs_rejected"],
        "surviving": result["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
