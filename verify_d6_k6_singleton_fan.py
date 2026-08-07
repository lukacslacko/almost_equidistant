#!/usr/bin/env python3
"""Independent verifier for the K6 singleton--two-arm layer.

The checker imports neither the new producer nor its discovery probe.  It
extends the frozen independent tight-Hall checker with a fresh transcription
of the support-fan predicate and replays every one of the 649 graph decisions.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import verify_d6_k6_tight_hall_support as parent


ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-singleton-fan-v1"
CERTIFICATE_SCHEMA = "d6-k6-singleton-fan-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-singleton-fan-verification-v1"
EXPECTED_INPUT = 649
EXPECTED_INPUT_SHA256 = "6460baec9c378c1da3ba79f9ba6645e688158f2d9910bc6b88f4578c6e28e8d4"
PARENT_REPORT_SHA256 = "73b7fdd5b641d52a7eb3f2daaeea9d8e94057913e38f2bde06ef62691d01e864"
PARENT_VERIFICATION_SHA256 = "5fa924e47f50b23a0fb6289d1beba7a6f445ea599d2078a8a8925fd23802523c"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
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
        "d6_k6_singleton_fan",
        "probe_d6_k6_tight_support_singleton_fan",
    }
    if imports & forbidden:
        raise AssertionError(f"independent checker imports new code: {imports & forbidden}")


PARENT_SUPPORT_OBSTRUCTED = parent.support_obstructed


def singleton_fan_obstructed(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> bool:
    """Independent index-pair transcription of the forbidden support fan."""

    reduced = tuple(defects[vertex] & ~removed for vertex in absolute)
    for anchor in range(len(absolute)):
        singleton = reduced[anchor]
        if singleton.bit_count() != 1:
            continue
        arms = []
        for candidate in range(len(absolute)):
            mask = reduced[candidate]
            if (
                candidate != anchor
                and mask.bit_count() == 2
                and mask & singleton
                and graph[anchor] & (1 << candidate)
            ):
                arms.append(candidate)
        for left_position, left in enumerate(arms):
            for right in arms[:left_position]:
                if (reduced[left] ^ singleton) == (reduced[right] ^ singleton):
                    continue
                if graph[left] & (1 << right):
                    return True
    return False


def strengthened_support_obstructed(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> bool:
    return PARENT_SUPPORT_OBSTRUCTED(
        graph, absolute, defects, removed
    ) or singleton_fan_obstructed(graph, absolute, defects, removed)


def activate_kernel() -> None:
    parent.support_obstructed = strengthened_support_obstructed


def load_input() -> tuple[list[dict], list[int]]:
    parent_report_path = ROOT / "d6_k6_tight_hall_support_report.json"
    parent_verification_path = ROOT / "d6_k6_tight_hall_support_verification.json"
    if (
        sha256(parent_report_path) != PARENT_REPORT_SHA256
        or sha256(parent_verification_path) != PARENT_VERIFICATION_SHA256
    ):
        raise ValueError("frozen parent artifacts changed")
    report = json.loads(parent_report_path.read_text(encoding="utf-8"))
    verification = json.loads(parent_verification_path.read_text(encoding="utf-8"))
    if (
        report.get("graphs_rejected") != 107
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 107
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen parent contents changed")
    all_records, _ = parent.load_input()
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("independent 649 residue reconstruction failed")
    return [by_index[index] for index in indices], indices


def evaluate_record(record: dict) -> dict:
    activate_kernel()
    return parent.evaluate_record(record)


def positive_control() -> dict:
    activate_kernel()
    return parent.positive_control()


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
    source_path = ROOT / "d6_k6_singleton_fan.py"
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
    ):
        raise ValueError("singleton-fan production boundary mismatch")
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
        raise AssertionError("independent singleton-fan decisions differ")
    report_rows = {int(row["index"]): row for row in report["graph_results"]}
    archived = {int(row["index"]): row for row in archive["rejected_graphs"]}
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
        raise AssertionError("independent singleton-fan positive control failed")
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
            "independent_fan_index_transcription": True,
            "all_graph_decisions": True,
            "certificate_Z0_coverage": True,
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
        default=ROOT / "d6_k6_singleton_fan_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_singleton_fan_certificates.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_singleton_fan_verification.json",
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
