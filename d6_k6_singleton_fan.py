#!/usr/bin/env python3
"""Exact singleton--two-arm algebra on the frozen 649-graph K6 residue.

This incremental layer preserves the frozen 107-hit tight-Hall result.  After
one of its exact tight coordinate deletions, it rejects a generic Lorentz-side
branch containing required points with reduced supports ``{i}``, ``{i,j}``,
``{i,k}`` for distinct ``j,k``.  The shared Lorentz ratio would have squared
value -5, impossible over the reals.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import d6_k6_tight_hall_support as parent


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_tight_hall_support_report.json"
PARENT_CERTIFICATES = ROOT / "d6_k6_tight_hall_support_certificates.json"
PARENT_VERIFICATION = ROOT / "d6_k6_tight_hall_support_verification.json"
REPORT_SCHEMA = "d6-k6-singleton-fan-v1"
CERTIFICATE_SCHEMA = "d6-k6-singleton-fan-certificates-v1"
EXPECTED_INPUT = 649
EXPECTED_INPUT_SHA256 = "6460baec9c378c1da3ba79f9ba6645e688158f2d9910bc6b88f4578c6e28e8d4"
DEFAULT_WORKERS = 11

EXPECTED_DEPENDENCIES = {
    "d6_k6_tight_hall_support.py": "6932918dac59c042bf22b505b6af84f923d4065ea9e3b94dd47d02e246e4fa1e",
    "d6_k6_tight_hall_support_report.json": "73b7fdd5b641d52a7eb3f2daaeea9d8e94057913e38f2bde06ef62691d01e864",
    "d6_k6_tight_hall_support_certificates.json": "a3b29e28db706de5b1d63ba4f82e2b527a53d6b0ac479d098c3bb2cccb66198f",
    "d6_k6_tight_hall_support_verification.json": "5fa924e47f50b23a0fb6289d1beba7a6f445ea599d2078a8a8925fd23802523c",
}


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


def singleton_fan_obstruction(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> dict | None:
    """Find a required ``{i}``, ``{i,j}``, ``{i,k}`` fan after deletion."""

    reduced = [defects[vertex] & ~removed for vertex in absolute]
    for anchor, anchor_mask in enumerate(reduced):
        if anchor_mask.bit_count() != 1:
            continue
        coordinate = parent.vertices(anchor_mask)[0]
        arms = []
        for arm, arm_mask in enumerate(reduced):
            if arm == anchor or arm_mask.bit_count() != 2:
                continue
            if not (arm_mask & anchor_mask):
                continue
            if not (graph[anchor] & (1 << arm)):
                continue
            other = parent.vertices(arm_mask ^ anchor_mask)[0]
            arms.append((arm, other))
        for position, (first, first_other) in enumerate(arms):
            for second, second_other in arms[:position]:
                if first_other == second_other:
                    continue
                if not (graph[first] & (1 << second)):
                    continue
                return {
                    "kind": "singleton_two_distinct_two_support_arms",
                    "vertices": [
                        absolute[anchor], absolute[second], absolute[first]
                    ],
                    "supports": [
                        parent.vertices(reduced[anchor]),
                        parent.vertices(reduced[second]),
                        parent.vertices(reduced[first]),
                    ],
                    "anchor_coordinate": coordinate,
                    "exact_contradiction": "(r_squared+5)^2=0",
                }
    return None


PARENT_SUPPORT_OBSTRUCTION = parent.support_obstruction


def strengthened_support_obstruction(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> dict | None:
    old = PARENT_SUPPORT_OBSTRUCTION(graph, absolute, defects, removed)
    if old is not None:
        return old
    return singleton_fan_obstruction(graph, absolute, defects, removed)


def activate_kernel() -> None:
    parent.support_obstruction = strengthened_support_obstruction


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"singleton-fan dependency boundary changed: {observed}")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    if (
        report.get("schema") != "d6-k6-tight-hall-support-v1"
        or report.get("status") != "COMPLETE"
        or report.get("graphs_rejected") != 107
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
        or verification.get("schema")
        != "d6-k6-tight-hall-support-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 107
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen tight-Hall result boundary changed")
    all_records, _, _ = parent.load_input()
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("649-graph residue changed")
    return [by_index[index] for index in indices], indices, observed


def evaluate_record(record: dict) -> dict:
    activate_kernel()
    return parent.evaluate_record(record)


def positive_control() -> dict:
    activate_kernel()
    result = parent.positive_control()
    if not result["passed"] or result["K6_seeds"] != 32:
        raise AssertionError("singleton-fan layer rejects positive 18 control")
    return result


def run(workers: int, output: Path, certificates: Path) -> dict:
    records, indices, dependencies = load_input()
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
    counts: Counter[str] = Counter()
    for row in results:
        counts.update(row["counts"])
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "production_source_sha256": sha256(Path(__file__)),
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "rejected_graphs": [
            {"index": row["index"], **row["certificate"]}
            for row in results if row["rejected"]
        ],
    }
    atomic_json(certificates, archive)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": "exact generic singleton--two-arm Lorentz contradiction",
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_coordinates_may_be_zero",
            "distinct_points": "required",
            "arithmetic": "exact_integer_bitmask_and_polynomial_identity",
        },
        "production_source_sha256": sha256(Path(__file__)),
        "dependencies": dependencies,
        "ordered_input_indices": indices,
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "input_graphs": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices": residue,
        "ordered_residue_indices_sha256": stable_hash(residue),
        "totals": {
            "seeds_checked": sum(row["seeds_checked"] for row in results),
            **dict(counts),
        },
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "rejected_graphs": len(archive["rejected_graphs"]),
        },
        "positive_18_control": positive_control(),
        "runtime": {
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "graph_results": [
            {key: value for key, value in row.items() if key != "certificate"}
            for row in results
        ],
        "nonclaims": [
            "a survivor is not a realization",
            "this incremental layer does not settle the K6 or dimension-six problem",
        ],
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_singleton_fan_report.json"
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_singleton_fan_certificates.json",
    )
    args = parser.parse_args()
    report = run(args.workers, args.output, args.certificates)
    print(json.dumps({
        "status": report["status"],
        "input": report["input_graphs"],
        "rejected": report["graphs_rejected"],
        "surviving": report["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
