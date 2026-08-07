#!/usr/bin/env python3
"""Exact repeated-two-support-arm collision on the frozen 634 K6 residue.

This incremental layer preserves the frozen tight-Hall and singleton-fan
decisions.  After one of their exact tight coordinate deletions, it rejects a
generic Lorentz-side branch containing three required points with reduced
allowed supports ``{i}``, ``{i,j}``, ``{i,j}``, with both arms required
adjacent to the anchor.  The singleton anchor and the K6 identities uniquely
reconstruct either arm, so two distinct arms collide.  Their mutual distance
is unused.
"""

from __future__ import annotations

import argparse
import gzip
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

import d6_k6_singleton_fan as fan_parent
import d6_k6_tight_hall_support as kernel


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_singleton_fan_report.json"
PARENT_CERTIFICATES_GZ = ROOT / "d6_k6_singleton_fan_certificates.json.gz"
PARENT_VERIFICATION = ROOT / "d6_k6_singleton_fan_verification.json"
REPORT_SCHEMA = "d6-k6-repeated-two-support-arm-v1"
CERTIFICATE_SCHEMA = "d6-k6-repeated-two-support-arm-certificates-v1"
EXPECTED_INPUT = 634
EXPECTED_INPUT_SHA256 = (
    "bafeef8b2f16419bbbfd1dbd8756081be78d14d0dc4187a1d954738a7edb9d88"
)
EXPECTED_PARENT_RAW_CERTIFICATE_SHA256 = (
    "b2fbf39d443679fee52d4984f3a54e4302acf416760087ed75e0953eb384d7bf"
)
DEFAULT_WORKERS = 11

EXPECTED_DEPENDENCIES = {
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


def repeated_two_support_arm_obstruction(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> dict | None:
    """Find an anchored ``{i}``, ``{i,j}``, ``{i,j}`` collision motif.

    The masks are only upper supports.  The proof does not assume that the
    second coordinate is nonzero: the diagonal identity forces its normalized
    value to ``(R+5)/12``, which is nonzero for real ``R=r^2``.
    """

    reduced = [defects[vertex] & ~removed for vertex in absolute]
    for anchor, singleton in enumerate(reduced):
        if singleton.bit_count() != 1:
            continue
        repeated_by_mask: dict[int, list[int]] = {}
        for arm, mask in enumerate(reduced):
            if (
                arm != anchor
                and mask.bit_count() == 2
                and bool(mask & singleton)
                and bool(graph[anchor] & (1 << arm))
            ):
                repeated_by_mask.setdefault(mask, []).append(arm)
        for mask, arms in repeated_by_mask.items():
            if len(arms) < 2:
                continue
            first, second = arms[:2]
            coordinate = kernel.vertices(singleton)[0]
            other = kernel.vertices(mask ^ singleton)[0]
            return {
                "kind": "singleton_two_identical_two_support_arms",
                "vertices": [
                    absolute[anchor], absolute[first], absolute[second]
                ],
                "supports": [
                    kernel.vertices(singleton),
                    kernel.vertices(mask),
                    kernel.vertices(mask),
                ],
                "anchor_coordinate": coordinate,
                "arm_coordinate": other,
                "required_anchor_pairs": [
                    [absolute[anchor], absolute[first]],
                    [absolute[anchor], absolute[second]],
                ],
                "arm_pair_distance": "unused",
                "optional_zero_resolution": (
                    "arm_second_normalized_coordinate=(r_squared+5)/12"
                ),
                "exact_conclusion": (
                    "unique_normalized_defect_and_finite_scale_force_collision"
                ),
            }
    return None


FAN_SUPPORT_OBSTRUCTION = fan_parent.strengthened_support_obstruction


def strengthened_support_obstruction(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> dict | None:
    old = FAN_SUPPORT_OBSTRUCTION(graph, absolute, defects, removed)
    if old is not None:
        return old
    return repeated_two_support_arm_obstruction(
        graph, absolute, defects, removed
    )


def activate_kernel() -> None:
    kernel.support_obstruction = strengthened_support_obstruction


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"repeated-arm dependency boundary changed: {observed}")
    if gunzipped_sha256(PARENT_CERTIFICATES_GZ) != (
        EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
    ):
        raise ValueError("singleton-fan compressed certificate payload changed")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
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
        raise ValueError("frozen singleton-fan result boundary changed")
    all_records, _, _ = fan_parent.load_input()
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("634-graph residue changed")
    if any(index not in by_index for index in indices):
        raise ValueError("singleton-fan residue is outside its 649-graph input")
    return [by_index[index] for index in indices], indices, observed


def evaluate_record(record: dict) -> dict:
    activate_kernel()
    return kernel.evaluate_record(record)


def positive_control() -> dict:
    activate_kernel()
    result = kernel.positive_control()
    if not result["passed"] or result["K6_seeds"] != 32:
        raise AssertionError("repeated-arm layer rejects positive 18 control")
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
        "description": "exact repeated two-support arm collision",
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_coordinates_may_be_zero",
            "optional_zero_branch": (
                "diagonal_equation_forces_arm_second_coordinate_nonzero"
            ),
            "distinct_points": "required",
            "arithmetic": "exact_integer_bitmask_and_rational_identity",
        },
        "production_source_sha256": sha256(Path(__file__)),
        "dependencies": dependencies,
        "parent_raw_certificate_payload_sha256": (
            EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
        ),
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
        "--output", type=Path,
        default=ROOT / "d6_k6_repeated_two_support_arm_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_repeated_two_support_arm_certificates.json",
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
