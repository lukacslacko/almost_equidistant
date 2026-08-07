#!/usr/bin/env python3
"""Independent complete replay of the full K7 joint-support campaign.

This verifier imports neither the production full runner, its pilot evaluator,
nor the support-capacity kernel.  It validates every archive binding and then
uses the previously independent support/sparse-value implementations to
exhaust every recorded failing K7 seed against the inherited exact cover
certificates.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import gzip
import hashlib
import io
import json
import os
import platform
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import d6_k7_rank_reference as reference
import verify_d6_k7_support_capacity_pilot as independent
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
RANK_INPUT = ROOT / ".runs/d6_k7_rank_survivors.json"
UNION = ROOT / "d6_k7_rankone_pattern_union.json"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
QOS_CLASS_USER_INITIATED = 0x19
EXPECTED_INDEPENDENT_SOURCE_SHA256 = (
    "d16e7e9037d14aea797ff90df51df2ca1bb583b412a3f65992d995964263da6b"
)
DECISION_FIELDS = (
    "ordinal", "index", "status", "seeds", "covers",
    "current_passing_covers", "joint_pre_capacity_failing_covers",
    "pre_capacity_passing_covers", "capacity_passing_covers",
    "capacity_incremental_failing_covers", "capacity_unresolved_covers",
    "labeled_z_families", "pre_capacity_passing_families",
    "capacity_feasible_families", "capacity_infeasible_families",
    "capacity_unresolved_families", "capacity_nodes", "error_type",
)
INTEGER_FIELDS = DECISION_FIELDS[0:2] + DECISION_FIELDS[3:-1]


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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def resolve(report_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    local = report_path.parent / path
    return local if local.exists() else ROOT / path


def read_gzip(path: Path) -> bytes:
    with gzip.open(path, "rb") as stream:
        return stream.read()


def initialize_worker() -> None:
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    if sys.platform != "darwin":
        return
    libc = ctypes.CDLL(None)
    setter = libc.pthread_set_qos_class_self_np
    setter.argtypes = (ctypes.c_uint, ctypes.c_int)
    setter.restype = ctypes.c_int
    error = int(setter(QOS_CLASS_USER_INITIATED, 0))
    if error:
        raise OSError(error, "pthread_set_qos_class_self_np failed")


def replay_task(payload: tuple[dict, int, list, list]) -> dict:
    graph, seed_mask, prior_keys, tetrad_keys = payload
    started = time.perf_counter()
    try:
        result = independent.replay_failing_seed(
            graph,
            seed_mask,
            {(tuple(seed), int(zmask)) for seed, zmask in prior_keys},
            {(tuple(seed), int(zmask)) for seed, zmask in tetrad_keys},
        )
        return {
            "status": "PASS",
            "result": result,
            "elapsed_seconds": time.perf_counter() - started,
        }
    except BaseException as error:  # noqa: BLE001 - verification record
        return {
            "status": "ERROR",
            "result": None,
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_joint_support_full_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_joint_support_full_verification.json",
    )
    args = parser.parse_args()
    if args.workers <= 0:
        raise SystemExit("workers must be positive")
    report_hash = sha256(args.report)
    if report_hash != args.report_sha256:
        raise SystemExit(f"report hash {report_hash} != explicit pin")
    if sha256(ROOT / "verify_d6_k7_support_capacity_pilot.py") != (
        EXPECTED_INDEPENDENT_SOURCE_SHA256
    ):
        raise ValueError("independent bounded verifier source changed")
    report = json.loads(args.report.read_text(encoding="utf-8"))
    if report.get("status") != "COMPLETE" or report.get("schema") != 1:
        raise ValueError("full production report is not complete schema one")
    selection = report["selection"]
    union = json.loads(UNION.read_text(encoding="utf-8"))
    indices = [int(value) for value in union["sets"]["exact_residue"]["indices"]]
    if (
        len(indices) != 258
        or selection["graphs"] != 258
        or selection["indices_sha256"] != stable_hash(indices)
        or selection["source_sha256"] != sha256(UNION)
    ):
        raise ValueError("full selection does not reconstruct")
    artifacts = report["artifacts"]
    decisions_path = resolve(args.report, artifacts["decisions"])
    certificates_path = resolve(args.report, artifacts["certificates"])
    checkpoint_path = resolve(args.report, artifacts["checkpoint_copy"])
    for path, key in (
        (decisions_path, "decisions_sha256"),
        (certificates_path, "certificates_sha256"),
        (checkpoint_path, "checkpoint_copy_sha256"),
    ):
        if sha256(path) != artifacts[key]:
            raise ValueError(f"artifact hash mismatch: {path}")
    decisions_raw = read_gzip(decisions_path)
    certificates_raw = read_gzip(certificates_path)
    if hashlib.sha256(decisions_raw).hexdigest() != artifacts[
        "decisions_uncompressed_sha256"
    ] or hashlib.sha256(certificates_raw).hexdigest() != artifacts[
        "certificates_uncompressed_sha256"
    ]:
        raise ValueError("uncompressed archive hash mismatch")
    reader = csv.DictReader(
        io.StringIO(decisions_raw.decode("ascii"), newline=""), delimiter="\t"
    )
    if tuple(reader.fieldnames or ()) != DECISION_FIELDS:
        raise ValueError("decision schema mismatch")
    rows = list(reader)
    if len(rows) != 258:
        raise ValueError("decision population mismatch")
    statuses: Counter = Counter()
    column_totals: Counter = Counter()
    for ordinal, (row, index) in enumerate(zip(rows, indices)):
        if int(row["ordinal"]) != ordinal or int(row["index"]) != index:
            raise ValueError("decision order mismatch")
        if row["status"] not in {
            "SURVIVOR", "JOINT_PRE_CAPACITY_REJECTED",
            "CAPACITY_REJECTED", "UNRESOLVED",
        } or row["error_type"]:
            raise ValueError(f"bad decision status at graph {index}")
        statuses[row["status"]] += 1
        for name in INTEGER_FIELDS[2:]:
            value = int(row[name])
            if value < 0:
                raise ValueError("negative decision count")
            column_totals[name] += value
    if dict(sorted(statuses.items())) != report["summary"]["status_counts"]:
        raise ValueError("status summary mismatch")
    reported_rejected = [
        int(row["index"]) for row in rows
        if row["status"] == "JOINT_PRE_CAPACITY_REJECTED"
    ]
    if reported_rejected != report["summary"]["joint_rejected_indices"]:
        raise ValueError("joint rejection index list mismatch")
    if statuses.get("CAPACITY_REJECTED") or statuses.get("UNRESOLVED"):
        raise ValueError("full campaign contains a capacity-only or unresolved nonclaim")

    certificates = [
        json.loads(line) for line in certificates_raw.decode("ascii").splitlines()
        if line
    ]
    if (
        len(certificates) != artifacts["certificate_entries"]
        or [int(item["index"]) for item in certificates] != reported_rejected
    ):
        raise ValueError("certificate/rejection population mismatch")
    certificate_by_index = {}
    tasks = []
    selected_set = set(indices)
    prior = independent.load_witness_keys(
        PRIOR_CERTIFICATES, selected_set, "dual_failure_witnesses"
    )
    tetrad = independent.load_witness_keys(
        TETRAD_CERTIFICATES, selected_set, "tetrad_failure_witnesses"
    )
    rank_payload = json.loads(RANK_INPUT.read_text(encoding="utf-8"))
    graph_by_index = {int(g["index"]): g for g in rank_payload["graphs"]}
    for item in certificates:
        index = int(item["index"])
        if item.get("kind") != "joint_k7_cover_support_rejection":
            raise ValueError("bad certificate kind")
        if index in certificate_by_index or not item.get("failing_seeds"):
            raise ValueError("bad certificate identity/seed population")
        certificate_by_index[index] = item
        for seed in item["failing_seeds"]:
            if not seed.get("current_cover_failures"):
                raise ValueError("failing seed lacks archived cover failures")
            for cover in seed["current_cover_failures"]:
                if cover.get("first_pre_capacity_failure") is None:
                    raise ValueError("cover lacks first failure diagnostic")
            tasks.append((
                graph_by_index[index], int(seed["seed_mask"]),
                [[list(key[0]), key[1]] for key in sorted(prior.get(index, set()))],
                [[list(key[0]), key[1]] for key in sorted(tetrad.get(index, set()))],
            ))

    positive = tuple(lower_bound_18_graph())
    reference.validate_graph(positive)
    if tuple(reference.clique_masks(positive, 7)):
        raise AssertionError("18-point lower-bound control unexpectedly has K7")

    started = time.monotonic()
    with ProcessPoolExecutor(
        max_workers=args.workers, initializer=initialize_worker
    ) as executor:
        replayed = list(executor.map(replay_task, tasks, chunksize=1))
    errors = [item for item in replayed if item["status"] != "PASS"]
    if errors:
        raise RuntimeError(f"independent replay errors: {errors[:3]}")
    replay_by_key = {
        (int(item["result"]["index"]), int(item["result"]["seed_mask"])):
        item["result"]
        for item in replayed
    }
    if len(replay_by_key) != len(tasks):
        raise ValueError("duplicate independent replay key")
    for index, certificate in certificate_by_index.items():
        for seed in certificate["failing_seeds"]:
            key = (index, int(seed["seed_mask"]))
            replay = replay_by_key.get(key)
            if replay is None:
                raise ValueError("missing independent seed replay")
            if replay["counts"]["current_passing_covers"] != len(
                seed["current_cover_failures"]
            ):
                raise ValueError("independent current-cover count mismatch")

    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if (
        checkpoint.get("status") != "COMPLETE"
        or checkpoint.get("completed") != 258
        or checkpoint.get("status_counts") != dict(sorted(statuses.items()))
        or checkpoint.get("config_sha256") != report["configuration_sha256"]
        or checkpoint.get("config") != report["configuration"]
    ):
        raise ValueError("checkpoint-copy binding mismatch")
    output = {
        "schema": 1,
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": report_hash},
        "verifier_source_sha256": sha256(Path(__file__)),
        "independent_kernel_source_sha256": EXPECTED_INDEPENDENT_SOURCE_SHA256,
        "selection": {"graphs": 258, "indices_sha256": stable_hash(indices)},
        "status_counts": dict(sorted(statuses.items())),
        "joint_rejections_verified": len(reported_rejected),
        "failing_seeds_replayed": len(replayed),
        "cover_failures_reconstructed": sum(
            result["result"]["counts"]["current_passing_covers"]
            for result in replayed
        ),
        "capacity_only_hits": 0,
        "unresolved": 0,
        "positive_18_control": {
            "status": "NOT_APPLICABLE_NO_K7", "vertices": 18, "K7_seeds": 0,
        },
        "environment": {
            "workers": args.workers,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "elapsed_seconds": time.monotonic() - started,
        },
        "replay_totals": dict(sorted(sum(
            (Counter(item["result"]["counts"]) for item in replayed), Counter()
        ).items())),
    }
    atomic_json(args.output, output)
    print(
        f"PASS: {len(reported_rejected)} joint graph rejections, "
        f"{len(replayed)} failing seeds; SHA-256 {sha256(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
