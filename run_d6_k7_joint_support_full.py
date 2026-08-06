#!/usr/bin/env python3
"""Restartable full exact conjunction of K7 cover and support layers.

The target is the frozen 258-graph complement of the independently verified
rank-one-tetrad/pattern-954 union.  For every K7 seed, this runner intersects
the exact cover survivors of the strict-H, degree-one, and tetrad layers with
the exact labeled-support propagation and sparse-value survivors.  A graph is
``JOINT_PRE_CAPACITY_REJECTED`` only when one required K7 seed has no jointly
surviving cover/family.

The support-type capacity CSP remains measured as a separate adjunct.  It is
never credited to a joint rejection.  Node-cap exhaustion is recorded as
``UNRESOLVED`` and cannot become a mathematical rejection.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import fcntl
import gzip
import hashlib
import io
import json
import os
import platform
import sys
import time
import traceback
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import run_d6_k7_support_capacity_pilot as pilot


ROOT = Path(__file__).resolve().parent
RANK_INPUT = ROOT / ".runs/d6_k7_rank_survivors.json"
UNION = ROOT / "d6_k7_rankone_pattern_union.json"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"
EXPECTED_PILOT_SOURCE_SHA256 = (
    "9af542217d178bec2a71cb4c30faabc3dea8c8279053d729d924d9ca35eb128a"
)
QOS_CLASS_USER_INITIATED = 0x19

DECISION_FIELDS = (
    "ordinal", "index", "status", "seeds", "covers",
    "current_passing_covers", "joint_pre_capacity_failing_covers",
    "pre_capacity_passing_covers", "capacity_passing_covers",
    "capacity_incremental_failing_covers", "capacity_unresolved_covers",
    "labeled_z_families", "pre_capacity_passing_families",
    "capacity_feasible_families", "capacity_infeasible_families",
    "capacity_unresolved_families", "capacity_nodes", "error_type",
)


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


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


def fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def deterministic_gzip(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with temporary.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0
        ) as stream:
            stream.write(content)
        raw.flush()
        os.fsync(raw.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)


@contextmanager
def exclusive_lock(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "campaign.lock"
    with path.open("a+", encoding="ascii") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"checkpoint directory is locked: {directory}") from error
        stream.seek(0)
        stream.truncate()
        stream.write(f"pid={os.getpid()} started={utc_now()}\n")
        stream.flush()
        os.fsync(stream.fileno())
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def git_provenance() -> dict:
    import subprocess

    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        ).stdout.strip()

    try:
        status = git("status", "--porcelain=v1", "--untracked-files=all")
        return {
            "available": True,
            "commit": git("rev-parse", "HEAD"),
            "branch": git("branch", "--show-current"),
            "dirty": bool(status),
            "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
            "porcelain_lines": status.splitlines(),
        }
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}


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


def evaluate(payload: tuple[dict, list, list, dict, int]) -> dict:
    graph, prior_keys, tetrad_keys, row, node_limit = payload
    started = time.perf_counter()
    try:
        result = pilot.analyze_graph(
            graph,
            {((tuple(seed)), int(zmask)) for seed, zmask in prior_keys},
            {((tuple(seed)), int(zmask)) for seed, zmask in tetrad_keys},
            row,
            node_limit,
        )
        return {
            "schema": 1,
            "ordinal": int(graph["ordinal"]),
            "index": int(graph["index"]),
            "status": result["decision"],
            "result": result,
            "error_type": "",
            "error": "",
        }
    except BaseException as error:  # noqa: BLE001 - durable infra record
        return {
            "schema": 1,
            "ordinal": int(graph["ordinal"]),
            "index": int(graph["index"]),
            "status": "INFRA_ERROR",
            "result": None,
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
        }


def checkpoint_path(directory: Path, ordinal: int, index: int) -> Path:
    return directory / f"graph_{ordinal:03d}_{index}.json"


def load_checkpoints(
    directory: Path, config_sha256: str, selected: Sequence[dict]
) -> dict[int, dict]:
    expected = {int(g["ordinal"]): int(g["index"]) for g in selected}
    output = {}
    for path in sorted(directory.glob("graph_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != 1 or payload.get("config_sha256") != config_sha256:
            raise ValueError(f"incompatible checkpoint: {path}")
        record = payload.get("record")
        ordinal = int(record["ordinal"])
        index = int(record["index"])
        if expected.get(ordinal) != index or ordinal in output:
            raise ValueError(f"bad checkpoint identity: {path}")
        if record["status"] not in {
            "SURVIVOR", "JOINT_PRE_CAPACITY_REJECTED",
            "CAPACITY_REJECTED", "UNRESOLVED", "INFRA_ERROR",
        }:
            raise ValueError(f"bad checkpoint status: {path}")
        output[ordinal] = record
    return output


def decision_row(record: dict) -> tuple:
    result = record.get("result") or {}
    counts = result.get("counts", {})
    return (
        int(record["ordinal"]), int(record["index"]), record["status"],
        int(counts.get("seeds", 0)), int(counts.get("covers", 0)),
        int(counts.get("cover_passing", 0)),
        int(counts.get("joint_pre_capacity_failing_covers", 0)),
        int(counts.get("pre_capacity_passing_covers", 0)),
        int(counts.get("capacity_passing_covers", 0)),
        int(counts.get("capacity_incremental_failing_covers", 0)),
        int(counts.get("capacity_unresolved_covers", 0)),
        int(counts.get("labeled_z_families", 0)),
        int(counts.get("pre_capacity_passing_families", 0)),
        int(counts.get("capacity_feasible_families", 0)),
        int(counts.get("capacity_infeasible_families", 0)),
        int(counts.get("capacity_unresolved_families", 0)),
        int(counts.get("capacity_nodes", 0)), record.get("error_type", ""),
    )


def certificate_entry(record: dict) -> dict | None:
    if record["status"] != "JOINT_PRE_CAPACITY_REJECTED":
        return None
    result = record["result"]
    failing = set(int(value) for value in result["pre_capacity_rejected_seeds"])
    seeds = []
    for seed in result["seed_records"]:
        if int(seed["seed_mask"]) not in failing:
            continue
        seeds.append({
            "seed_mask": int(seed["seed_mask"]),
            "seed": list(seed["seed"]),
            "current_cover_failures": [
                {
                    "zmask": int(cover["zmask"]),
                    "counts": cover["counts"],
                    "first_pre_capacity_failure": cover[
                        "first_pre_capacity_failure"
                    ],
                }
                for cover in seed["current_cover_records"]
            ],
        })
    if not seeds:
        raise ValueError("joint rejection lacks a failing-seed certificate")
    return {
        "schema": 1,
        "ordinal": int(record["ordinal"]),
        "index": int(record["index"]),
        "kind": "joint_k7_cover_support_rejection",
        "failing_seeds": seeds,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--max-inflight", type=int, default=22)
    parser.add_argument("--node-limit", type=int, default=200_000)
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=ROOT / ".runs/d6_k7_joint_support_full",
    )
    parser.add_argument(
        "--decisions", type=Path,
        default=ROOT / "d6_k7_joint_support_full_decisions.tsv.gz",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k7_joint_support_full_certificates.jsonl.gz",
    )
    parser.add_argument(
        "--checkpoint-copy", type=Path,
        default=ROOT / "d6_k7_joint_support_full_checkpoint.json",
    )
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_joint_support_full_report.json",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--progress-every", type=int, default=16)
    args = parser.parse_args()
    if args.workers <= 0 or args.max_inflight < args.workers:
        raise SystemExit("workers must be positive and max-inflight >= workers")
    if args.node_limit < 0:
        raise SystemExit("node limit must be nonnegative")

    pilot_hash = sha256(ROOT / "run_d6_k7_support_capacity_pilot.py")
    if EXPECTED_PILOT_SOURCE_SHA256 == "TO_BE_PINNED":
        raise SystemExit("pin EXPECTED_PILOT_SOURCE_SHA256 before production")
    if pilot_hash != EXPECTED_PILOT_SOURCE_SHA256:
        raise SystemExit(f"pilot evaluator hash mismatch: {pilot_hash}")
    dependency_hashes = pilot.verify_inputs()
    union = json.loads(UNION.read_text(encoding="utf-8"))
    indices = [int(value) for value in union["sets"]["exact_residue"]["indices"]]
    if len(indices) != 258 or stable_hash(indices) != (
        "55ab329dbe3ae4dbfa10378ff023a168fc6ab306ffc7af14bbd9790a68108a09"
    ):
        raise ValueError("exact union residue mismatch")
    rank_payload = json.loads(RANK_INPUT.read_text(encoding="utf-8"))
    graph_by_index = {int(g["index"]): g for g in rank_payload["graphs"]}
    selected = [
        {**graph_by_index[index], "ordinal": ordinal}
        for ordinal, index in enumerate(indices)
    ]
    selected_set = set(indices)
    prior = pilot.load_jsonl_witness_keys(
        PRIOR_CERTIFICATES, selected_set, "dual_failure_witnesses"
    )
    tetrad = pilot.load_jsonl_witness_keys(
        TETRAD_CERTIFICATES, selected_set, "tetrad_failure_witnesses"
    )
    rows = pilot.load_tetrad_rows(selected_set)
    pilot.load_sparse_survivors(selected_set)

    source_hashes = {
        "run_d6_k7_joint_support_full.py": sha256(Path(__file__)),
        "run_d6_k7_support_capacity_pilot.py": pilot_hash,
        "d6_k7_support_capacity.py": sha256(ROOT / "d6_k7_support_capacity.py"),
        **dependency_hashes,
    }
    config = {
        "schema": 1,
        "selection_indices_sha256": stable_hash(indices),
        "selection_count": len(indices),
        "node_limit": args.node_limit,
        "source_sha256": source_hashes,
        "decision_fields": list(DECISION_FIELDS),
    }
    config_sha256 = stable_hash(config)
    command = [sys.executable, *sys.argv]
    started_unix = time.time()
    started = time.monotonic()
    with exclusive_lock(args.checkpoint_dir):
        existing = load_checkpoints(args.checkpoint_dir, config_sha256, selected)
        if existing and not args.resume:
            raise SystemExit(
                f"{len(existing)} checkpoints exist; pass --resume or use a new directory"
            )
        missing = deque(
            graph for graph in selected if int(graph["ordinal"]) not in existing
        )
        records = dict(existing)
        processed_session = 0
        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=initialize_worker
        ) as executor:
            active = {}
            while missing or active:
                while missing and len(active) < args.max_inflight:
                    graph = missing.popleft()
                    index = int(graph["index"])
                    payload = (
                        graph,
                        [[list(seed), zmask] for seed, zmask in sorted(prior.get(index, set()))],
                        [[list(seed), zmask] for seed, zmask in sorted(tetrad.get(index, set()))],
                        rows[index], args.node_limit,
                    )
                    future = executor.submit(evaluate, payload)
                    active[future] = graph
                done, _ = wait(tuple(active), return_when=FIRST_COMPLETED)
                for future in done:
                    graph = active.pop(future)
                    record = future.result()
                    ordinal = int(graph["ordinal"])
                    index = int(graph["index"])
                    if record["ordinal"] != ordinal or record["index"] != index:
                        raise ValueError("worker returned wrong graph identity")
                    atomic_json(
                        checkpoint_path(args.checkpoint_dir, ordinal, index),
                        {"schema": 1, "config_sha256": config_sha256, "record": record},
                    )
                    records[ordinal] = record
                    processed_session += 1
                    if processed_session % args.progress_every == 0:
                        statuses = Counter(item["status"] for item in records.values())
                        print(
                            f"progress {len(records)}/{len(selected)}; "
                            f"session {processed_session}; {dict(sorted(statuses.items()))}",
                            flush=True,
                        )

    ordered = [records[ordinal] for ordinal in range(len(selected))]
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(DECISION_FIELDS)
    writer.writerows(decision_row(record) for record in ordered)
    decisions_raw = buffer.getvalue().encode("ascii")
    deterministic_gzip(args.decisions, decisions_raw)
    certificate_lines = []
    for record in ordered:
        entry = certificate_entry(record)
        if entry is not None:
            certificate_lines.append(json.dumps(
                entry, sort_keys=True, separators=(",", ":"), ensure_ascii=True
            ))
    certificate_raw = (
        ("\n".join(certificate_lines) + "\n").encode("ascii")
        if certificate_lines else b""
    )
    deterministic_gzip(args.certificates, certificate_raw)
    statuses = Counter(record["status"] for record in ordered)
    totals: Counter = Counter()
    for record in ordered:
        if record["result"] is not None:
            totals.update(record["result"]["counts"])
    checkpoint_manifest = {
        "schema": 1,
        "status": "COMPLETE" if not statuses.get("INFRA_ERROR") else "INFRA_ERROR",
        "config": config,
        "config_sha256": config_sha256,
        "completed": len(ordered),
        "status_counts": dict(sorted(statuses.items())),
        "checkpoint_directory": str(args.checkpoint_dir),
        "checkpoint_files_sha256": stable_hash([
            [path.name, sha256(path)]
            for path in sorted(args.checkpoint_dir.glob("graph_*.json"))
        ]),
    }
    atomic_json(args.checkpoint_copy, checkpoint_manifest)
    report = {
        "schema": 1,
        "kind": "d6_k7_joint_support_full",
        "status": checkpoint_manifest["status"],
        "claim": (
            "JOINT_PRE_CAPACITY_REJECTED only: one K7 seed has no cover "
            "surviving both the inherited cover certificates and frozen "
            "labeled-support/sparse-value quantifiers. Capacity-only and "
            "unresolved statuses are explicit nonclaims."
        ),
        "selection": {
            "graphs": len(indices),
            "indices_sha256": stable_hash(indices),
            "source": UNION.name,
            "source_sha256": sha256(UNION),
        },
        "configuration": config,
        "configuration_sha256": config_sha256,
        "source_sha256": source_hashes,
        "command": command,
        "git": git_provenance(),
        "environment": {
            "workers": args.workers,
            "max_inflight": args.max_inflight,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
        "run": {
            "started_unix": started_unix,
            "finished_utc": utc_now(),
            "elapsed_seconds": time.monotonic() - started,
            "resumed_records": len(existing),
            "processed_records": processed_session,
        },
        "summary": {
            "status_counts": dict(sorted(statuses.items())),
            "joint_rejected_indices": [
                int(record["index"]) for record in ordered
                if record["status"] == "JOINT_PRE_CAPACITY_REJECTED"
            ],
            "totals": dict(sorted(totals.items())),
        },
        "artifacts": {
            "decisions": str(args.decisions),
            "decisions_sha256": sha256(args.decisions),
            "decisions_uncompressed_sha256": hashlib.sha256(decisions_raw).hexdigest(),
            "certificates": str(args.certificates),
            "certificates_sha256": sha256(args.certificates),
            "certificates_uncompressed_sha256": hashlib.sha256(certificate_raw).hexdigest(),
            "checkpoint_copy": str(args.checkpoint_copy),
            "checkpoint_copy_sha256": sha256(args.checkpoint_copy),
            "certificate_entries": len(certificate_lines),
        },
        "positive_18_control": {
            "status": "NOT_APPLICABLE_NO_K7",
            "vertices": 18,
            "K7_seeds": 0,
        },
    }
    atomic_json(args.report, report)
    print(
        f"complete: {dict(sorted(statuses.items()))}; "
        f"report SHA-256 {sha256(args.report)}",
        flush=True,
    )
    if statuses.get("INFRA_ERROR"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
