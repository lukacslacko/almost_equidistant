#!/usr/bin/env python3
"""Restartable full K7 sparse-value campaign on support-layer survivors.

The mathematical evaluator is the independently audited sample profiler in
``profile_d6_k7_small_support_value.py``.  This wrapper leaves that package
byte-for-byte unchanged and supplies only production selection, deterministic
process parallelism, atomic prefix checkpoints, provenance, and a compact
decision archive.

No graph is silently skipped.  A worker exception leaves the last fsynced
prefix restartable and records ``INFRA_ERROR`` in the checkpoint; it is never
converted into a mathematical rejection.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterator, Sequence, TextIO

import d6_k7_support_propagation as propagation
import profile_d6_k7_small_support_value as sample_profiler


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_SUPPORT_ARCHIVE_SHA256 = (
    "3c35b228c4768b06f881dc9046329196b683d390f4763f83816c6a3464ecf938"
)
EXPECTED_SUPPORT_DECISIONS_SHA256 = (
    "117e057a79135fd48cac316ee4f1e08e11a2be84e00d9d10393bc17d4a518a31"
)
EXPECTED_SUPPORT_REPORT_SHA256 = (
    "5138d207543967b393baae9b14b47b49a276273cbc4468ee843856891917691e"
)
EXPECTED_SUPPORT_CHECKPOINT_SHA256 = (
    "daeee836f1f2c632656bbae812ef884cb5f2eca882fb2a5d5131f1c18edbdb5b"
)
EXPECTED_STRICT_H_REPORT_SHA256 = (
    "59be9c6a2d4cd2e8e8214d28d3198487846e58697046415ae1e15f4f69dc3add"
)
EXPECTED_DEPENDENCY_SHA256 = {
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_support_propagation.py": (
        "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
    ),
    "d6_k7_small_support_value.py": (
        "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f"
    ),
    "d6_k7_small_support_value.md": (
        "e8ec6a012f9ea824e9cddb76510840537dba0e603e207cd3220c1b06866cb42a"
    ),
    "profile_d6_k7_small_support_value.py": (
        "c8231ff47263b348e3daba2b32cd3eb693ff70f3db5b10aa808b1b5296589573"
    ),
    "test_d6_k7_small_support_value.py": (
        "6a2ee7a97844b5922a5467e1c3a8cf4c8a49240014b4154f45307335dfea5c44"
    ),
}
EXPECTED_RANK_SURVIVORS = 17_764
EXPECTED_SUPPORT_REJECTED = 1_536
EXPECTED_SELECTED = 16_228

BASE_COUNT_KEYS = (
    "seeds",
    "covers",
    "baseline_passing_covers",
    "labeled_z_families",
    "propagation_passing_families",
    "small_support_intersection_failures",
    "small_support_intersection_passes",
    "sparse_value_failures",
    "sparse_value_passes",
    "value_passing_covers",
)
PROPAGATION_REASONS = (
    "empty_propagated_mask",
    "disjoint_required_edge",
    "clique",
    "degree",
    "basis",
    "mask",
    "subspace_K",
    "component_B",
)
VALUE_REASONS = (
    "disjoint_required_small_support",
    "three_one_defects",
    "duplicate_one_defect",
    "three_same_two_defects",
    "two_defects_repeat_at_one_defect",
    "forbidden_two_defect_cycle",
    "parallel_two_defect_type_touches_another",
    "one_defects_joined_by_two_defects",
    "one_defect_component_branches",
    "incompatible_two_defect_branch_distance",
    "inconsistent_two_defect_branch_signs",
)
DECISION_FIELDS = (
    "index",
    "decision",
    "first_failing_seed",
    *BASE_COUNT_KEYS,
    *(f"propagation_failure_{reason}" for reason in PROPAGATION_REASONS),
    *(f"sparse_value_branch_{reason}" for reason in VALUE_REASONS),
)


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def decompressed_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ).stdout.strip()

    try:
        commit = git("rev-parse", "HEAD")
        branch = git("branch", "--show-current")
        status = git("status", "--porcelain=v1", "--untracked-files=all")
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}
    return {
        "available": True,
        "commit": commit,
        "branch": branch,
        "dirty": bool(status),
        "porcelain_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "porcelain_lines": status.splitlines(),
    }


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


def atomic_decisions(path: Path, rows: Sequence[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with temporary.open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(DECISION_FIELDS)
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)


def deterministic_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with source.open("rb") as input_stream, temporary.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=9,
            fileobj=raw_output,
            mtime=0,
        ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream, length=1 << 20)
        raw_output.flush()
        os.fsync(raw_output.fileno())
    os.replace(temporary, destination)
    fsync_directory(destination.parent)


@contextmanager
def exclusive_lock(checkpoint: Path):
    lock_path = checkpoint.with_name(checkpoint.name + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="ascii") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                f"checkpoint is locked by another process: {checkpoint}"
            ) from error
        stream.seek(0)
        stream.truncate()
        stream.write(f"pid={os.getpid()} started={utc_now()}\n")
        stream.flush()
        os.fsync(stream.fileno())
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def verify_dependencies() -> dict[str, str]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCY_SHA256
    }
    if observed != EXPECTED_DEPENDENCY_SHA256:
        raise ValueError(
            f"exact dependency hash mismatch: observed {observed}, "
            f"expected {EXPECTED_DEPENDENCY_SHA256}"
        )
    return observed


def open_decisions(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open("r", encoding="ascii", newline="")


def load_selection(
    input_path: Path,
    support_archive: Path,
    support_report_path: Path,
    support_checkpoint_path: Path,
) -> tuple[list[dict], dict]:
    hashes = {
        "input": sha256(input_path),
        "support_archive": sha256(support_archive),
        "support_decisions_uncompressed": decompressed_sha256(support_archive),
        "support_report": sha256(support_report_path),
        "support_checkpoint": sha256(support_checkpoint_path),
    }
    expected = {
        "input": EXPECTED_INPUT_SHA256,
        "support_archive": EXPECTED_SUPPORT_ARCHIVE_SHA256,
        "support_decisions_uncompressed": EXPECTED_SUPPORT_DECISIONS_SHA256,
        "support_report": EXPECTED_SUPPORT_REPORT_SHA256,
        "support_checkpoint": EXPECTED_SUPPORT_CHECKPOINT_SHA256,
    }
    if hashes != expected:
        raise ValueError(f"selection artifact hash mismatch: {hashes} != {expected}")

    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != EXPECTED_RANK_SURVIVORS:
        raise ValueError("rank-survivor input has an unexpected graph count")
    indices = [int(graph["index"]) for graph in graphs]
    if len(indices) != len(set(indices)):
        raise ValueError("rank-survivor input has duplicate indices")

    support_report = json.loads(support_report_path.read_text(encoding="utf-8"))
    if (
        int(support_report.get("graphs", -1)) != EXPECTED_RANK_SURVIVORS
        or int(support_report.get("refined_rejected", -1))
        != EXPECTED_SUPPORT_REJECTED
        or int(support_report.get("refined_survivors", -1)) != EXPECTED_SELECTED
        or support_report.get("input_sha256") != EXPECTED_INPUT_SHA256
        or support_report.get("decisions_sha256")
        != EXPECTED_SUPPORT_DECISIONS_SHA256
        or support_report.get("source_sha256")
        != EXPECTED_DEPENDENCY_SHA256["d6_k7_support_propagation.py"]
        or support_report.get("reference_sha256")
        != EXPECTED_DEPENDENCY_SHA256["d6_k7_rank_reference.py"]
    ):
        raise ValueError("support report metadata is inconsistent")
    report_rejected = set(map(int, support_report["refined_rejected_indices"]))

    support_checkpoint = json.loads(
        support_checkpoint_path.read_text(encoding="utf-8")
    )
    compatibility = support_checkpoint.get("compatibility", {})
    if (
        support_checkpoint.get("status") != "COMPLETE"
        or int(support_checkpoint.get("committed_rows", -1))
        != EXPECTED_RANK_SURVIVORS
        or support_checkpoint.get("decisions_sha256")
        != EXPECTED_SUPPORT_DECISIONS_SHA256
        or compatibility.get("input_sha256") != EXPECTED_INPUT_SHA256
        or compatibility.get("source_sha256")
        != EXPECTED_DEPENDENCY_SHA256["d6_k7_support_propagation.py"]
        or compatibility.get("reference_sha256")
        != EXPECTED_DEPENDENCY_SHA256["d6_k7_rank_reference.py"]
    ):
        raise ValueError("support checkpoint metadata is inconsistent")

    selected: list[dict] = []
    observed_rejected: set[int] = set()
    with open_decisions(support_archive) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != propagation.DECISION_FIELDS:
            raise ValueError("support decision schema mismatch")
        rows = list(reader)
    if len(rows) != EXPECTED_RANK_SURVIVORS:
        raise ValueError("support decision archive has an unexpected row count")
    for position, (graph, row) in enumerate(zip(graphs, rows)):
        index = int(row["index"])
        if index != int(graph["index"]):
            raise ValueError(f"support row/input mismatch at position {position}")
        decision = row["refined_decision"]
        if decision == "SURVIVOR":
            selected.append(graph)
        elif decision == "REJECTED":
            observed_rejected.add(index)
        else:
            raise ValueError(f"bad support decision {decision!r}")
    if observed_rejected != report_rejected:
        raise ValueError("support archive rejection set disagrees with report")
    if len(selected) != EXPECTED_SELECTED:
        raise ValueError(f"unexpected support survivors: {len(selected)}")
    selected_indices = [int(graph["index"]) for graph in selected]
    return selected, {
        "hashes": hashes,
        "rank_survivors": len(graphs),
        "support_rejected": len(observed_rejected),
        "selected": len(selected),
        "selected_indices_sha256": stable_hash(selected_indices),
    }


def load_strict_h(path: Path, rank_indices: set[int]) -> tuple[set[int], dict]:
    actual = sha256(path)
    if actual != EXPECTED_STRICT_H_REPORT_SHA256:
        raise ValueError(f"strict-H report hash mismatch: {actual}")
    report = json.loads(path.read_text(encoding="utf-8"))
    rejected = set(map(int, report.get("decisions", {}).get("strict_H_rejected", ())))
    survivors = set(map(int, report.get("decisions", {}).get("strict_H_survivors", ())))
    if (
        rejected & survivors
        or rejected | survivors != rank_indices
        or len(rejected) != 603
        or len(survivors) != 17_161
    ):
        raise ValueError("strict-H report is not the expected rank-residue partition")
    return rejected, {
        "path": str(path),
        "sha256": actual,
        "rejected": len(rejected),
        "survivors": len(survivors),
    }


def initialize_worker() -> None:
    propagation.set_worker_qos()


def evaluate_graph(graph: dict) -> dict:
    return sample_profiler.analyze_graph(graph)


def flatten_result(result: dict) -> tuple[object, ...]:
    counts = result["counts"]
    allowed = set(BASE_COUNT_KEYS)
    allowed.update(f"propagation_failure:{reason}" for reason in PROPAGATION_REASONS)
    allowed.update(f"sparse_value_branch:{reason}" for reason in VALUE_REASONS)
    unknown = set(counts) - allowed
    if unknown:
        raise ValueError(f"worker returned unknown counters: {sorted(unknown)}")
    decision = result["decision"]
    if decision not in ("REJECTED", "SURVIVOR"):
        raise ValueError(f"worker returned bad decision {decision!r}")
    return (
        int(result["index"]),
        decision,
        int(result["first_failing_seed"]),
        *(int(counts.get(key, 0)) for key in BASE_COUNT_KEYS),
        *(int(counts.get(f"propagation_failure:{reason}", 0))
          for reason in PROPAGATION_REASONS),
        *(int(counts.get(f"sparse_value_branch:{reason}", 0))
          for reason in VALUE_REASONS),
    )


def validate_row(row: Sequence[object], expected_index: int) -> tuple[str, ...]:
    values = tuple(map(str, row))
    if len(values) != len(DECISION_FIELDS):
        raise ValueError("decision row has a wrong field count")
    if int(values[0]) != expected_index:
        raise ValueError(f"decision row index {values[0]} != {expected_index}")
    if values[1] not in ("REJECTED", "SURVIVOR"):
        raise ValueError(f"bad decision {values[1]!r}")
    for value in (values[0], *values[2:]):
        if int(value) < 0:
            raise ValueError("negative integer in decision row")
    if (values[1] == "REJECTED") != (int(values[2]) != 0):
        raise ValueError("decision and failing-seed witness disagree")
    return values


def load_completed(
    partial: Path,
    expected_indices: Sequence[int],
    committed_rows: int,
) -> list[tuple[str, ...]]:
    if committed_rows < 0 or committed_rows > len(expected_indices):
        raise ValueError("checkpoint committed-row count is out of range")
    if not partial.exists():
        if committed_rows:
            raise ValueError("checkpoint names rows but partial TSV is absent")
        return []
    with partial.open("r", encoding="ascii", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, None)
        if header != list(DECISION_FIELDS):
            raise ValueError("partial decision schema mismatch")
        rows = [tuple(row) for row in reader]
    if len(rows) < committed_rows:
        raise ValueError("partial TSV is shorter than its committed prefix")
    completed = [
        validate_row(row, expected_indices[position])
        for position, row in enumerate(rows[:committed_rows])
    ]
    if len(rows) != committed_rows:
        atomic_decisions(partial, completed)
    return completed


def batches(values: Sequence[dict], size: int) -> Iterator[Sequence[dict]]:
    for start in range(0, len(values), size):
        yield values[start:start + size]


def summarize_rows(rows: Sequence[Sequence[object]]) -> tuple[dict, set[int]]:
    positions = {name: DECISION_FIELDS.index(name) for name in DECISION_FIELDS}
    rejected = {int(row[0]) for row in rows if row[1] == "REJECTED"}
    totals = {
        name: sum(int(row[position]) for row in rows)
        for name, position in positions.items()
        if name not in ("index", "decision", "first_failing_seed")
    }
    return {
        "graphs": len(rows),
        "rejected": len(rejected),
        "survivors": len(rows) - len(rejected),
        "rejected_indices": sorted(rejected),
        "survivor_indices_sha256": stable_hash(
            [int(row[0]) for row in rows if row[1] == "SURVIVOR"]
        ),
        "totals": totals,
    }, rejected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(".runs/d6_k7_rank_survivors.json")
    )
    parser.add_argument(
        "--support-decisions", type=Path,
        default=ROOT / "d6_k7_support_rank_survivors_decisions.tsv.gz",
    )
    parser.add_argument(
        "--support-report", type=Path,
        default=ROOT / "d6_k7_support_rank_survivors_report.json",
    )
    parser.add_argument(
        "--support-checkpoint", type=Path,
        default=ROOT / "d6_k7_support_rank_survivors_checkpoint.json",
    )
    parser.add_argument(
        "--strict-h-report", type=Path,
        default=ROOT / "d6_k7_strict_h_full_report.json",
    )
    parser.add_argument(
        "--decisions", type=Path,
        default=Path(".runs/d6_k7_small_support_value_full_decisions.tsv"),
    )
    parser.add_argument(
        "--archive", type=Path,
        default=ROOT / "d6_k7_small_support_value_full_decisions.tsv.gz",
    )
    parser.add_argument(
        "--checkpoint", type=Path,
        default=Path(".runs/d6_k7_small_support_value_full_checkpoint.json"),
    )
    parser.add_argument(
        "--checkpoint-copy", type=Path,
        default=ROOT / "d6_k7_small_support_value_full_checkpoint.json",
    )
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_small_support_value_full_report.json",
    )
    parser.add_argument(
        "--pid-file", type=Path,
        default=Path(".runs/d6_k7_small_support_value_full.pid.json"),
    )
    parser.add_argument(
        "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)
    )
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--map-chunksize", type=int, default=1)
    parser.add_argument("--checkpoint-every", type=int, default=32)
    parser.add_argument("--progress-every", type=int, default=64)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--stop-after-new", type=int,
        help="test hook: checkpoint and pause after this many new rows",
    )
    parser.add_argument(
        "--outer-command",
        default=os.environ.get("D6_OUTER_LAUNCH_COMMAND"),
        help="record the taskpolicy/caffeinate wrapper used to launch Python",
    )
    args = parser.parse_args()
    for name in (
        "workers", "batch_size", "map_chunksize", "checkpoint_every",
        "progress_every",
    ):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.start < 0:
        parser.error("--start must be nonnegative")
    if args.stop_after_new is not None and args.stop_after_new < 1:
        parser.error("--stop-after-new must be positive")

    dependencies = verify_dependencies()
    selected_full, selection = load_selection(
        args.input,
        args.support_decisions,
        args.support_report,
        args.support_checkpoint,
    )
    rank_payload = json.loads(args.input.read_text(encoding="utf-8"))
    rank_indices = {int(graph["index"]) for graph in rank_payload["graphs"]}
    strict_rejected, strict_provenance = load_strict_h(
        args.strict_h_report, rank_indices
    )
    stop = (
        args.start + args.limit
        if args.limit is not None
        else len(selected_full)
    )
    if stop > len(selected_full):
        parser.error("the requested selection slice exceeds the full input")
    selected = selected_full[args.start:stop]
    if not selected:
        parser.error("the requested selection slice is empty")
    indices = [int(graph["index"]) for graph in selected]
    source_hash = sha256(Path(__file__))
    git_state = git_provenance()
    command = [sys.executable, *sys.argv]
    outer_command = args.outer_command or " ".join(command)
    compatibility = {
        "schema": 1,
        "runner_source_sha256": source_hash,
        "dependencies_sha256": dependencies,
        "selection_artifact_hashes": selection["hashes"],
        "full_selected_indices_sha256": selection["selected_indices_sha256"],
        "run_indices_sha256": stable_hash(indices),
        "run_graphs": len(selected),
        "start": args.start,
        "limit": args.limit,
        "decision_fields": list(DECISION_FIELDS),
    }
    partial = args.decisions.with_name(args.decisions.name + ".partial")

    with exclusive_lock(args.checkpoint):
        previous = None
        if args.resume:
            if not args.checkpoint.exists():
                raise SystemExit(f"resume checkpoint is absent: {args.checkpoint}")
            previous = json.loads(args.checkpoint.read_text(encoding="utf-8"))
            if previous.get("compatibility") != compatibility:
                raise SystemExit("resume checkpoint belongs to another campaign")
            completed = load_completed(
                partial, indices, int(previous["committed_rows"])
            )
        else:
            existing = [
                path for path in (
                    partial, args.checkpoint, args.decisions, args.archive,
                    args.checkpoint_copy, args.report,
                ) if path.exists()
            ]
            if existing:
                raise SystemExit(
                    "refusing to overwrite existing campaign artifacts: "
                    + ", ".join(map(str, existing))
                )
            completed = []
            atomic_decisions(partial, ())

        initial_started = (
            previous["initial_started_utc"] if previous else utc_now()
        )
        initial_command = previous["initial_command"] if previous else command
        initial_outer_command = (
            previous["initial_outer_command"] if previous else outer_command
        )
        resume_count = int(previous.get("resume_count", 0)) + 1 if previous else 0
        session_started_utc = utc_now()
        session_started_unix = time.time()
        session_started = time.perf_counter()

        def checkpoint_value(
            status: str,
            committed_rows: int,
            *,
            error: dict | None = None,
            decisions_sha: str | None = None,
            archive_sha: str | None = None,
        ) -> dict:
            value = {
                "schema": 1,
                "status": status,
                "compatibility": compatibility,
                "committed_rows": committed_rows,
                "partial_decisions": str(partial),
                "final_decisions": str(args.decisions),
                "archive": str(args.archive),
                "initial_started_utc": initial_started,
                "initial_command": initial_command,
                "initial_outer_command": initial_outer_command,
                "last_session_started_utc": session_started_utc,
                "last_command": command,
                "last_outer_command": outer_command,
                "resume_count": resume_count,
                "workers": args.workers,
                "batch_size": args.batch_size,
                "map_chunksize": args.map_chunksize,
                "checkpoint_every": args.checkpoint_every,
                "progress_every": args.progress_every,
                "git": git_state,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "python_version": sys.version,
            }
            if error is not None:
                value["error"] = error
            if decisions_sha is not None:
                value["decisions_sha256"] = decisions_sha
            if archive_sha is not None:
                value["archive_sha256"] = archive_sha
            return value

        atomic_json(
            args.checkpoint,
            checkpoint_value("RUNNING", len(completed)),
        )
        atomic_json(
            args.pid_file,
            {
                "schema": 1,
                "status": "RUNNING",
                "pid": os.getpid(),
                "started_utc": session_started_utc,
                "started_unix": session_started_unix,
                "command": command,
                "outer_command": outer_command,
                "checkpoint": str(args.checkpoint),
                "workers": args.workers,
                "compatibility": compatibility,
            },
        )
        print(
            f"PID {os.getpid()}; {len(selected)} selected graphs; "
            f"{args.workers} workers; resume position {len(completed)}",
            file=sys.stderr,
            flush=True,
        )

        all_rows: list[Sequence[object]] = list(completed)
        durable_rows = len(completed)
        newly_completed = 0
        paused = False
        try:
            remaining = selected[len(completed):]
            with partial.open("a", encoding="ascii", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                with ProcessPoolExecutor(
                    max_workers=args.workers,
                    initializer=initialize_worker,
                ) as executor:
                    for group in batches(remaining, args.batch_size):
                        for result in executor.map(
                            evaluate_graph, group, chunksize=args.map_chunksize
                        ):
                            row = flatten_result(result)
                            expected_index = indices[len(all_rows)]
                            validate_row(row, expected_index)
                            writer.writerow(row)
                            all_rows.append(row)
                            newly_completed += 1
                            done = len(all_rows)
                            must_sync = (
                                done % args.checkpoint_every == 0
                                or done == len(selected)
                                or (
                                    args.stop_after_new is not None
                                    and newly_completed >= args.stop_after_new
                                )
                            )
                            if must_sync:
                                stream.flush()
                                os.fsync(stream.fileno())
                                atomic_json(
                                    args.checkpoint,
                                    checkpoint_value("RUNNING", done),
                                )
                                durable_rows = done
                            if done % args.progress_every == 0 or done == len(selected):
                                elapsed = time.perf_counter() - session_started
                                rate = newly_completed / max(elapsed, 1e-9)
                                print(
                                    f"progress {done}/{len(selected)} "
                                    f"({done/len(selected):.2%}); "
                                    f"{rate:.2f} graph/s; {elapsed:.1f}s",
                                    file=sys.stderr,
                                    flush=True,
                                )
                            if (
                                args.stop_after_new is not None
                                and newly_completed >= args.stop_after_new
                            ):
                                paused = True
                                break
                        if paused:
                            break
        except BaseException as error:
            status = "ABORT" if isinstance(error, KeyboardInterrupt) else "INFRA_ERROR"
            error_record = {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            }
            atomic_json(
                args.checkpoint,
                checkpoint_value(status, durable_rows, error=error_record),
            )
            atomic_json(
                args.pid_file,
                {
                    "schema": 1,
                    "status": status,
                    "pid": os.getpid(),
                    "finished_utc": utc_now(),
                    "checkpoint": str(args.checkpoint),
                    "error": error_record,
                },
            )
            raise

        if paused:
            atomic_json(
                args.checkpoint,
                checkpoint_value("PAUSED_TEST_HOOK", len(all_rows)),
            )
            atomic_json(
                args.pid_file,
                {
                    "schema": 1,
                    "status": "PAUSED_TEST_HOOK",
                    "pid": os.getpid(),
                    "finished_utc": utc_now(),
                    "checkpoint": str(args.checkpoint),
                    "committed_rows": len(all_rows),
                },
            )
            print(
                f"paused by test hook after {newly_completed} new rows; "
                f"resume with --resume",
                file=sys.stderr,
                flush=True,
            )
            return

        if len(all_rows) != len(selected):
            raise RuntimeError(
                f"completed {len(all_rows)} rows for {len(selected)} graphs"
            )
        os.replace(partial, args.decisions)
        fsync_directory(args.decisions.parent)
        deterministic_gzip(args.decisions, args.archive)
        decisions_hash = sha256(args.decisions)
        archive_hash = sha256(args.archive)
        final_checkpoint = checkpoint_value(
            "COMPLETE",
            len(all_rows),
            decisions_sha=decisions_hash,
            archive_sha=archive_hash,
        )
        atomic_json(args.checkpoint, final_checkpoint)
        atomic_json(args.checkpoint_copy, final_checkpoint)
        checkpoint_hash = sha256(args.checkpoint_copy)

        summary, value_rejected = summarize_rows(all_rows)
        all_support_selected = {
            int(graph["index"]) for graph in selected_full
        }
        strict_in_support = strict_rejected & all_support_selected
        strict_in_run = strict_rejected & set(indices)
        combined_in_run = strict_in_run | value_rejected
        report = {
            "schema": 1,
            "description": (
                "Exact K7 one/two-defect sparse-value campaign on the "
                "labeled-support survivors."
            ),
            "configuration": compatibility,
            "configuration_sha256": stable_hash(compatibility),
            "selection": selection,
            "strict_H": strict_provenance,
            "run": {
                "graphs": len(selected),
                "start": args.start,
                "limit": args.limit,
                "workers": args.workers,
                "batch_size": args.batch_size,
                "map_chunksize": args.map_chunksize,
                "checkpoint_every": args.checkpoint_every,
                "command": command,
                "outer_command": outer_command,
                "initial_started_utc": initial_started,
                "finished_utc": utc_now(),
                "resume_count": resume_count,
                "last_session_wall_seconds": time.perf_counter() - session_started,
            },
            "environment": {
                "git": git_state,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "python_version": sys.version,
            },
            "artifacts": {
                "decisions": str(args.decisions),
                "decisions_sha256": decisions_hash,
                "archive": str(args.archive),
                "archive_sha256": archive_hash,
                "checkpoint": str(args.checkpoint_copy),
                "checkpoint_sha256": checkpoint_hash,
            },
            "decisions": summary,
            "comparison_with_strict_H": {
                "strict_H_rejected_among_all_support_survivors": len(
                    strict_in_support
                ),
                "strict_H_rejected_in_run": len(strict_in_run),
                "value_rejected_in_run": len(value_rejected),
                "intersection_in_run": len(strict_in_run & value_rejected),
                "strict_H_only_in_run": len(strict_in_run - value_rejected),
                "value_only_in_run": len(value_rejected - strict_in_run),
                "union_rejected_in_run": len(combined_in_run),
                "union_survivors_in_run": len(selected) - len(combined_in_run),
                "intersection_indices": sorted(strict_in_run & value_rejected),
            },
        }
        atomic_json(args.report, report)
        report_hash = sha256(args.report)
        atomic_json(
            args.pid_file,
            {
                "schema": 1,
                "status": "COMPLETE",
                "pid": os.getpid(),
                "finished_utc": utc_now(),
                "checkpoint": str(args.checkpoint_copy),
                "checkpoint_sha256": checkpoint_hash,
                "archive": str(args.archive),
                "archive_sha256": archive_hash,
                "report": str(args.report),
                "report_sha256": report_hash,
                "graphs": len(selected),
                "rejected": summary["rejected"],
                "survivors": summary["survivors"],
            },
        )
        print(
            f"complete: {len(selected)} graphs, {summary['rejected']} rejected, "
            f"{summary['survivors']} survivors; report SHA-256 {report_hash}",
            file=sys.stderr,
            flush=True,
        )


if __name__ == "__main__":
    main()
