#!/usr/bin/env python3
"""Restartable 11-worker full K6 normal-inertia campaign.

The mathematical evaluator is frozen in ``d6_k6_normal_inertia.py``.  This
wrapper selects the 1,097 survivors of the prior exact K6 campaign and adds
hash pinning, deterministic process parallelism, atomic prefix checkpoints,
and immutable decision/report artifacts.  Worker failures are infrastructure
errors and are never converted into mathematical rejections.
"""

from __future__ import annotations

import argparse
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
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterator, Sequence

from d6_k6_normal_inertia import evaluate_normal_graph


ROOT = Path(__file__).resolve().parent
EXPECTED_SELECTED = 1_097
EXPECTED_INPUT_SHA256 = (
    "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
)
EXPECTED_PRIOR_REPORT_SHA256 = (
    "ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc"
)
EXPECTED_DEPENDENCY_SHA256 = {
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "d6_k6_support_reference.py": (
        "d6481137ef49d88154882285660dd755eb7cb652ab87280f05c3792de9742577"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_normal_inertia.py": (
        "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023"
    ),
    "d6_k6_normal_inertia_test.py": (
        "4c107490befec2919c0a7eabe8f79f61eb8b11e56e40800cc630b5fdd3db7e12"
    ),
    "d6_k6_normal_coordinates.md": (
        "dd755349da45bbaca2f6c85dcdd5d56e4c08042f5d137cb349204e81b9686ca7"
    ),
    "d6_k6_normal_inertia_sample.json": (
        "3bc04d4541e11597302631fb30466af19d4c8b18333d4dffabca2f059b06128d"
    ),
}
DECISION_SCHEMA = "d6-k6-normal-inertia-full-decision-v1"
CHECKPOINT_SCHEMA = "d6-k6-normal-inertia-full-checkpoint-v1"
REPORT_SCHEMA = "d6-k6-normal-inertia-full-report-v1"


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


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


def atomic_jsonl(path: Path, records: Sequence[dict]) -> None:
    content = "".join(
        json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
        for record in records
    )
    atomic_text(path, content)


def read_jsonl(path: Path) -> list[dict]:
    records = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            try:
                value = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"invalid JSONL at {path}:{line_number}: {error}"
                ) from error
            if not isinstance(value, dict):
                raise ValueError(f"non-object JSONL record at {path}:{line_number}")
            records.append(value)
    return records


def deterministic_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with source.open("rb") as input_stream, temporary.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=9, fileobj=raw_output, mtime=0
        ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream, length=1 << 20)
        raw_output.flush()
        os.fsync(raw_output.fileno())
    os.replace(temporary, destination)
    fsync_directory(destination.parent)


@contextmanager
def exclusive_lock(checkpoint: Path) -> Iterator[None]:
    lock_path = checkpoint.with_name(checkpoint.name + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="ascii") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"checkpoint is locked: {checkpoint}") from error
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


def verify_dependencies(input_path: Path, prior_report_path: Path) -> dict[str, str]:
    observed = {
        input_path.name: sha256(input_path),
        prior_report_path.name: sha256(prior_report_path),
        **{name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCY_SHA256},
    }
    expected = {
        input_path.name: EXPECTED_INPUT_SHA256,
        prior_report_path.name: EXPECTED_PRIOR_REPORT_SHA256,
        **EXPECTED_DEPENDENCY_SHA256,
    }
    if observed != expected:
        raise ValueError(f"dependency hash mismatch: observed={observed}, expected={expected}")
    return observed


def select_residue(input_path: Path, prior_report_path: Path) -> list[dict]:
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    prior = json.loads(prior_report_path.read_text(encoding="utf-8"))
    if prior.get("input_graphs") != 1_098 or prior.get("graphs_surviving") != 1_097:
        raise ValueError("unexpected prior K6 population")
    rejected = {
        int(item["index"])
        for item in prior["graph_results"]
        if item["decision"]["rejected"]
    }
    if rejected != {461_363}:
        raise ValueError(f"unexpected prior rejection set: {sorted(rejected)}")
    records = [record for record in payload["graphs"] if record["index"] not in rejected]
    if len(records) != EXPECTED_SELECTED:
        raise ValueError(f"expected {EXPECTED_SELECTED} selected graphs, got {len(records)}")
    if len({record["index"] for record in records}) != len(records):
        raise ValueError("duplicate corpus index in selected records")
    return records


def _worker(record: dict) -> dict:
    decision = evaluate_normal_graph(record["adjacency"], scan_all=True)
    if not decision.applicable:
        raise ValueError(f"selected graph {record['index']} is not K6-only applicable")
    return {
        "schema": DECISION_SCHEMA,
        "index": int(record["index"]),
        "decision": decision.jsonable(),
    }


def validate_prefix(rows: Sequence[dict], selected: Sequence[dict]) -> None:
    if len(rows) > len(selected):
        raise ValueError("decision prefix is longer than the selected population")
    for position, row in enumerate(rows):
        if row.get("schema") != DECISION_SCHEMA:
            raise ValueError(f"decision schema mismatch at row {position}")
        if row.get("index") != selected[position]["index"]:
            raise ValueError(f"decision prefix index mismatch at row {position}")
        decision = row.get("decision")
        if not isinstance(decision, dict) or not isinstance(
            decision.get("rejected"), bool
        ):
            raise ValueError(f"malformed decision at row {position}")


def aggregate(rows: Sequence[dict]) -> dict:
    rejected = [row["index"] for row in rows if row["decision"]["rejected"]]
    sums = {}
    for name in (
        "seeds_checked",
        "impossible_seeds",
        "z0_subsets_considered",
        "z0_support_matchable",
        "bipartite_components_checked",
        "inertia_failures",
        "inertia_cache_entries",
        "inertia_cache_hits",
    ):
        sums[name] = sum(int(row["decision"][name]) for row in rows)
    return {
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(rows) - len(rejected),
        "rejected_indices": rejected,
        **sums,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k6_bipartite_rank_input.json"
    )
    parser.add_argument(
        "--prior-report",
        type=Path,
        default=ROOT / "d6_k6_bipartite_rank_report.json",
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=Path(".runs/d6_k6_normal_inertia_full_decisions.jsonl"),
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_full_decisions.jsonl.gz",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=Path(".runs/d6_k6_normal_inertia_full_checkpoint.json"),
    )
    parser.add_argument(
        "--checkpoint-copy",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_full_checkpoint.json",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_full_report.json",
    )
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--chunksize", type=int, default=1)
    parser.add_argument("--checkpoint-every", type=int, default=32)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--stop-after",
        type=int,
        default=0,
        help="test hook: checkpoint and pause after this many newly completed rows",
    )
    args = parser.parse_args()
    if min(args.workers, args.chunksize, args.checkpoint_every) < 1:
        parser.error("workers, chunksize, and checkpoint-every must be positive")
    if args.stop_after < 0:
        parser.error("stop-after cannot be negative")

    input_path = args.input.resolve()
    prior_report_path = args.prior_report.resolve()
    decisions_path = args.decisions.resolve()
    partial_path = decisions_path.with_name(decisions_path.name + ".partial")
    archive_path = args.archive.resolve()
    checkpoint_path = args.checkpoint.resolve()
    checkpoint_copy_path = args.checkpoint_copy.resolve()
    report_path = args.report.resolve()

    hashes = verify_dependencies(input_path, prior_report_path)
    selected = select_residue(input_path, prior_report_path)
    selection_hash = stable_hash(
        [{"index": row["index"], "adjacency": row["adjacency"]} for row in selected]
    )
    runner_hash = sha256(Path(__file__).resolve())
    compatibility = {
        "schema": CHECKPOINT_SCHEMA,
        "selected_graphs": len(selected),
        "selection_sha256": selection_hash,
        "dependency_sha256": hashes,
        "runner_sha256": runner_hash,
        "decision_schema": DECISION_SCHEMA,
        "workers": args.workers,
        "chunksize": args.chunksize,
        "checkpoint_every": args.checkpoint_every,
    }

    with exclusive_lock(checkpoint_path):
        previous = None
        rows: list[dict]
        if args.resume:
            if not checkpoint_path.exists():
                raise SystemExit("resume requires an existing checkpoint")
            previous = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            if previous.get("compatibility") != compatibility:
                raise SystemExit("resume checkpoint is incompatible with this campaign")
            if previous.get("status") == "COMPLETE":
                if not decisions_path.exists() or not archive_path.exists():
                    raise SystemExit("COMPLETE checkpoint is missing final artifacts")
                rows = read_jsonl(decisions_path)
                validate_prefix(rows, selected)
                if len(rows) != len(selected):
                    raise SystemExit("COMPLETE decision artifact has the wrong length")
                if previous.get("decisions_sha256") != sha256(decisions_path):
                    raise SystemExit("COMPLETE decision hash mismatch")
                if previous.get("archive_sha256") != sha256(archive_path):
                    raise SystemExit("COMPLETE archive hash mismatch")
                print(
                    f"already COMPLETE: {len(rows)}/{len(selected)} decisions; "
                    "zero recomputation",
                    flush=True,
                )
                return
            if not partial_path.exists():
                raise SystemExit("incomplete resume requires partial decisions")
            committed = int(previous.get("committed_rows", -1))
            loaded = read_jsonl(partial_path)
            if committed < 0 or committed > len(loaded):
                raise SystemExit("invalid committed-row count in checkpoint")
            rows = loaded[:committed]
            validate_prefix(rows, selected)
            if len(loaded) != committed:
                atomic_jsonl(partial_path, rows)
        else:
            existing = [
                path
                for path in (checkpoint_path, partial_path, decisions_path)
                if path.exists()
            ]
            if existing:
                raise SystemExit(
                    "campaign state already exists; use --resume after checking it: "
                    + ", ".join(map(str, existing))
                )
            rows = []
            atomic_jsonl(partial_path, rows)

        resume_count = int(previous.get("resume_count", 0)) + 1 if previous else 0
        started_at = previous.get("started_at", utc_now()) if previous else utc_now()

        def checkpoint(status: str, error: dict | None = None) -> dict:
            value = {
                "schema": CHECKPOINT_SCHEMA,
                "status": status,
                "started_at": started_at,
                "updated_at": utc_now(),
                "committed_rows": len(rows),
                "total_rows": len(selected),
                "resume_count": resume_count,
                "partial_decisions": str(partial_path),
                "compatibility": compatibility,
                "command": " ".join(sys.argv),
                "pid": os.getpid(),
            }
            if error is not None:
                value["infrastructure_error"] = error
            return value

        atomic_json(checkpoint_path, checkpoint("RUNNING"))
        print(
            f"K6 normal inertia: {len(selected)} graphs, {args.workers} workers, "
            f"resume position {len(rows)}",
            flush=True,
        )
        wall_start = time.time()
        new_rows = 0
        try:
            pending = selected[len(rows) :]
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                for record in executor.map(_worker, pending, chunksize=args.chunksize):
                    rows.append(record)
                    new_rows += 1
                    if len(rows) % args.checkpoint_every == 0:
                        atomic_jsonl(partial_path, rows)
                        atomic_json(checkpoint_path, checkpoint("RUNNING"))
                        print(f"completed {len(rows)}/{len(selected)}", flush=True)
                    if args.stop_after and new_rows >= args.stop_after:
                        break
        except BaseException as error:
            atomic_jsonl(partial_path, rows)
            detail = {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            }
            atomic_json(checkpoint_path, checkpoint("INFRA_ERROR", detail))
            raise

        atomic_jsonl(partial_path, rows)
        validate_prefix(rows, selected)
        if len(rows) != len(selected):
            atomic_json(checkpoint_path, checkpoint("PAUSED_TEST_HOOK"))
            print(
                f"paused after {len(rows)}/{len(selected)}; resume with --resume",
                flush=True,
            )
            return

        os.replace(partial_path, decisions_path)
        fsync_directory(decisions_path.parent)
        deterministic_gzip(decisions_path, archive_path)
        decisions_hash = sha256(decisions_path)
        archive_hash = sha256(archive_path)
        summary = aggregate(rows)
        final_checkpoint = checkpoint("COMPLETE")
        final_checkpoint.update(
            {
                "completed_at": utc_now(),
                "decisions_sha256": decisions_hash,
                "archive_sha256": archive_hash,
                "final_decisions": str(decisions_path),
                "archive": str(archive_path),
            }
        )
        atomic_json(checkpoint_path, final_checkpoint)
        atomic_json(checkpoint_copy_path, final_checkpoint)
        checkpoint_copy_hash = sha256(checkpoint_copy_path)

        report = {
            "schema": REPORT_SCHEMA,
            "status": "COMPLETE",
            "filter": "K6_bipartite_component_I_plus_adjacency_exact_inertia",
            "input_graphs": len(selected),
            **summary,
            "selection": {
                "description": "survivors of frozen exact K6 bipartite-rank report",
                "sha256": selection_hash,
            },
            "inputs": hashes,
            "sources": {
                "runner": Path(__file__).name,
                "runner_sha256": runner_hash,
                "decision_schema": DECISION_SCHEMA,
            },
            "artifacts": {
                "decisions": str(decisions_path),
                "decisions_sha256": decisions_hash,
                "archive": str(archive_path),
                "archive_sha256": archive_hash,
                "checkpoint": str(checkpoint_copy_path),
                "checkpoint_sha256": checkpoint_copy_hash,
            },
            "runtime": {
                "workers": args.workers,
                "chunksize": args.chunksize,
                "checkpoint_every": args.checkpoint_every,
                "resume_count": resume_count,
                "wall_seconds_this_invocation": time.time() - wall_start,
                "python": sys.version,
                "platform": platform.platform(),
                "command": " ".join(sys.argv),
                "started_at": started_at,
                "completed_at": utc_now(),
            },
            "git": git_provenance(),
            "arithmetic": (
                "finite graph operations and exact rational symmetric congruence; "
                "floating point is used only for elapsed time"
            ),
            "semantics": (
                "candidate nonedges are unconstrained and allowed defect "
                "coordinates may be zero"
            ),
        }
        atomic_json(report_path, report)
        print(
            json.dumps(
                {
                    "status": "COMPLETE",
                    "graphs_rejected": summary["graphs_rejected"],
                    "graphs_surviving": summary["graphs_surviving"],
                    "report": str(report_path),
                },
                sort_keys=True,
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
