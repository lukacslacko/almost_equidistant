#!/usr/bin/env python3
"""Restartable full K7 positive-polynomial dual campaign.

This wrapper leaves the frozen mathematical locator unchanged.  It supplies
hash-pinned selection, deterministic 11-process execution, one atomic JSON
checkpoint per graph, exact independent checking of every emitted rational
certificate, and compact deterministic gzip decision/certificate archives.

``REJECTED`` is a mathematical claim only when ``marginal_dual_rejected`` is
true and its retained cover certificates pass the independent checker.
``SURVIVOR`` means only that this degree-bounded locator found no graph-level
rejection.  ``INFRA_ERROR`` makes no claim and makes the campaign exit 2.
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
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import d6_k7_positive_polynomial_dual as locator
import d6_k7_rank_reference as prior
import verify_d6_k7_positive_polynomial_dual as independent


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_INPUT_GRAPHS = 17_764
EXPECTED_DEPENDENCY_SHA256 = {
    "d6_k7_positive_polynomial_dual.py": (
        "190b6df17b90d6d19e4ce16829cd3dfebbfa87b9e7eb3d9d5e341b5a6469008b"
    ),
    "verify_d6_k7_positive_polynomial_dual.py": (
        "0f9c277cd7324cf887faa3d161a9328730833ceabca63c819f7197d423c6953a"
    ),
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_special_h_reference.py": (
        "e649633d7262941c9c4608ea694ae0a5d6f29c146d47060c3dae71a036a9d190"
    ),
    "d6_k7_positive_polynomial_dual.md": (
        "5fba6ad821528c7b1f2236282c77e65efd8d26dae99c3bfae34d8cf4f7f18e61"
    ),
    "test_d6_k7_positive_polynomial_dual.py": (
        "94cb523db058170d4ae546051b046942dec544220fcda1ceb17209460811a5d5"
    ),
}
QOS_CLASS_USER_INITIATED = 0x19
DECISION_FIELDS = (
    "ordinal",
    "index",
    "status",
    "applicable",
    "seeds",
    "covers",
    "enhanced_passing_covers",
    "strict_h_passing_covers",
    "dual_passing_covers",
    "strict_h_passing_cliques",
    "dual_failing_cliques",
    "strict_h_rejected",
    "dual_rejected",
    "marginal_dual_rejected",
    "error_type",
)
COUNT_FIELDS = (
    "seeds",
    "covers",
    "enhanced_passing_covers",
    "strict_h_passing_covers",
    "dual_passing_covers",
    "strict_h_passing_cliques",
    "dual_failing_cliques",
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


def atomic_deterministic_gzip(path: Path, content: bytes) -> None:
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


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
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


def verify_dependencies() -> dict[str, str]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCY_SHA256
    }
    if observed != EXPECTED_DEPENDENCY_SHA256:
        raise ValueError(
            f"frozen dual dependency hash mismatch: {observed} != "
            f"{EXPECTED_DEPENDENCY_SHA256}"
        )
    return observed


def resolve_manifest_source(manifest_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    local = manifest_path.parent / path
    return local if local.exists() else ROOT / path


def load_selection(
    input_path: Path,
    selection_report_path: Path,
    expected_selection_sha256: str,
) -> tuple[list[dict], dict]:
    input_hash = sha256(input_path)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"rank-survivor input hash mismatch: {input_hash}")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != EXPECTED_INPUT_GRAPHS:
        raise ValueError("rank-survivor input has an unexpected graph count")
    graph_by_index: dict[int, dict] = {}
    for graph in graphs:
        index = int(graph["index"])
        if index in graph_by_index:
            raise ValueError(f"duplicate rank-survivor graph {index}")
        prior.validate_graph(tuple(graph["adjacency"]))
        graph_by_index[index] = graph

    report_hash = sha256(selection_report_path)
    if report_hash != expected_selection_sha256:
        raise ValueError(
            f"selection report hash {report_hash} != explicit pin "
            f"{expected_selection_sha256}"
        )
    report = json.loads(selection_report_path.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind") != "d6_k7_positive_dual_selection"
        or report.get("input_sha256") != input_hash
    ):
        raise ValueError("selection report has the wrong schema/input")
    selected_values = report.get("selected_indices")
    if not isinstance(selected_values, list):
        raise ValueError("selection report has no selected_indices list")
    selected_indices = [int(value) for value in selected_values]
    if len(selected_indices) != len(set(selected_indices)):
        raise ValueError("selection report repeats an index")
    if report.get("selected_indices_sha256") != stable_hash(selected_indices):
        raise ValueError("selection index hash mismatch")
    if any(index not in graph_by_index for index in selected_indices):
        raise ValueError("selection contains a non-rank-survivor index")
    sources = report.get("source_artifacts")
    if not isinstance(sources, list) or not sources:
        raise ValueError("selection report has no source_artifacts")
    observed_sources = []
    for source in sources:
        path = resolve_manifest_source(selection_report_path, source["path"])
        actual = sha256(path)
        if actual != source["sha256"]:
            raise ValueError(f"selection source hash mismatch: {path}")
        observed_sources.append({
            "path": source["path"], "sha256": actual,
        })
    selected = [graph_by_index[index] for index in selected_indices]
    return selected, {
        "input": str(input_path),
        "input_sha256": input_hash,
        "rank_survivors": len(graphs),
        "selection_report": str(selection_report_path),
        "selection_report_sha256": report_hash,
        "selected": len(selected),
        "selected_indices_sha256": stable_hash(selected_indices),
        "source_artifacts": observed_sources,
    }


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


def evaluate_graph(graph: dict, multiplier_degree: int) -> dict:
    started = time.perf_counter()
    try:
        analysis = asdict(locator.analyze_graph(graph, multiplier_degree))
        # first_certificate duplicates the first item in the complete witness
        # list and is omitted from every production checkpoint.
        analysis.pop("first_certificate", None)
        return {
            "ordinal": graph["ordinal"],
            "index": int(graph["index"]),
            "status": "COMPLETE",
            "analysis": analysis,
            "elapsed_seconds": time.perf_counter() - started,
            "error_type": None,
            "error_message": None,
            "error_traceback": None,
        }
    except Exception as error:
        return {
            "ordinal": graph["ordinal"],
            "index": int(graph["index"]),
            "status": "INFRA_ERROR",
            "analysis": None,
            "elapsed_seconds": time.perf_counter() - started,
            "error_type": type(error).__name__,
            "error_message": str(error),
            "error_traceback": traceback.format_exc(),
        }


def evaluate_graph_star(arguments: tuple[dict, int]) -> dict:
    return evaluate_graph(*arguments)


def verify_witnesses_independently(graph: dict, analysis: dict) -> None:
    adjacency = tuple(graph["adjacency"])
    for witness in analysis["dual_failure_witnesses"]:
        nvertices = witness["nvertices"]
        if (
            not isinstance(nvertices, list)
            or len(nvertices) != len(set(nvertices))
            or any(not 0 <= int(vertex) < len(adjacency) for vertex in nvertices)
        ):
            raise ValueError("certificate witness has invalid induced vertices")
        graph_n = prior.induced_graph(adjacency, nvertices)
        certificate = witness["certificate"]
        clique_mask = sum(1 << vertex for vertex in certificate["clique"])
        independent.verify_clique_certificate(
            tuple(graph_n), clique_mask, certificate
        )


def validate_result(result: dict, graph: dict) -> None:
    if (
        result.get("ordinal") != graph["ordinal"]
        or result.get("index") != int(graph["index"])
        or result.get("status") not in ("COMPLETE", "INFRA_ERROR")
        or not isinstance(result.get("elapsed_seconds"), (int, float))
        or result["elapsed_seconds"] < 0
    ):
        raise ValueError("worker result identity/status/timing mismatch")
    if result["status"] == "INFRA_ERROR":
        if not result.get("error_type") or result.get("analysis") is not None:
            raise ValueError("invalid INFRA_ERROR worker result")
        return
    analysis = result.get("analysis")
    if not isinstance(analysis, dict) or int(analysis.get("index")) != int(
        graph["index"]
    ):
        raise ValueError("COMPLETE worker result has bad analysis")
    for field in COUNT_FIELDS:
        if not isinstance(analysis.get(field), int) or analysis[field] < 0:
            raise ValueError(f"invalid analysis counter {field}")
    for field in (
        "applicable", "strict_h_rejected", "dual_rejected",
        "marginal_dual_rejected",
    ):
        if not isinstance(analysis.get(field), bool):
            raise ValueError(f"invalid analysis flag {field}")
    if analysis["marginal_dual_rejected"] != (
        analysis["dual_rejected"] and not analysis["strict_h_rejected"]
    ):
        raise ValueError("marginal dual decision is inconsistent")
    witnesses = analysis.get("dual_failure_witnesses")
    if not isinstance(witnesses, list) or len(witnesses) != analysis[
        "dual_failing_cliques"
    ]:
        raise ValueError("certificate witness count mismatch")
    verify_witnesses_independently(graph, analysis)


def result_path(directory: Path, graph: dict) -> Path:
    return directory / "results" / (
        f"result_{graph['ordinal']:06d}_{int(graph['index']):07d}.json"
    )


def checkpoint_value(config_hash: str, result: dict) -> dict:
    return {
        "schema": 1,
        "config_sha256": config_hash,
        "completed_utc": utc_now(),
        "result": result,
    }


def load_checkpoint(path: Path, config_hash: str, graph: dict) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema") != 1 or value.get("config_sha256") != config_hash:
        raise ValueError(f"checkpoint schema/config mismatch: {path}")
    result = value.get("result")
    validate_result(result, graph)
    return result


def decision_row(result: dict) -> tuple[object, ...]:
    if result["status"] == "INFRA_ERROR":
        return (
            result["ordinal"], result["index"], "INFRA_ERROR",
            "", "", "", "", "", "", "", "", "", "", "",
            result["error_type"],
        )
    analysis = result["analysis"]
    decision = "REJECTED" if analysis["marginal_dual_rejected"] else "SURVIVOR"
    return (
        result["ordinal"], result["index"], decision,
        int(analysis["applicable"]),
        *(analysis[field] for field in COUNT_FIELDS),
        int(analysis["strict_h_rejected"]),
        int(analysis["dual_rejected"]),
        int(analysis["marginal_dual_rejected"]),
        "",
    )


def render_decisions(results: Sequence[dict]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(DECISION_FIELDS)
    writer.writerows(decision_row(result) for result in results)
    return buffer.getvalue().encode("ascii")


def render_certificates(results: Sequence[dict]) -> bytes:
    lines = []
    for result in results:
        if result["status"] != "COMPLETE":
            continue
        witnesses = result["analysis"]["dual_failure_witnesses"]
        if not witnesses:
            continue
        value = {
            "schema": 1,
            "ordinal": result["ordinal"],
            "index": result["index"],
            "dual_failure_witnesses": witnesses,
        }
        lines.append(json.dumps(value, sort_keys=True, separators=(",", ":")))
    return (("\n".join(lines) + "\n") if lines else "").encode("ascii")


def aggregate(results: Sequence[dict]) -> dict:
    complete = [result for result in results if result["status"] == "COMPLETE"]
    analyses = [result["analysis"] for result in complete]
    rejected = sorted(
        result["index"] for result in complete
        if result["analysis"]["marginal_dual_rejected"]
    )
    strict_rejected = sorted(
        result["index"] for result in complete
        if result["analysis"]["strict_h_rejected"]
    )
    return {
        "graphs": len(results),
        "complete": len(complete),
        "infra_errors": len(results) - len(complete),
        "marginal_dual_rejected": len(rejected),
        "survivors": len(complete) - len(rejected),
        "strict_h_rejected": len(strict_rejected),
        "marginal_dual_rejected_indices": rejected,
        "strict_h_rejected_indices": strict_rejected,
        "certificate_graphs": sum(
            bool(analysis["dual_failure_witnesses"]) for analysis in analyses
        ),
        "certificate_witnesses": sum(
            len(analysis["dual_failure_witnesses"]) for analysis in analyses
        ),
        "totals": {
            field: sum(analysis[field] for analysis in analyses)
            for field in COUNT_FIELDS
        },
    }


def run_campaign(
    graphs: list[dict],
    *,
    configuration: dict,
    workers: int,
    max_inflight: int,
    retry_infra_errors: bool,
    checkpoint_dir: Path,
    decisions_archive: Path,
    certificate_archive: Path,
    checkpoint_copy: Path,
    report_path: Path,
    pid_file: Path,
) -> tuple[dict, bool]:
    config_hash = stable_hash(configuration)
    campaign_path = checkpoint_dir / "campaign.json"
    campaign = {
        "schema": 1,
        "configuration": configuration,
        "config_sha256": config_hash,
    }
    with exclusive_lock(checkpoint_dir):
        if campaign_path.exists():
            existing = json.loads(campaign_path.read_text(encoding="utf-8"))
            if existing != campaign:
                raise RuntimeError("checkpoint directory belongs to another campaign")
        else:
            preexisting_outputs = [
                str(path) for path in (
                    decisions_archive, certificate_archive, checkpoint_copy,
                    report_path,
                ) if path.exists()
            ]
            if preexisting_outputs:
                raise RuntimeError(
                    "refusing to overwrite outputs without a compatible "
                    f"campaign checkpoint: {preexisting_outputs}"
                )
            atomic_json(campaign_path, campaign)
        completed: dict[int, dict] = {}
        pending: deque[dict] = deque()
        for graph in graphs:
            path = result_path(checkpoint_dir, graph)
            if not path.exists():
                pending.append(graph)
                continue
            result = load_checkpoint(path, config_hash, graph)
            if retry_infra_errors and result["status"] == "INFRA_ERROR":
                pending.append(graph)
            else:
                completed[graph["ordinal"]] = result

        started = time.perf_counter()
        initially_completed = len(completed)
        atomic_json(pid_file, {
            "schema": 1, "status": "RUNNING", "pid": os.getpid(),
            "started_utc": utc_now(), "checkpoint_dir": str(checkpoint_dir),
            "config_sha256": config_hash, "workers": workers,
        })
        atomic_json(checkpoint_dir / "run_state.json", {
            "schema": 1, "status": "RUNNING", "pid": os.getpid(),
            "updated_utc": utc_now(), "config_sha256": config_hash,
            "completed_graphs": len(completed),
        })

        executor = None
        if pending:
            executor = ProcessPoolExecutor(
                max_workers=workers, initializer=initialize_worker
            )
        active: dict = {}
        last_reported = len(completed)
        try:
            while pending or active:
                while pending and len(active) < max_inflight:
                    graph = pending.popleft()
                    assert executor is not None
                    future = executor.submit(
                        evaluate_graph_star,
                        (graph, configuration["multiplier_degree"]),
                    )
                    active[future] = graph
                if not active:
                    continue
                finished, _ = wait(active, return_when=FIRST_COMPLETED)
                for future in finished:
                    graph = active.pop(future)
                    try:
                        result = future.result()
                    except Exception as error:
                        result = {
                            "ordinal": graph["ordinal"],
                            "index": int(graph["index"]),
                            "status": "INFRA_ERROR",
                            "analysis": None,
                            "elapsed_seconds": 0.0,
                            "error_type": type(error).__name__,
                            "error_message": str(error),
                            "error_traceback": traceback.format_exc(),
                        }
                    validate_result(result, graph)
                    atomic_json(
                        result_path(checkpoint_dir, graph),
                        checkpoint_value(config_hash, result),
                    )
                    completed[graph["ordinal"]] = result
                if len(completed) - last_reported >= 32 or len(completed) == len(graphs):
                    print(
                        f"progress {len(completed)}/{len(graphs)}",
                        file=sys.stderr, flush=True,
                    )
                    atomic_json(checkpoint_dir / "progress.json", {
                        "schema": 1, "config_sha256": config_hash,
                        "completed_graphs": len(completed),
                        "total_graphs": len(graphs), "updated_utc": utc_now(),
                        "status_counts": dict(Counter(
                            result["status"] for result in completed.values()
                        )),
                    })
                    last_reported = len(completed)
            if executor is not None:
                executor.shutdown(wait=True)
        except BaseException:
            for future in active:
                future.cancel()
            if executor is not None:
                executor.shutdown(wait=True, cancel_futures=True)
            atomic_json(checkpoint_dir / "run_state.json", {
                "schema": 1, "status": "INTERRUPTED", "pid": os.getpid(),
                "updated_utc": utc_now(), "config_sha256": config_hash,
                "completed_graphs": len(completed),
            })
            raise

        results = [completed[ordinal] for ordinal in range(len(graphs))]
        decision_bytes = render_decisions(results)
        certificate_bytes = render_certificates(results)
        atomic_deterministic_gzip(decisions_archive, decision_bytes)
        atomic_deterministic_gzip(certificate_archive, certificate_bytes)
        checkpoint_entries = [
            [result_path(checkpoint_dir, graph).name,
             sha256(result_path(checkpoint_dir, graph))]
            for graph in graphs
        ]
        summary = aggregate(results)
        report = {
            "schema": 1,
            "description": (
                "Restartable exact degree-one K7 positive-polynomial dual "
                "campaign; numerical LP is only a locator and every retained "
                "certificate passed the independent rational checker."
            ),
            "claim": (
                "Only marginal_dual_rejected indices are new exact "
                "non-realizability conclusions from this layer."
            ),
            "configuration": configuration,
            "config_sha256": config_hash,
            "checkpoint_directory": str(checkpoint_dir),
            "checkpoint_result_index_sha256": stable_hash(checkpoint_entries),
            "summary": summary,
            "artifacts": {
                "decisions_archive": str(decisions_archive),
                "decisions_archive_sha256": sha256(decisions_archive),
                "decisions_uncompressed_sha256": hashlib.sha256(
                    decision_bytes
                ).hexdigest(),
                "decision_rows": len(results),
                "certificate_archive": str(certificate_archive),
                "certificate_archive_sha256": sha256(certificate_archive),
                "certificates_uncompressed_sha256": hashlib.sha256(
                    certificate_bytes
                ).hexdigest(),
            },
            "runtime": {
                "resume_wall_seconds": time.perf_counter() - started,
                "initially_completed_graphs": initially_completed,
                "newly_completed_graphs": len(results) - initially_completed,
                "sum_graph_seconds": sum(
                    result["elapsed_seconds"] for result in results
                ),
                "maximum_graph_seconds": max(
                    (result["elapsed_seconds"] for result in results), default=0.0
                ),
            },
            "environment": {
                "platform": platform.platform(),
                "machine": platform.machine(),
                "python": sys.version,
                "cpu_count": os.cpu_count(),
            },
        }
        atomic_json(report_path, report)
        final_status = (
            "COMPLETE" if summary["infra_errors"] == 0
            else "COMPLETE_WITH_INFRA_ERRORS"
        )
        final_manifest = {
            "schema": 1,
            "status": final_status,
            "config_sha256": config_hash,
            "completed_graphs": len(results),
            "checkpoint_result_index_sha256": stable_hash(checkpoint_entries),
            "report": str(report_path),
            "report_sha256": sha256(report_path),
            "artifacts": report["artifacts"],
            "finished_utc": utc_now(),
        }
        atomic_json(checkpoint_dir / "final_manifest.json", final_manifest)
        atomic_json(checkpoint_copy, final_manifest)
        atomic_json(pid_file, {
            "schema": 1, "status": final_status, "pid": os.getpid(),
            "finished_utc": utc_now(), "report": str(report_path),
            "report_sha256": sha256(report_path),
            "checkpoint": str(checkpoint_copy),
            "checkpoint_sha256": sha256(checkpoint_copy),
        })
        atomic_json(checkpoint_dir / "run_state.json", {
            "schema": 1, "status": final_status, "pid": os.getpid(),
            "updated_utc": utc_now(), "config_sha256": config_hash,
            "completed_graphs": len(results), "summary": summary,
        })
        return report, summary["infra_errors"] != 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(".runs/d6_k7_rank_survivors.json")
    )
    parser.add_argument("--selection-report", type=Path, required=True)
    parser.add_argument("--selection-report-sha256", required=True)
    parser.add_argument("--multiplier-degree", type=int, default=1)
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--max-inflight", type=int, default=22)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--retry-infra-errors", action="store_true")
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=Path(".runs/d6_k7_positive_dual_full"),
    )
    parser.add_argument(
        "--decisions-archive", type=Path,
        default=ROOT / "d6_k7_positive_dual_full_decisions.tsv.gz",
    )
    parser.add_argument(
        "--certificate-archive", type=Path,
        default=ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz",
    )
    parser.add_argument(
        "--checkpoint-copy", type=Path,
        default=ROOT / "d6_k7_positive_dual_full_checkpoint.json",
    )
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_positive_dual_full_report.json",
    )
    parser.add_argument(
        "--pid-file", type=Path,
        default=Path(".runs/d6_k7_positive_dual_full.pid.json"),
    )
    parser.add_argument(
        "--outer-command",
        default=os.environ.get("D6_OUTER_LAUNCH_COMMAND"),
    )
    args = parser.parse_args()
    if args.multiplier_degree != 1:
        parser.error("the frozen audited campaign currently pins degree one")
    if args.workers != 11:
        parser.error("production campaign pins exactly 11 workers")
    if args.max_inflight < args.workers:
        parser.error("--max-inflight must be at least the worker count")
    if args.start < 0 or (args.limit is not None and args.limit < 1):
        parser.error("invalid start/limit")

    dependencies = verify_dependencies()
    selected_full, selection = load_selection(
        args.input, args.selection_report, args.selection_report_sha256
    )
    stop = len(selected_full) if args.limit is None else args.start + args.limit
    if args.start >= len(selected_full) or stop > len(selected_full):
        parser.error("requested bounded stratum is outside the selection")
    selected = selected_full[args.start:stop]
    for ordinal, graph in enumerate(selected):
        graph = dict(graph)
        graph["ordinal"] = ordinal
        selected[ordinal] = graph
    indices = [int(graph["index"]) for graph in selected]
    python_command = [sys.executable, *sys.argv]
    recorded_outer_command = args.outer_command or shlex.join(python_command)
    configuration = {
        "schema": 1,
        "runner_source_sha256": sha256(Path(__file__).resolve()),
        "dependencies_sha256": dependencies,
        "selection": selection,
        "full_selection_graphs": len(selected_full),
        "full_selection_indices_sha256": selection["selected_indices_sha256"],
        "run_graphs": len(selected),
        "run_indices_sha256": stable_hash(indices),
        "start": args.start,
        "limit": args.limit,
        "multiplier_degree": args.multiplier_degree,
        "workers": args.workers,
        "max_inflight": args.max_inflight,
        "thread_limits": {
            name: "1" for name in (
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
            )
        },
        "git": git_provenance(),
        "command": {
            "python_argv": python_command,
            "python_command": shlex.join(python_command),
            "outer_command": recorded_outer_command,
            "outer_command_explicit": args.outer_command is not None,
        },
        "claim_semantics": {
            "REJECTED": (
                "Every eligible cover for at least one K7 seed fails after "
                "an exact rational degree-one positive-polynomial certificate."
            ),
            "SURVIVOR": "No rejection found by this layer; no realizability claim.",
            "INFRA_ERROR": "Worker/infrastructure failure; no claim.",
        },
    }
    print(
        f"positive dual: {len(selected)} graphs from offset {args.start}; "
        f"11 workers; checkpoint {args.checkpoint_dir}",
        file=sys.stderr, flush=True,
    )
    report, has_infra = run_campaign(
        selected,
        configuration=configuration,
        workers=args.workers,
        max_inflight=args.max_inflight,
        retry_infra_errors=args.retry_infra_errors,
        checkpoint_dir=args.checkpoint_dir,
        decisions_archive=args.decisions_archive,
        certificate_archive=args.certificate_archive,
        checkpoint_copy=args.checkpoint_copy,
        report_path=args.report,
        pid_file=args.pid_file,
    )
    print(
        f"complete: {report['summary']['marginal_dual_rejected']} rejected; "
        f"{report['summary']['survivors']} survivors; "
        f"{report['summary']['infra_errors']} infrastructure errors",
        file=sys.stderr, flush=True,
    )
    if has_infra:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
