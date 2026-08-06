#!/usr/bin/env python3
"""Restartable full degree-four K7 near-saturating tetrad campaign.

The input is the hash-pinned 12,839-graph survivor set of the completed
degree-one positive-polynomial campaign.  Each worker independently replays
the prior seed/cover/rank/strict-H quantifiers, consumes the already verified
degree-one cover certificates, and applies only degree-four rank-one Schur
tetrads to remaining covers without a saturating clique.

HiGHS is an untrusted locator.  ``REJECTED`` is emitted only when exact
rational tetrad identities eliminate every prior-passing cover for at least
one K7 seed.  ``SURVIVOR`` is not a realizability claim.  ``INFRA_ERROR``
makes no mathematical claim and makes finalization exit with status two.
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
import subprocess
import sys
import time
import traceback
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import d6_k7_positive_polynomial_dual as prior_dual
import d6_k7_rank_reference as prior
import d6_k7_rankone_tetrad_pilot as tetrad
import d6_k7_special_h_reference as strict_h


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_INPUT_GRAPHS = 17_764
EXPECTED_SELECTION_KIND = "d6_k7_rankone_tetrad_full_selection"
EXPECTED_SELECTED_GRAPHS = 12_839
EXPECTED_PRIOR_DECISIONS_SHA256 = (
    "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
)
EXPECTED_PRIOR_CERTIFICATES_SHA256 = (
    "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
)
EXPECTED_DEPENDENCY_SHA256 = {
    "d6_k7_rankone_tetrad_pilot.py": (
        "fe5a964b64817a7347776f4b1d5807b1f72d43bddfdd9934a0a324e13bfd260e"
    ),
    "d6_k7_rankone_tetrad_verify.py": (
        "ca12cfc05d745ebb0529100273a8ed8afddd1177bda435101d66786e88af71eb"
    ),
    "d6_k7_rankone_tetrad_test.py": (
        "007b924835a6298a97bc207ec9e17332a3347f1a95650758b34cba6f6b940bb4"
    ),
    "d6_k7_rankone_tetrad.md": (
        "0efd08d0b93b710d82d83c4359fd6e51ce29eb3aa75c8bb88df3e93f31f3e8d7"
    ),
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_special_h_reference.py": (
        "e649633d7262941c9c4608ea694ae0a5d6f29c146d47060c3dae71a036a9d190"
    ),
    "d6_k7_positive_polynomial_dual.py": (
        "190b6df17b90d6d19e4ce16829cd3dfebbfa87b9e7eb3d9d5e341b5a6469008b"
    ),
}
QOS_CLASS_USER_INITIATED = 0x19
DECISION_FIELDS = (
    "ordinal", "index", "status", "seeds", "covers",
    "enhanced_passing_covers", "strict_h_passing_covers",
    "prior_dual_failed_covers", "prior_dual_passing_covers",
    "prior_passing_no_saturating_clique_covers",
    "covers_with_near_clique", "covers_without_near_clique",
    "near_cliques_tested", "tetrad_failed_covers",
    "tetrad_passing_covers", "prior_dual_rejected",
    "tetrad_rejected", "marginal_tetrad_rejected", "error_type",
)
COUNT_FIELDS = tuple(
    field for field in DECISION_FIELDS[3:15]
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


@contextmanager
def exclusive_lock(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "campaign.lock"
    with path.open("a+", encoding="ascii") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                f"checkpoint directory is locked: {directory}"
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
        "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
        "porcelain_lines": status.splitlines(),
    }


def verify_dependencies() -> dict[str, str]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCY_SHA256
    }
    if observed != EXPECTED_DEPENDENCY_SHA256:
        raise ValueError(
            f"frozen dependency hash mismatch: {observed} != "
            f"{EXPECTED_DEPENDENCY_SHA256}"
        )
    return observed


def resolve_relative(manifest: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    local = manifest.parent / path
    return local if local.exists() else ROOT / path


def load_inputs(
    input_path: Path,
    selection_path: Path,
    expected_selection_sha256: str,
) -> tuple[list[dict], dict, dict[int, dict], dict[int, list[dict]]]:
    input_hash = sha256(input_path)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"rank input hash mismatch: {input_hash}")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != EXPECTED_INPUT_GRAPHS:
        raise ValueError("rank input graph count mismatch")
    graph_by_index = {}
    for graph in graphs:
        index = int(graph["index"])
        if index in graph_by_index:
            raise ValueError("rank input repeats an index")
        prior.validate_graph(tuple(graph["adjacency"]))
        graph_by_index[index] = graph

    selection_hash = sha256(selection_path)
    if selection_hash != expected_selection_sha256:
        raise ValueError(
            f"selection hash {selection_hash} != explicit pin "
            f"{expected_selection_sha256}"
        )
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if (
        selection.get("schema") != 1
        or selection.get("kind") != EXPECTED_SELECTION_KIND
        or selection["rank_input"]["sha256"] != input_hash
    ):
        raise ValueError("selection schema/input mismatch")
    selected_indices = [int(value) for value in selection["selected_indices"]]
    if (
        len(selected_indices) != EXPECTED_SELECTED_GRAPHS
        or len(selected_indices) != len(set(selected_indices))
        or selection["selected_indices_sha256"] != stable_hash(selected_indices)
    ):
        raise ValueError("selection index list mismatch")
    if any(index not in graph_by_index for index in selected_indices):
        raise ValueError("selection index absent from rank input")

    prior_decision_path = resolve_relative(
        selection_path, selection["prior_degree_one"]["decisions"]["path"]
    )
    prior_certificate_path = resolve_relative(
        selection_path,
        selection["prior_degree_one"]["certificates"]["path"],
    )
    if (
        sha256(prior_decision_path) != EXPECTED_PRIOR_DECISIONS_SHA256
        or sha256(prior_certificate_path) != EXPECTED_PRIOR_CERTIFICATES_SHA256
    ):
        raise ValueError("prior degree-one artifact hash mismatch")
    prior_decisions = {}
    with gzip.open(
        prior_decision_path, "rt", encoding="utf-8", newline=""
    ) as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            prior_decisions[int(row["index"])] = row
    prior_certificates = {}
    with gzip.open(prior_certificate_path, "rt", encoding="utf-8") as stream:
        for line in stream:
            entry = json.loads(line)
            prior_certificates[int(entry["index"])] = entry[
                "dual_failure_witnesses"
            ]
    selected = []
    for ordinal, index in enumerate(selected_indices):
        row = prior_decisions.get(index)
        if row is None or row["status"] != "SURVIVOR":
            raise ValueError("tetrad selection includes a prior nonsurvivor")
        selected.append({**graph_by_index[index], "ordinal": ordinal})
    provenance = {
        "input": str(input_path),
        "input_sha256": input_hash,
        "rank_survivors": len(graphs),
        "selection_report": str(selection_path),
        "selection_report_sha256": selection_hash,
        "selected": len(selected),
        "selected_indices_sha256": stable_hash(selected_indices),
        "prior_decisions": str(prior_decision_path),
        "prior_decisions_sha256": sha256(prior_decision_path),
        "prior_certificates": str(prior_certificate_path),
        "prior_certificates_sha256": sha256(prior_certificate_path),
    }
    return selected, provenance, prior_decisions, prior_certificates


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


@dataclass
class GraphResult:
    index: int
    seeds: int = 0
    covers: int = 0
    enhanced_passing_covers: int = 0
    strict_h_passing_covers: int = 0
    prior_dual_failed_covers: int = 0
    prior_dual_passing_covers: int = 0
    prior_passing_no_saturating_clique_covers: int = 0
    covers_with_near_clique: int = 0
    covers_without_near_clique: int = 0
    near_cliques_tested: int = 0
    tetrad_failed_covers: int = 0
    tetrad_passing_covers: int = 0
    prior_dual_rejected: bool = False
    tetrad_rejected: bool = False
    marginal_tetrad_rejected: bool = False
    tetrad_failure_witnesses: list[dict] = field(default_factory=list)
    elapsed_seconds: float = 0.0


def _prior_int(row: dict, name: str) -> int:
    return int(row[name])


def analyze_graph(
    graph: dict, prior_row: dict, prior_witnesses: Sequence[dict]
) -> GraphResult:
    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    result = GraphResult(index=int(graph["index"]))
    witness_by_key = {
        (tuple(item["seed"]), int(item["zmask"])): item
        for item in prior_witnesses
    }
    if len(witness_by_key) != len(prior_witnesses):
        raise ValueError("duplicate prior degree-one cover witness")
    consumed = set()
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    prior_rejected = tetrad_rejected = False
    strict_cliques = 0
    for seed_mask in prior.clique_masks(adj, 7):
        result.seeds += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(
            adj, seed_mask
        )
        total_term_rank = prior.matching_size(defects)
        seed_prior_passes = seed_tetrad_passes = 0
        for zmask in prior.eligible_covers(ladj, eligible):
            result.covers += 1
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            result.enhanced_passing_covers += 1
            nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            saturating = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            if any(
                strict_h.assess_saturating_clique(graph_n, mask).failed
                for mask in saturating
            ):
                continue
            strict_cliques += len(saturating)
            result.strict_h_passing_covers += 1
            key = (tuple(seed), zmask)
            old_witness = witness_by_key.get(key)
            if old_witness is not None:
                if old_witness["nvertices"] != nvertices:
                    raise ValueError("prior witness induced vertices mismatch")
                old_clique = sum(
                    1 << vertex
                    for vertex in old_witness["certificate"]["clique"]
                )
                if old_clique not in saturating:
                    raise ValueError("prior witness clique is not saturating")
                consumed.add(key)
                result.prior_dual_failed_covers += 1
                continue
            result.prior_dual_passing_covers += 1
            seed_prior_passes += 1
            if saturating:
                seed_tetrad_passes += 1
                result.tetrad_passing_covers += 1
                continue
            result.prior_passing_no_saturating_clique_covers += 1
            near_size = analysis.k_rank_upper - 1
            near_cliques = tuple(prior.clique_masks(
                graph_n, near_size
            )) if near_size > 0 else ()
            if not near_cliques:
                result.covers_without_near_clique += 1
                seed_tetrad_passes += 1
                result.tetrad_passing_covers += 1
                continue
            result.covers_with_near_clique += 1
            failure = None
            for near_mask in near_cliques:
                result.near_cliques_tested += 1
                certificate = tetrad.find_clique_certificate(
                    graph_n, near_mask, output_degree=4
                )
                if certificate is not None:
                    failure = {
                        "seed": seed,
                        "zmask": zmask,
                        "nvertices": nvertices,
                        "k_rank_upper": analysis.k_rank_upper,
                        "certificate": certificate,
                    }
                    break
            if failure is None:
                seed_tetrad_passes += 1
                result.tetrad_passing_covers += 1
            else:
                result.tetrad_failed_covers += 1
                result.tetrad_failure_witnesses.append(failure)
        if not seed_prior_passes:
            prior_rejected = True
        if not seed_tetrad_passes:
            tetrad_rejected = True
    if consumed != set(witness_by_key):
        raise ValueError("orphan prior degree-one cover witness")
    expected_prior = {
        "seeds": result.seeds,
        "covers": result.covers,
        "enhanced_passing_covers": result.enhanced_passing_covers,
        "strict_h_passing_covers": result.strict_h_passing_covers,
        "dual_passing_covers": result.prior_dual_passing_covers,
        "strict_h_passing_cliques": strict_cliques,
        "dual_failing_cliques": result.prior_dual_failed_covers,
    }
    for name, value in expected_prior.items():
        if _prior_int(prior_row, name) != value:
            raise ValueError(
                f"prior decision field {name}={prior_row[name]} != {value}"
            )
    if prior_row["status"] != "SURVIVOR" or prior_rejected:
        raise ValueError("full tetrad input is not a prior survivor")
    result.prior_dual_rejected = prior_rejected
    result.tetrad_rejected = tetrad_rejected
    result.marginal_tetrad_rejected = tetrad_rejected and not prior_rejected
    result.elapsed_seconds = time.perf_counter() - started
    return result


def evaluate_graph(payload: tuple[dict, dict, list[dict]]) -> dict:
    graph, prior_row, prior_witnesses = payload
    started = time.perf_counter()
    try:
        result = asdict(analyze_graph(graph, prior_row, prior_witnesses))
        status = "REJECTED" if result["marginal_tetrad_rejected"] else "SURVIVOR"
        return {
            "schema": 1,
            "ordinal": int(graph["ordinal"]),
            "index": int(graph["index"]),
            "status": status,
            "result": result,
        }
    except BaseException as error:  # noqa: BLE001 - infrastructure record
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
    return directory / f"graph_{ordinal:05d}_{index}.json"


def load_checkpoints(
    directory: Path, config_sha256: str, selected: Sequence[dict]
) -> dict[int, dict]:
    expected = {
        int(graph["ordinal"]): int(graph["index"]) for graph in selected
    }
    completed = {}
    for path in sorted(directory.glob("graph_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("config_sha256") != config_sha256:
            raise ValueError(f"checkpoint config mismatch: {path}")
        result = payload.get("record")
        ordinal = int(result["ordinal"])
        if expected.get(ordinal) != int(result["index"]):
            raise ValueError(f"checkpoint graph mismatch: {path}")
        if ordinal in completed:
            raise ValueError(f"duplicate checkpoint ordinal {ordinal}")
        completed[ordinal] = result
    return completed


def decision_row(record: dict) -> dict:
    result = record.get("result") or {}
    return {
        "ordinal": record["ordinal"],
        "index": record["index"],
        "status": record["status"],
        **{name: result.get(name, "") for name in COUNT_FIELDS},
        "prior_dual_rejected": result.get("prior_dual_rejected", ""),
        "tetrad_rejected": result.get("tetrad_rejected", ""),
        "marginal_tetrad_rejected": result.get(
            "marginal_tetrad_rejected", ""
        ),
        "error_type": record.get("error_type", ""),
    }


def render_tsv(records: Sequence[dict]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=DECISION_FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for record in records:
        writer.writerow(decision_row(record))
    return stream.getvalue().encode("utf-8")


def render_certificates(records: Sequence[dict]) -> bytes:
    lines = []
    for record in records:
        result = record.get("result") or {}
        witnesses = result.get("tetrad_failure_witnesses", [])
        if witnesses:
            lines.append(json.dumps({
                "schema": 1,
                "ordinal": record["ordinal"],
                "index": record["index"],
                "tetrad_failure_witnesses": witnesses,
            }, sort_keys=True, separators=(",", ":")))
    return (("\n".join(lines) + "\n") if lines else "").encode("utf-8")


def write_pid(path: Path, command: Sequence[str], config_sha256: str) -> None:
    atomic_json(path, {
        "schema": 1,
        "pid": os.getpid(),
        "started": utc_now(),
        "command": list(command),
        "config_sha256": config_sha256,
        "source_sha256": sha256(Path(__file__)),
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    parser.add_argument(
        "--selection-report", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_selection.json"),
    )
    parser.add_argument("--selection-report-sha256", required=True)
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--max-inflight", type=int, default=22)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=Path(".runs/d6_k7_rankone_tetrad_full"),
    )
    parser.add_argument(
        "--decisions-archive", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_decisions.tsv.gz"),
    )
    parser.add_argument(
        "--certificate-archive", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_certificates.jsonl.gz"),
    )
    parser.add_argument(
        "--checkpoint-copy", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_checkpoint.json"),
    )
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_report.json"),
    )
    parser.add_argument(
        "--pid-file", type=Path,
        default=Path(".runs/d6_k7_rankone_tetrad_full.pid.json"),
    )
    parser.add_argument("--outer-command")
    parser.add_argument("--progress-every", type=int, default=64)
    args = parser.parse_args()
    if args.workers <= 0 or args.max_inflight < args.workers:
        raise ValueError("workers must be positive and max-inflight >= workers")
    dependencies = verify_dependencies()
    selected, selection_provenance, prior_rows, prior_certificates = load_inputs(
        args.input, args.selection_report, args.selection_report_sha256
    )
    if args.limit is not None:
        if not 0 <= args.limit <= len(selected):
            raise ValueError("invalid limit")
        selected = selected[:args.limit]
    config = {
        "schema": 1,
        "input_sha256": selection_provenance["input_sha256"],
        "selection_report_sha256": selection_provenance[
            "selection_report_sha256"
        ],
        "selected_indices": [int(graph["index"]) for graph in selected],
        "prior_decisions_sha256": selection_provenance[
            "prior_decisions_sha256"
        ],
        "prior_certificates_sha256": selection_provenance[
            "prior_certificates_sha256"
        ],
        "output_degree": 4,
        "dependencies_sha256": dependencies,
    }
    config_sha256 = stable_hash(config)
    args.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(args.checkpoint_dir):
        existing_paths = list(args.checkpoint_dir.glob("graph_*.json"))
        if existing_paths and not args.resume:
            raise ValueError(
                "checkpoint directory is nonempty; pass --resume after audit"
            )
        completed = load_checkpoints(
            args.checkpoint_dir, config_sha256, selected
        ) if args.resume else {}
        command = [sys.executable, *sys.argv]
        write_pid(args.pid_file, command, config_sha256)
        started_wall = time.perf_counter()
        pending_graphs = [
            graph for graph in selected if int(graph["ordinal"]) not in completed
        ]
        submitted = completed_count = len(completed)
        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=initialize_worker
        ) as executor:
            inflight = {}
            queue = deque(pending_graphs)
            while queue or inflight:
                while queue and len(inflight) < args.max_inflight:
                    graph = queue.popleft()
                    index = int(graph["index"])
                    future = executor.submit(
                        evaluate_graph,
                        (
                            graph,
                            prior_rows[index],
                            prior_certificates.get(index, []),
                        ),
                    )
                    inflight[future] = graph
                    submitted += 1
                done, _pending = wait(
                    inflight, return_when=FIRST_COMPLETED
                )
                for future in done:
                    graph = inflight.pop(future)
                    record = future.result()
                    ordinal = int(graph["ordinal"])
                    path = checkpoint_path(
                        args.checkpoint_dir, ordinal, int(graph["index"])
                    )
                    atomic_json(path, {
                        "schema": 1,
                        "config_sha256": config_sha256,
                        "record": record,
                    })
                    completed[ordinal] = record
                    completed_count += 1
                    if (
                        completed_count % args.progress_every == 0
                        or completed_count == len(selected)
                    ):
                        elapsed = time.perf_counter() - started_wall
                        print(
                            f"[{utc_now()}] {completed_count}/{len(selected)} "
                            f"elapsed={elapsed:.1f}s",
                            flush=True,
                        )
        records = [completed[ordinal] for ordinal in range(len(selected))]
        if len(records) != len(selected):
            raise AssertionError("campaign finalization is incomplete")
        decision_bytes = render_tsv(records)
        certificate_bytes = render_certificates(records)
        atomic_deterministic_gzip(args.decisions_archive, decision_bytes)
        atomic_deterministic_gzip(
            args.certificate_archive, certificate_bytes
        )
        statuses = Counter(record["status"] for record in records)
        rejected_indices = [
            int(record["index"]) for record in records
            if record["status"] == "REJECTED"
        ]
        totals = {
            field: sum(
                int((record.get("result") or {}).get(field, 0))
                for record in records
            )
            for field in COUNT_FIELDS
        }
        summary = {
            "schema": 1,
            "config_sha256": config_sha256,
            "complete": len(records),
            "graphs": len(selected),
            "status_counts": dict(statuses),
            "marginal_tetrad_rejected": len(rejected_indices),
            "marginal_tetrad_rejected_indices": rejected_indices,
            "survivors": statuses.get("SURVIVOR", 0),
            "infra_errors": statuses.get("INFRA_ERROR", 0),
            "totals": totals,
        }
        atomic_json(args.checkpoint_copy, summary)
        report = {
            "schema": 1,
            "description": (
                "Full exact degree-four rank-one Schur tetrad campaign on "
                "the completed degree-one K7 residue."
            ),
            "claim": (
                "Only REJECTED rows whose rational identities pass the "
                "independent verifier are new exact non-realizability claims."
            ),
            "configuration": {
                "config_sha256": config_sha256,
                "workers": args.workers,
                "max_inflight": args.max_inflight,
                "limit": args.limit,
                "output_degree": 4,
                "selection": selection_provenance,
                "dependencies_sha256": dependencies,
                "runner_source_sha256": sha256(Path(__file__)),
                "command": {
                    "python_argv": command,
                    "python_command": shlex.join(command),
                    "outer_command": args.outer_command,
                },
                "git": git_provenance(),
            },
            "environment": {
                "python": sys.version,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "thread_limits": {
                    name: os.environ.get(name)
                    for name in (
                        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                        "NUMEXPR_NUM_THREADS",
                    )
                },
            },
            "runtime": {
                "finished": utc_now(),
                "wall_seconds": time.perf_counter() - started_wall,
                "resumed_records": len(records) - len(pending_graphs),
                "computed_records": len(pending_graphs),
            },
            "artifacts": {
                "decisions_archive": str(args.decisions_archive),
                "decisions_archive_sha256": sha256(args.decisions_archive),
                "decisions_uncompressed_sha256": hashlib.sha256(
                    decision_bytes
                ).hexdigest(),
                "certificate_archive": str(args.certificate_archive),
                "certificate_archive_sha256": sha256(args.certificate_archive),
                "certificates_uncompressed_sha256": hashlib.sha256(
                    certificate_bytes
                ).hexdigest(),
                "checkpoint_copy": str(args.checkpoint_copy),
                "checkpoint_copy_sha256": sha256(args.checkpoint_copy),
            },
            "summary": summary,
        }
        atomic_json(args.report, report)
        if summary["infra_errors"]:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
