#!/usr/bin/env python3
"""Independent restartable verifier for the full K7 tetrad campaign.

This program imports neither the numerical tetrad locator nor the production
runner.  For every selected graph it rebuilds every K7 seed, eligible cover,
prior exact rank/strict-H decision, saturating and near-saturating clique,
then independently expands every retained degree-one and tetrad rational
polynomial identity.  It compares the complete reconstructed quantifiers to
the compact production decision row.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import fcntl
import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as prior
import d6_k7_rankone_tetrad_verify as tetrad_checker
import d6_k7_special_h_reference as strict_h
import verify_d6_k7_positive_polynomial_dual as prior_checker


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_INPUT_GRAPHS = 17_764
EXPECTED_SELECTION_SHA256 = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
EXPECTED_SELECTED_GRAPHS = 12_839
EXPECTED_PRIOR_DECISIONS_SHA256 = (
    "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
)
EXPECTED_PRIOR_CERTIFICATES_SHA256 = (
    "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
)
EXPECTED_DEPENDENCY_SHA256 = {
    "d6_k7_rankone_tetrad_verify.py": (
        "ca12cfc05d745ebb0529100273a8ed8afddd1177bda435101d66786e88af71eb"
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
}
QOS_CLASS_USER_INITIATED = 0x19


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


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    content = json.dumps(value, indent=2, sort_keys=True) + "\n"
    with temporary.open("w", encoding="utf-8") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)


@contextmanager
def exclusive_lock(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "verification.lock"
    with path.open("a+", encoding="ascii") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                f"verification directory is locked: {directory}"
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


def verify_dependencies() -> dict[str, str]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCY_SHA256
    }
    if observed != EXPECTED_DEPENDENCY_SHA256:
        raise ValueError("independent verifier dependency hash mismatch")
    return observed


def resolve_relative(manifest: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    local = manifest.parent / path
    return local if local.exists() else ROOT / path


def read_decisions(path: Path) -> tuple[list[dict], dict[int, dict]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    mapping = {int(row["index"]): row for row in rows}
    if len(mapping) != len(rows):
        raise ValueError(f"decision archive repeats an index: {path}")
    return rows, mapping


def read_certificates(path: Path, field: str) -> dict[int, list[dict]]:
    output = {}
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            entry = json.loads(line)
            index = int(entry["index"])
            if index in output:
                raise ValueError(f"certificate archive repeats graph {index}")
            output[index] = entry[field]
    return output


def load_inputs(
    input_path: Path,
    selection_path: Path,
    report_path: Path,
    expected_report_sha256: str,
) -> tuple[
    list[dict], dict[int, dict], dict[int, list[dict]],
    dict[int, dict], dict[int, list[dict]], dict,
]:
    if sha256(input_path) != EXPECTED_INPUT_SHA256:
        raise ValueError("rank input hash mismatch")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if len(payload.get("graphs", [])) != EXPECTED_INPUT_GRAPHS:
        raise ValueError("rank input count mismatch")
    graph_by_index = {
        int(graph["index"]): graph for graph in payload["graphs"]
    }
    if len(graph_by_index) != EXPECTED_INPUT_GRAPHS:
        raise ValueError("rank input repeats an index")
    if sha256(selection_path) != EXPECTED_SELECTION_SHA256:
        raise ValueError("full tetrad selection hash mismatch")
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    selected_indices = [int(value) for value in selection["selected_indices"]]
    if (
        len(selected_indices) != EXPECTED_SELECTED_GRAPHS
        or selection["selected_indices_sha256"] != stable_hash(selected_indices)
    ):
        raise ValueError("full tetrad selection index mismatch")
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
        raise ValueError("prior degree-one archive hash mismatch")
    _prior_rows, prior_by_index = read_decisions(prior_decision_path)
    prior_certificates = read_certificates(
        prior_certificate_path, "dual_failure_witnesses"
    )

    report_hash = sha256(report_path)
    if report_hash != expected_report_sha256:
        raise ValueError(
            f"tetrad report hash {report_hash} != explicit pin "
            f"{expected_report_sha256}"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    run_graphs = int(report["summary"]["graphs"])
    if (
        report.get("schema") != 1
        or not 0 < run_graphs <= EXPECTED_SELECTED_GRAPHS
        or report["summary"]["complete"] != run_graphs
        or report["summary"]["infra_errors"] != 0
    ):
        raise ValueError("tetrad report is incomplete")
    decision_path = resolve_relative(
        report_path, report["artifacts"]["decisions_archive"]
    )
    certificate_path = resolve_relative(
        report_path, report["artifacts"]["certificate_archive"]
    )
    if (
        sha256(decision_path)
        != report["artifacts"]["decisions_archive_sha256"]
        or sha256(certificate_path)
        != report["artifacts"]["certificate_archive_sha256"]
    ):
        raise ValueError("tetrad output archive hash mismatch")
    decision_rows, decision_by_index = read_decisions(decision_path)
    run_indices = selected_indices[:run_graphs]
    if [int(row["index"]) for row in decision_rows] != run_indices:
        raise ValueError("tetrad decision order differs from selection")
    status_counts = Counter(row["status"] for row in decision_rows)
    rejected_indices = [
        int(row["index"]) for row in decision_rows
        if row["status"] == "REJECTED"
    ]
    if (
        dict(status_counts) != report["summary"]["status_counts"]
        or rejected_indices
        != report["summary"]["marginal_tetrad_rejected_indices"]
        or len(rejected_indices)
        != report["summary"]["marginal_tetrad_rejected"]
    ):
        raise ValueError("tetrad report summary differs from decisions")
    tetrad_certificates = read_certificates(
        certificate_path, "tetrad_failure_witnesses"
    )
    selected = []
    for ordinal, index in enumerate(run_indices):
        if index not in graph_by_index:
            raise ValueError("selected graph absent from rank input")
        if prior_by_index[index]["status"] != "SURVIVOR":
            raise ValueError("selected graph was not a prior survivor")
        selected.append({**graph_by_index[index], "ordinal": ordinal})
    provenance = {
        "input": str(input_path), "input_sha256": sha256(input_path),
        "selection": str(selection_path),
        "selection_sha256": sha256(selection_path),
        "prior_decisions": str(prior_decision_path),
        "prior_decisions_sha256": sha256(prior_decision_path),
        "prior_certificates": str(prior_certificate_path),
        "prior_certificates_sha256": sha256(prior_certificate_path),
        "tetrad_report": str(report_path),
        "tetrad_report_sha256": report_hash,
        "tetrad_decisions": str(decision_path),
        "tetrad_decisions_sha256": sha256(decision_path),
        "tetrad_certificates": str(certificate_path),
        "tetrad_certificates_sha256": sha256(certificate_path),
    }
    return (
        selected, prior_by_index, prior_certificates,
        decision_by_index, tetrad_certificates, provenance,
    )


def _parse_bool(value: str) -> bool:
    if value in ("True", "1", "true"):
        return True
    if value in ("False", "0", "false"):
        return False
    raise ValueError(f"invalid serialized boolean {value!r}")


def verify_graph(payload: tuple[
    dict, dict, list[dict], dict, list[dict]
]) -> dict:
    graph, prior_row, prior_witnesses, decision, tetrad_witnesses = payload
    started = time.perf_counter()
    adj = tuple(graph["adjacency"])
    prior.validate_graph(adj)
    prior_by_key = {
        (tuple(item["seed"]), int(item["zmask"])): item
        for item in prior_witnesses
    }
    tetrad_by_key = {
        (tuple(item["seed"]), int(item["zmask"])): item
        for item in tetrad_witnesses
    }
    if (
        len(prior_by_key) != len(prior_witnesses)
        or len(tetrad_by_key) != len(tetrad_witnesses)
    ):
        raise ValueError("duplicate cover witness")
    consumed_prior = set()
    consumed_tetrad = set()
    counts = Counter()
    strict_cliques = 0
    prior_rejected = tetrad_rejected = False
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()
    for seed_mask in prior.clique_masks(adj, 7):
        counts["seeds"] += 1
        seed, outside, defects, ladj, eligible = prior.seed_instance(
            adj, seed_mask
        )
        total_term_rank = prior.matching_size(defects)
        seed_prior_passes = seed_tetrad_passes = 0
        covers = prior.eligible_covers(ladj, eligible)
        counts["covers"] += len(covers)
        for zmask in covers:
            analysis = prior.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            counts["enhanced_passing_covers"] += 1
            nvertices = [
                outside[index] for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            graph_n = tuple(prior.induced_graph(adj, nvertices))
            saturating = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            if any(
                strict_h.assess_saturating_clique(graph_n, mask).failed
                for mask in saturating
            ):
                continue
            strict_cliques += len(saturating)
            counts["strict_h_passing_covers"] += 1
            key = (tuple(seed), zmask)
            old_witness = prior_by_key.get(key)
            if old_witness is not None:
                if old_witness["nvertices"] != nvertices:
                    raise ValueError("prior witness induced vertices mismatch")
                certificate = old_witness["certificate"]
                clique_mask = sum(
                    1 << vertex for vertex in certificate["clique"]
                )
                if clique_mask not in saturating:
                    raise ValueError("prior certificate clique is not saturating")
                prior_checker.verify_clique_certificate(
                    graph_n, clique_mask, certificate
                )
                consumed_prior.add(key)
                counts["prior_dual_failed_covers"] += 1
                continue
            counts["prior_dual_passing_covers"] += 1
            seed_prior_passes += 1
            if saturating:
                seed_tetrad_passes += 1
                counts["tetrad_passing_covers"] += 1
                continue
            counts["prior_passing_no_saturating_clique_covers"] += 1
            near_size = analysis.k_rank_upper - 1
            near_cliques = tuple(prior.clique_masks(
                graph_n, near_size
            )) if near_size > 0 else ()
            if not near_cliques:
                counts["covers_without_near_clique"] += 1
                seed_tetrad_passes += 1
                counts["tetrad_passing_covers"] += 1
                continue
            counts["covers_with_near_clique"] += 1
            witness = tetrad_by_key.get(key)
            if witness is None:
                counts["near_cliques_tested"] += len(near_cliques)
                seed_tetrad_passes += 1
                counts["tetrad_passing_covers"] += 1
                continue
            if (
                witness["nvertices"] != nvertices
                or int(witness["k_rank_upper"]) != analysis.k_rank_upper
            ):
                raise ValueError("tetrad witness cover context mismatch")
            certificate = witness["certificate"]
            if certificate.get("output_degree") != 4:
                raise ValueError("production tetrad certificate is not degree four")
            clique_mask = sum(
                1 << vertex for vertex in certificate["clique"]
            )
            try:
                position = near_cliques.index(clique_mask)
            except ValueError as error:
                raise ValueError(
                    "tetrad certificate clique is not near-saturating"
                ) from error
            counts["near_cliques_tested"] += position + 1
            tetrad_checker.verify_clique_certificate(
                graph_n, clique_mask, certificate
            )
            consumed_tetrad.add(key)
            counts["tetrad_failed_covers"] += 1
        if not seed_prior_passes:
            prior_rejected = True
        if not seed_tetrad_passes:
            tetrad_rejected = True
    if consumed_prior != set(prior_by_key):
        raise ValueError("orphan prior certificate")
    if consumed_tetrad != set(tetrad_by_key):
        raise ValueError("orphan tetrad certificate")
    marginal = tetrad_rejected and not prior_rejected
    expected_status = "REJECTED" if marginal else "SURVIVOR"
    if decision["status"] != expected_status:
        raise ValueError("production status differs from independent replay")
    expected_fields = {
        **{name: counts[name] for name in (
            "seeds", "covers", "enhanced_passing_covers",
            "strict_h_passing_covers", "prior_dual_failed_covers",
            "prior_dual_passing_covers",
            "prior_passing_no_saturating_clique_covers",
            "covers_with_near_clique", "covers_without_near_clique",
            "near_cliques_tested", "tetrad_failed_covers",
            "tetrad_passing_covers",
        )},
        "prior_dual_rejected": prior_rejected,
        "tetrad_rejected": tetrad_rejected,
        "marginal_tetrad_rejected": marginal,
    }
    for name, expected in expected_fields.items():
        if isinstance(expected, bool):
            observed = _parse_bool(decision[name])
        else:
            observed = int(decision[name])
        if observed != expected:
            raise ValueError(
                f"decision field {name}: {observed!r} != {expected!r}"
            )
    prior_expected = {
        "seeds": counts["seeds"],
        "covers": counts["covers"],
        "enhanced_passing_covers": counts["enhanced_passing_covers"],
        "strict_h_passing_covers": counts["strict_h_passing_covers"],
        "dual_passing_covers": counts["prior_dual_passing_covers"],
        "strict_h_passing_cliques": strict_cliques,
        "dual_failing_cliques": counts["prior_dual_failed_covers"],
    }
    for name, expected in prior_expected.items():
        if int(prior_row[name]) != expected:
            raise ValueError(f"prior decision field {name} mismatch")
    return {
        "schema": 1,
        "ordinal": int(graph["ordinal"]),
        "index": int(graph["index"]),
        "status": "PASS",
        "production_status": expected_status,
        "prior_certificates_checked": len(consumed_prior),
        "tetrad_certificates_checked": len(consumed_tetrad),
        "elapsed_seconds": time.perf_counter() - started,
    }


def evaluate(payload: tuple) -> dict:
    graph = payload[0]
    started = time.perf_counter()
    try:
        return verify_graph(payload)
    except BaseException as error:  # noqa: BLE001 - verification record
        return {
            "schema": 1,
            "ordinal": int(graph["ordinal"]),
            "index": int(graph["index"]),
            "status": "ERROR",
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
    output = {}
    for path in sorted(directory.glob("graph_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("config_sha256") != config_sha256:
            raise ValueError(f"verification checkpoint config mismatch: {path}")
        record = payload["record"]
        ordinal = int(record["ordinal"])
        if expected.get(ordinal) != int(record["index"]):
            raise ValueError(f"verification checkpoint graph mismatch: {path}")
        if ordinal in output:
            raise ValueError("duplicate verification checkpoint")
        output[ordinal] = record
    return output


def git_provenance() -> dict:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=ROOT, check=True, capture_output=True, text=True,
        ).stdout
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}
    return {
        "available": True, "commit": commit, "dirty": bool(status),
        "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    parser.add_argument(
        "--selection", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_selection.json"),
    )
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_report.json"),
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--max-inflight", type=int, default=22)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=Path(".runs/d6_k7_rankone_tetrad_full_verification"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_verification_report.json"),
    )
    parser.add_argument("--progress-every", type=int, default=64)
    args = parser.parse_args()
    dependencies = verify_dependencies()
    (
        selected, prior_rows, prior_certificates,
        decisions, tetrad_certificates, provenance,
    ) = load_inputs(
        args.input, args.selection, args.report, args.report_sha256
    )
    if args.limit is not None:
        if not 0 <= args.limit <= len(selected):
            raise ValueError("invalid limit")
        selected = selected[:args.limit]
    config = {
        "schema": 1,
        "provenance": provenance,
        "indices": [int(graph["index"]) for graph in selected],
        "dependencies_sha256": dependencies,
    }
    config_sha256 = stable_hash(config)
    with exclusive_lock(args.checkpoint_dir):
        paths = list(args.checkpoint_dir.glob("graph_*.json"))
        if paths and not args.resume:
            raise ValueError(
                "verification checkpoint directory is nonempty; pass --resume"
            )
        completed = load_checkpoints(
            args.checkpoint_dir, config_sha256, selected
        ) if args.resume else {}
        pending = deque(
            graph for graph in selected
            if int(graph["ordinal"]) not in completed
        )
        initial_completed = len(completed)
        started = time.perf_counter()
        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=initialize_worker
        ) as executor:
            inflight = {}
            while pending or inflight:
                while pending and len(inflight) < args.max_inflight:
                    graph = pending.popleft()
                    index = int(graph["index"])
                    future = executor.submit(evaluate, (
                        graph,
                        prior_rows[index],
                        prior_certificates.get(index, []),
                        decisions[index],
                        tetrad_certificates.get(index, []),
                    ))
                    inflight[future] = graph
                done, _ = wait(inflight, return_when=FIRST_COMPLETED)
                for future in done:
                    graph = inflight.pop(future)
                    record = future.result()
                    ordinal = int(graph["ordinal"])
                    atomic_json(
                        checkpoint_path(
                            args.checkpoint_dir, ordinal, int(graph["index"])
                        ),
                        {
                            "schema": 1,
                            "config_sha256": config_sha256,
                            "record": record,
                        },
                    )
                    completed[ordinal] = record
                    if (
                        len(completed) % args.progress_every == 0
                        or len(completed) == len(selected)
                    ):
                        print(
                            f"[{utc_now()}] {len(completed)}/{len(selected)} "
                            f"elapsed={time.perf_counter()-started:.1f}s",
                            flush=True,
                        )
        records = [completed[ordinal] for ordinal in range(len(selected))]
        statuses = Counter(record["status"] for record in records)
        report = {
            "schema": 1,
            "description": (
                "Independent complete quantifier and rational-identity replay "
                "for the full K7 rank-one tetrad campaign."
            ),
            "config_sha256": config_sha256,
            "provenance": provenance,
            "dependencies_sha256": dependencies,
            "verifier_source_sha256": sha256(Path(__file__)),
            "git": git_provenance(),
            "environment": {
                "python": sys.version,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "workers": args.workers,
            },
            "runtime": {
                "finished": utc_now(),
                "wall_seconds": time.perf_counter() - started,
                "resumed_records": initial_completed,
                "computed_records": len(records) - initial_completed,
            },
            "summary": {
                "graphs": len(records),
                "status_counts": dict(statuses),
                "prior_certificates_checked": sum(
                    record.get("prior_certificates_checked", 0)
                    for record in records
                ),
                "tetrad_certificates_checked": sum(
                    record.get("tetrad_certificates_checked", 0)
                    for record in records
                ),
                "verified_rejections": sum(
                    record.get("production_status") == "REJECTED"
                    for record in records if record["status"] == "PASS"
                ),
            },
            "errors": [
                record for record in records if record["status"] != "PASS"
            ],
        }
        atomic_json(args.output, report)
        if statuses.get("ERROR", 0):
            raise SystemExit(2)


if __name__ == "__main__":
    main()
