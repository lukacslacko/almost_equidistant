#!/usr/bin/env python3
"""Restartable short-circuit runner for the exact strict-H K7 layer.

The mathematical kernels are frozen in ``d6_k7_rank_reference.py`` and
``d6_k7_special_h_reference.py``.  This file changes only quantifier and
campaign management:

* the input is restricted to graphs surviving the enhanced C rank pass;
* covers of size at most three are tried in deterministic order;
* a seed stops at its first strict-H-passing cover;
* a graph stops at its first seed with no passing cover;
* fixed ordered chunks are written atomically and can be resumed;
* the final report and TSV are reconstructed in original input order.

Every reported rejection retains one exact strict-H witness for every
enhanced-rank-passing cover of its rejecting seed.  Candidate nonedges keep
their optional-zero semantics through the imported frozen references.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import resource
import shlex
import sys
import time
import traceback
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


# Keep numerical discovery libraries single-threaded inside each exact worker.
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    # This runner owns its worker environment.  Inheriting, for example,
    # OPENBLAS_NUM_THREADS=8 would silently turn twelve exact workers into
    # ninety-six numerical threads when SciPy is imported by the frozen
    # strict-H reference.
    os.environ[_name] = "1"


ROOT = Path(__file__).resolve().parent
SPECIAL_SOURCE = ROOT / "d6_k7_special_h_reference.py"
RANK_SOURCE = ROOT / "d6_k7_rank_reference.py"
SAMPLE = ROOT / "d6_k7_rank_sample.json"
SAMPLE_RANK_REPORT = ROOT / "d6_k7_rank_report.json"
SAMPLE_STRICT_REPORT = ROOT / "d6_k7_strict_h_report.json"
FULL_INPUT = Path("/tmp/d6_k7_rank_full_residue.txt")
FULL_DECISIONS = Path("/tmp/d6_k7_rank_final_decisions.tsv")
FULL_C_REPORT = ROOT / "d6_k7_rank_c_full_report.json"

EXPECTED_SPECIAL_SHA256 = (
    "e649633d7262941c9c4608ea694ae0a5d6f29c146d47060c3dae71a036a9d190"
)
EXPECTED_RANK_SHA256 = (
    "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
)
EXPECTED_SAMPLE_SHA256 = (
    "bf699f36a1c90b6d6751498e52fab73e6b58db322c5b71b04accf32623c2c0dd"
)
EXPECTED_SAMPLE_RANK_REPORT_SHA256 = (
    "60af52ea66b37eb5beb47f68fa53ca2172612498feefa26e74e9854ec12d5a6f"
)
EXPECTED_SAMPLE_STRICT_REPORT_SHA256 = (
    "d2ff29f493d9997887a392d0e56f228b2395e48a324f3118dcbeeff9949f266f"
)
EXPECTED_FULL_INPUT_SHA256 = (
    "a1934445c9e3fa1b6e28168fbefe161fa6a5f33ea75e74a0aaf71553eceda16b"
)
EXPECTED_FULL_DECISIONS_SHA256 = (
    "eaa4d5061b44cc86e6ab98ac6315b29a509539e5ac22ab33d1b7a03a34358391"
)
EXPECTED_FULL_C_REPORT_SHA256 = (
    "ac3e73e728ca95624d0bf1da7b62777239289e052d08c62b1b863f4e0f3aa971"
)
EXPECTED_SAMPLE_GRAPHS = 86
EXPECTED_SAMPLE_REJECTIONS = (2215178, 3331579)
EXPECTED_FULL_GRAPHS = 17_764


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> str:
    actual = file_sha256(path)
    if actual != expected:
        raise RuntimeError(
            f"frozen input changed: {path}: {actual}, expected {expected}"
        )
    return actual


# Refuse to import mutable mathematical code before verifying both sources.
require_hash(SPECIAL_SOURCE, EXPECTED_SPECIAL_SHA256)
require_hash(RANK_SOURCE, EXPECTED_RANK_SHA256)
import d6_k7_rank_reference as rank  # noqa: E402
import d6_k7_special_h_reference as strict_h  # noqa: E402


@dataclass(frozen=True)
class InputGraph:
    ordinal: int
    index: int
    adjacency: tuple[int, ...]


@dataclass
class StrictGraphResult:
    ordinal: int
    index: int
    decision: str
    seeds_examined: int
    covers_examined: int
    prior_passing_covers_examined: int
    strict_h_cliques_examined: int
    prior_failure_counts: dict[str, int] = field(default_factory=dict)
    strict_h_reason_counts: dict[str, int] = field(default_factory=dict)
    seed_pass_witnesses: list[dict] = field(default_factory=list)
    rejecting_seed: list[int] | None = None
    rejecting_cover_witnesses: list[dict] = field(default_factory=list)
    elapsed_seconds: float = 0.0


def _prior_failure_kind(analysis: rank.CoverAnalysis) -> str:
    """Return a deterministic first prior rule for effort diagnostics only."""

    if analysis.direct_cap_failure is not None:
        return analysis.direct_cap_failure
    checks = (
        ("support", analysis.support_failed),
        ("subspace_K", analysis.subspace_k_failed),
        ("component_B", analysis.component_b_failed),
        ("PD_clique", analysis.pd_clique_failed),
        ("perpendicular_degree", analysis.perpendicular_degree_failed),
        ("basis_kernel", analysis.basis_kernel_failed),
        ("saturating_masks", analysis.saturating_mask_failed),
    )
    return next((name for name, failed in checks if failed), "unknown")


def _strict_cover_assessment(
    adj: Sequence[int],
    outside: Sequence[int],
    zmask: int,
    analysis: rank.CoverAnalysis,
) -> tuple[bool, int, Counter[str], dict | None]:
    """Stop at the first failing saturating clique of one prior-pass cover."""

    # The frozen exhaustive reference treats rank upper bound zero as an
    # empty family of nonempty saturating cliques.  ``clique_masks(G, 0)``
    # itself yields the empty mask, which is not a legal input to the H
    # checker, so preserve that convention explicitly here.
    if analysis.k_rank_upper == 0:
        return True, 0, Counter(), None

    nvertices = [
        outside[local]
        for local in range(len(outside))
        if not (zmask & (1 << local))
    ]
    graph_n = rank.induced_graph(adj, nvertices)
    checked = 0
    reasons: Counter[str] = Counter()
    for clique_mask in rank.clique_masks(graph_n, analysis.k_rank_upper):
        checked += 1
        result = strict_h.assess_saturating_clique(graph_n, clique_mask)
        reasons[result.reason] += 1
        if result.failed:
            witness = {
                "zmask": zmask,
                "N": nvertices,
                "k_rank_upper": analysis.k_rank_upper,
                "local_clique": list(rank.bits(clique_mask)),
                "reason": result.reason,
                "exact_witness": result.witness,
            }
            return False, checked, reasons, witness
    return True, checked, reasons, None


def analyze_graph_short_circuit(graph: InputGraph) -> StrictGraphResult:
    """Apply the nested existential cover quantifiers with safe early exits."""

    started = time.perf_counter()
    adj = graph.adjacency
    rank.validate_graph(adj)
    support_solver = rank.SupportSolver()
    zero_forcing = rank.ZeroForcingSolver()
    clique_solver = rank.CliqueStructureSolver()
    seeds_examined = covers_examined = prior_passing = strict_cliques = 0
    prior_reasons: Counter[str] = Counter()
    h_reasons: Counter[str] = Counter()
    seed_passes: list[dict] = []

    for seed_mask in rank.clique_masks(adj, 7):
        seeds_examined += 1
        seed_vertices, outside, defects, ladj, eligible = rank.seed_instance(
            adj, seed_mask
        )
        # The separately proved direct transitions force |Z|<=3 at n=19.
        covers = rank.eligible_covers(ladj, eligible, cap=3)
        total_term_rank = rank.matching_size(defects)
        seed_prior_passes = 0
        failed_h_covers: list[dict] = []
        seed_pass = None
        for zmask in covers:
            covers_examined += 1
            analysis = rank.analyze_cover(
                adj,
                outside,
                defects,
                zmask,
                support_solver,
                zero_forcing,
                total_term_rank,
                clique_solver,
            )
            if analysis.enhanced_joint_failed:
                prior_reasons[_prior_failure_kind(analysis)] += 1
                continue
            seed_prior_passes += 1
            prior_passing += 1
            passed, checked, reasons, witness = _strict_cover_assessment(
                adj, outside, zmask, analysis
            )
            strict_cliques += checked
            h_reasons.update(reasons)
            if passed:
                seed_pass = {
                    "seed": seed_vertices,
                    "zmask": zmask,
                    "k_rank_upper": analysis.k_rank_upper,
                    "strict_h_cliques_checked": checked,
                }
                seed_passes.append(seed_pass)
                break
            assert witness is not None
            failed_h_covers.append(witness)

        if seed_prior_passes == 0:
            raise RuntimeError(
                f"enhanced survivor {graph.index} has prior-rejected seed "
                f"{seed_vertices}; frozen C/Python decisions disagree"
            )
        if seed_pass is None:
            if len(failed_h_covers) != seed_prior_passes:
                raise AssertionError("strict-H cover witnesses are incomplete")
            return StrictGraphResult(
                ordinal=graph.ordinal,
                index=graph.index,
                decision="REJECTED",
                seeds_examined=seeds_examined,
                covers_examined=covers_examined,
                prior_passing_covers_examined=prior_passing,
                strict_h_cliques_examined=strict_cliques,
                prior_failure_counts=dict(prior_reasons),
                strict_h_reason_counts=dict(h_reasons),
                seed_pass_witnesses=seed_passes,
                rejecting_seed=seed_vertices,
                rejecting_cover_witnesses=failed_h_covers,
                elapsed_seconds=time.perf_counter() - started,
            )

    if not seeds_examined:
        raise RuntimeError(f"selected graph {graph.index} has no K7")
    return StrictGraphResult(
        ordinal=graph.ordinal,
        index=graph.index,
        decision="SURVIVOR",
        seeds_examined=seeds_examined,
        covers_examined=covers_examined,
        prior_passing_covers_examined=prior_passing,
        strict_h_cliques_examined=strict_cliques,
        prior_failure_counts=dict(prior_reasons),
        strict_h_reason_counts=dict(h_reasons),
        seed_pass_witnesses=seed_passes,
        elapsed_seconds=time.perf_counter() - started,
    )


def load_sample_graphs() -> tuple[list[InputGraph], dict[str, str]]:
    hashes = {
        "sample": require_hash(SAMPLE, EXPECTED_SAMPLE_SHA256),
        "rank_report": require_hash(
            SAMPLE_RANK_REPORT, EXPECTED_SAMPLE_RANK_REPORT_SHA256
        ),
        "strict_report": require_hash(
            SAMPLE_STRICT_REPORT, EXPECTED_SAMPLE_STRICT_REPORT_SHA256
        ),
    }
    sample = json.loads(SAMPLE.read_text(encoding="utf-8"))
    rank_report = json.loads(SAMPLE_RANK_REPORT.read_text(encoding="utf-8"))
    selected = set(
        rank_report["individual_graph_decisions"]["enhanced_survivors"]
    )
    graphs = [
        InputGraph(
            ordinal=ordinal,
            index=graph["index"],
            adjacency=tuple(graph["adjacency"]),
        )
        for ordinal, graph in enumerate(
            graph for graph in sample["graphs"] if graph["index"] in selected
        )
    ]
    if len(graphs) != EXPECTED_SAMPLE_GRAPHS:
        raise RuntimeError(f"sample selected {len(graphs)}, expected 86")
    return graphs, hashes


def load_full_graphs() -> tuple[list[InputGraph], dict[str, str]]:
    hashes = {
        "rank_input": require_hash(FULL_INPUT, EXPECTED_FULL_INPUT_SHA256),
        "rank_decisions": require_hash(
            FULL_DECISIONS, EXPECTED_FULL_DECISIONS_SHA256
        ),
        "rank_report": require_hash(
            FULL_C_REPORT, EXPECTED_FULL_C_REPORT_SHA256
        ),
    }
    graphs: list[InputGraph] = []
    with FULL_INPUT.open(encoding="ascii") as input_stream, (
        FULL_DECISIONS.open(newline="", encoding="ascii")
    ) as decision_stream:
        decisions = csv.DictReader(decision_stream, delimiter="\t")
        required_columns = {"index", "internal_error", "joint_rejected"}
        if not required_columns.issubset(decisions.fieldnames or ()):
            raise ValueError("full rank decisions have an unexpected header")
        for input_line, decision in zip(input_stream, decisions, strict=True):
            fields = input_line.split()
            if len(fields) != 21 or fields[1] != "19":
                raise ValueError("malformed full rank-residue input")
            index = int(fields[0])
            if index != int(decision["index"]):
                raise ValueError("rank input and decisions lost alignment")
            if int(decision["internal_error"]):
                raise ValueError(f"rank evaluator internal error at {index}")
            if not int(decision["joint_rejected"]):
                graphs.append(
                    InputGraph(
                        ordinal=len(graphs),
                        index=index,
                        adjacency=tuple(map(int, fields[2:])),
                    )
                )
    if len(graphs) != EXPECTED_FULL_GRAPHS:
        raise RuntimeError(
            f"full enhanced survivors {len(graphs)}, expected 17764"
        )
    return graphs, hashes


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(value, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        sink.write(value)
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def chunks(
    graphs: Sequence[InputGraph], size: int
) -> list[tuple[int, int, list[InputGraph]]]:
    return [
        (start, min(start + size, len(graphs)), list(graphs[start:start + size]))
        for start in range(0, len(graphs), size)
    ]


def chunk_path(directory: Path, start: int, stop: int) -> Path:
    return directory / f"chunk_{start:06d}_{stop:06d}.json"


def process_chunk(payload: tuple[str, int, int, list[InputGraph]]) -> dict:
    config_hash, start, stop, graphs = payload
    started = time.perf_counter()
    results = [analyze_graph_short_circuit(graph) for graph in graphs]
    maximum_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "schema": 1,
        "config_sha256": config_hash,
        "start": start,
        "stop": stop,
        "results": [asdict(result) for result in results],
        "worker_elapsed_seconds": time.perf_counter() - started,
        "worker_max_rss_native_units": maximum_rss,
    }


def validate_chunk(
    value: dict,
    config_hash: str,
    start: int,
    stop: int,
    graphs: Sequence[InputGraph],
) -> list[dict]:
    if value.get("schema") != 1 or value.get("config_sha256") != config_hash:
        raise ValueError("checkpoint schema/config mismatch")
    if value.get("start") != start or value.get("stop") != stop:
        raise ValueError("checkpoint range mismatch")
    results = value.get("results")
    if not isinstance(results, list) or len(results) != stop - start:
        raise ValueError("checkpoint result count mismatch")
    for expected, result in zip(graphs, results, strict=True):
        if (
            result.get("ordinal") != expected.ordinal
            or result.get("index") != expected.index
            or result.get("decision") not in ("REJECTED", "SURVIVOR")
        ):
            raise ValueError("checkpoint result identity/decision mismatch")
    return results


def progress_summary(results: Iterable[dict], total: int) -> dict:
    materialized = list(results)
    counts = Counter(result["decision"] for result in materialized)
    return {
        "completed_graphs": len(materialized),
        "total_graphs": total,
        "rejected": counts["REJECTED"],
        "survivors": counts["SURVIVOR"],
    }


def run_campaign(
    graphs: list[InputGraph],
    *,
    mode: str,
    input_hashes: dict[str, str],
    workers: int,
    chunk_size: int,
    checkpoint_dir: Path,
    report_path: Path,
    decisions_path: Path,
) -> dict:
    if workers < 1 or chunk_size < 1:
        raise ValueError("workers and chunk size must be positive")
    configuration = {
        "schema": 1,
        "mode": mode,
        "graphs": len(graphs),
        "indices_sha256": stable_hash([graph.index for graph in graphs]),
        "chunk_size": chunk_size,
        "runner_source_sha256": file_sha256(Path(__file__).resolve()),
        "special_H_source_sha256": EXPECTED_SPECIAL_SHA256,
        "rank_source_sha256": EXPECTED_RANK_SHA256,
        "input_hashes": input_hashes,
    }
    config_hash = stable_hash(configuration)
    campaign_path = checkpoint_dir / "campaign.json"
    if campaign_path.exists():
        existing = json.loads(campaign_path.read_text(encoding="utf-8"))
        if existing.get("configuration") != configuration:
            raise RuntimeError(
                f"checkpoint directory has a different campaign: {checkpoint_dir}"
            )
    else:
        atomic_json(
            campaign_path,
            {
                "schema": 1,
                "configuration": configuration,
                "config_sha256": config_hash,
                "created_utc": utc_now(),
                "argv": [sys.executable, *sys.argv],
                "command": shlex.join([sys.executable, *sys.argv]),
            },
        )

    work = chunks(graphs, chunk_size)
    completed: dict[int, dict] = {}
    missing: deque[tuple[int, int, list[InputGraph]]] = deque()
    for start, stop, part in work:
        path = chunk_path(checkpoint_dir, start, stop)
        if path.exists():
            value = json.loads(path.read_text(encoding="utf-8"))
            validate_chunk(value, config_hash, start, stop, part)
            completed[start] = value
        else:
            missing.append((start, stop, part))

    started = time.perf_counter()
    last_print = 0

    def ordered_results() -> list[dict]:
        return [
            result
            for start in sorted(completed)
            for result in completed[start]["results"]
        ]

    def checkpoint_progress() -> None:
        nonlocal last_print
        results = ordered_results()
        summary = progress_summary(results, len(graphs))
        atomic_json(
            checkpoint_dir / "progress.json",
            {
                "schema": 1,
                "config_sha256": config_hash,
                **summary,
                "completed_chunks": len(completed),
                "total_chunks": len(work),
                "updated_utc": utc_now(),
            },
        )
        completed_graphs = summary["completed_graphs"]
        if completed_graphs == len(graphs) or completed_graphs - last_print >= max(
            chunk_size, len(graphs) // 100
        ):
            print(
                f"progress {completed_graphs}/{len(graphs)}; "
                f"rejected {summary['rejected']}; survivors {summary['survivors']}",
                flush=True,
            )
            last_print = completed_graphs

    checkpoint_progress()
    if missing:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            active = {}
            while missing or active:
                while missing and len(active) < 2 * workers:
                    start, stop, part = missing.popleft()
                    future = executor.submit(
                        process_chunk, (config_hash, start, stop, part)
                    )
                    active[future] = (start, stop, part)
                finished, _ = wait(active, return_when=FIRST_COMPLETED)
                for future in finished:
                    start, stop, part = active.pop(future)
                    try:
                        value = future.result()
                    except BaseException as error:
                        failure_name = (
                            f"failure_{start:06d}_{stop:06d}_"
                            f"{time.time_ns()}.json"
                        )
                        atomic_json(
                            checkpoint_dir / "failures" / failure_name,
                            {
                                "schema": 1,
                                "config_sha256": config_hash,
                                "start": start,
                                "stop": stop,
                                "indices": [graph.index for graph in part],
                                "exception_type": type(error).__name__,
                                "exception": str(error),
                                "traceback": traceback.format_exc(),
                                "failed_utc": utc_now(),
                            },
                        )
                        checkpoint_progress()
                        raise
                    validate_chunk(value, config_hash, start, stop, part)
                    atomic_json(
                        chunk_path(checkpoint_dir, start, stop), value
                    )
                    completed[start] = value
                checkpoint_progress()

    results = ordered_results()
    if len(results) != len(graphs):
        raise AssertionError("campaign completed without every ordered result")
    if [result["ordinal"] for result in results] != list(range(len(graphs))):
        raise AssertionError("final results are not in input order")
    rejected = [
        result["index"] for result in results
        if result["decision"] == "REJECTED"
    ]
    survivors = [
        result["index"] for result in results
        if result["decision"] == "SURVIVOR"
    ]
    if mode == "sample":
        if tuple(sorted(rejected)) != EXPECTED_SAMPLE_REJECTIONS:
            raise AssertionError(
                f"sample strict-H decisions {sorted(rejected)} differ from "
                f"{list(EXPECTED_SAMPLE_REJECTIONS)}"
            )
        if len(survivors) != EXPECTED_SAMPLE_GRAPHS - len(
            EXPECTED_SAMPLE_REJECTIONS
        ):
            raise AssertionError("sample survivor count mismatch")

    prior_reasons: Counter[str] = Counter()
    strict_reasons: Counter[str] = Counter()
    for result in results:
        prior_reasons.update(result["prior_failure_counts"])
        strict_reasons.update(result["strict_h_reason_counts"])
    decision_lines = [
        "ordinal\tindex\tdecision\tseeds_examined\tcovers_examined\t"
        "prior_passing_covers_examined\tstrict_h_cliques_examined\t"
        "elapsed_seconds\n"
    ]
    decision_lines.extend(
        f"{result['ordinal']}\t{result['index']}\t{result['decision']}\t"
        f"{result['seeds_examined']}\t{result['covers_examined']}\t"
        f"{result['prior_passing_covers_examined']}\t"
        f"{result['strict_h_cliques_examined']}\t"
        f"{result['elapsed_seconds']:.9f}\n"
        for result in results
    )
    atomic_text(decisions_path, "".join(decision_lines))
    chunk_hashes = {
        chunk_path(checkpoint_dir, start, stop).name: file_sha256(
            chunk_path(checkpoint_dir, start, stop)
        )
        for start, stop, _ in work
    }
    report = {
        "schema": 1,
        "description": (
            "Restartable ordered exact strict-H pass over enhanced K7 rank "
            "survivors; early exits preserve the existential quantifiers."
        ),
        "configuration": configuration,
        "config_sha256": config_hash,
        "workers": workers,
        "checkpoint_directory": str(checkpoint_dir),
        "checkpoint_chunk_sha256": chunk_hashes,
        "graphs": len(results),
        "counts": {
            "strict_H_rejected": len(rejected),
            "strict_H_survivors": len(survivors),
            "seeds_examined": sum(result["seeds_examined"] for result in results),
            "covers_examined": sum(result["covers_examined"] for result in results),
            "prior_passing_covers_examined": sum(
                result["prior_passing_covers_examined"] for result in results
            ),
            "strict_H_cliques_examined": sum(
                result["strict_h_cliques_examined"] for result in results
            ),
        },
        "decisions": {
            "strict_H_rejected": rejected,
            "strict_H_survivors": survivors,
        },
        "prior_failure_counts": dict(prior_reasons),
        "strict_H_reason_counts": dict(strict_reasons),
        "runtime": {
            "resume_wall_seconds": time.perf_counter() - started,
            "sum_graph_seconds": sum(result["elapsed_seconds"] for result in results),
            "maximum_graph_seconds": max(
                (result["elapsed_seconds"] for result in results), default=0.0
            ),
            "worker_max_rss_native_units": max(
                (
                    value.get("worker_max_rss_native_units", 0)
                    for value in completed.values()
                ),
                default=0,
            ),
        },
        "decision_TSV": {
            "file": decisions_path.name,
            "path": str(decisions_path),
            "sha256": file_sha256(decisions_path),
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "thread_limits": {
                name: os.environ.get(name)
                for name in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
        },
        "results": results,
    }
    atomic_json(report_path, report)
    print(
        f"wrote {report_path}: rejected {len(rejected)}, "
        f"survivors {len(survivors)}; decisions {decisions_path}",
        flush=True,
    )
    return report


def self_test() -> None:
    strict_h.self_test()
    require_hash(SPECIAL_SOURCE, EXPECTED_SPECIAL_SHA256)
    require_hash(RANK_SOURCE, EXPECTED_RANK_SHA256)
    sample, _ = load_sample_graphs()
    assert len(sample) == EXPECTED_SAMPLE_GRAPHS
    assert sample[0].ordinal == 0
    assert len({graph.index for graph in sample}) == len(sample)
    print("runner self-test passed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("sample", "full"), required=True)
    parser.add_argument("--workers", type=int, default=min(12, os.cpu_count() or 1))
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--decisions", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()

    if args.mode == "sample":
        graphs, input_hashes = load_sample_graphs()
        checkpoint_dir = args.checkpoint_dir or Path(
            "/tmp/d6_k7_strict_h_sample_checkpoints"
        )
        report = args.report or Path("/tmp/d6_k7_strict_h_short_sample.json")
        decisions = args.decisions or Path(
            "/tmp/d6_k7_strict_h_short_sample.tsv"
        )
    else:
        graphs, input_hashes = load_full_graphs()
        checkpoint_dir = args.checkpoint_dir or Path(
            "/tmp/d6_k7_strict_h_full_checkpoints"
        )
        report = args.report or Path("d6_k7_strict_h_full_report.json")
        decisions = args.decisions or Path(
            "/tmp/d6_k7_strict_h_full_decisions.tsv"
        )
    run_campaign(
        graphs,
        mode=args.mode,
        input_hashes=input_hashes,
        workers=args.workers,
        chunk_size=args.chunk_size,
        checkpoint_dir=checkpoint_dir,
        report_path=report,
        decisions_path=decisions,
    )


if __name__ == "__main__":
    main()
