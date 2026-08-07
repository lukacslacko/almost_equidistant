#!/usr/bin/env python3
"""Exact labeled-support propagation refinement for the K7 cover layer.

This module deliberately imports, but never modifies, the frozen independent
rank reference in :mod:`d6_k7_rank_reference`.  For a cap-three zero-factor
cover ``Z`` it enumerates the *labeled* actual supports of the vectors in
``Z``.  Orthogonality to those fixed vectors can then shrink the allowed
support masks of every nonzero-factor vertex.  The resulting smaller term
ranks sharpen all existing K/B rank and saturating-clique tests.

The command-line profiler is parallel, deterministic, and restartable.  Its
TSV is a decision log; a ``.partial`` suffix is retained after interruption
and accepted only with ``--resume``.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import d6_k7_rank_reference as reference


ROOT = Path(__file__).resolve().parent
REFERENCE_PATH = ROOT / "d6_k7_rank_reference.py"
EXPECTED_REFERENCE_SHA256 = (
    "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
)
QOS_CLASS_USER_INITIATED = 0x19
CACHE_LIMIT = 50_000
REASONS = (
    "empty_propagated_mask",
    "disjoint_required_edge",
    "clique",
    "degree",
    "basis",
    "mask",
    "subspace_K",
    "component_B",
)
DECISION_FIELDS = (
    "index",
    "applicable",
    "baseline_decision",
    "refined_decision",
    "first_baseline_failing_seed",
    "first_refined_failing_seed",
    "seeds",
    "covers",
    "baseline_passing_covers",
    "refined_passing_covers",
    "newly_failed_covers",
    "support_families_checked",
    "propagated_coordinate_deletions",
    "refinement_cache_hits",
) + tuple(f"assignment_fail_{reason}" for reason in REASONS)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def set_worker_qos() -> None:
    """Request non-background QoS for workers on macOS.

    Codex-launched processes can otherwise inherit background throttling.  A
    failed request is fatal in a campaign run but may be explicitly tolerated
    by controls through ``D6_ALLOW_BACKGROUND_TEST_ONLY=1``.
    """

    if sys.platform != "darwin":
        return
    libc = ctypes.CDLL(None)
    setter = libc.pthread_set_qos_class_self_np
    setter.argtypes = (ctypes.c_uint, ctypes.c_int)
    setter.restype = ctypes.c_int
    error = int(setter(QOS_CLASS_USER_INITIATED, 0))
    if error:
        if os.environ.get("D6_ALLOW_BACKGROUND_TEST_ONLY") == "1":
            print(
                f"WARNING: worker {os.getpid()} retained background QoS: "
                f"error {error}",
                file=sys.stderr,
                flush=True,
            )
            return
        raise OSError(error, "pthread_set_qos_class_self_np failed")


def labeled_support_families(
    allowed_masks: Sequence[int],
) -> Iterator[tuple[int, ...]]:
    """Yield every labeled family satisfying the frozen support rules.

    No exchangeability quotient is taken: output position ``j`` always
    belongs to the ``j``-th zero-factor vertex.  Covers used here have size at
    most three, but retaining the final general validity check makes this
    helper agree with the frozen rules for every family size up to seven.
    """

    allowed = tuple(int(mask) for mask in allowed_masks)
    if len(allowed) > 7:
        return
    domains = tuple(reference.SUPPORT_DOMAINS[mask] for mask in allowed)
    if any(not domain for domain in domains):
        return
    assignment = [0] * len(allowed)
    # Assign the smallest domain first while restoring the original labels in
    # the yielded tuple.  The tie-break makes enumeration deterministic.
    order = tuple(sorted(range(len(allowed)), key=lambda i: (len(domains[i]), i)))

    def visit(depth: int, selected: tuple[int, ...]) -> Iterator[tuple[int, ...]]:
        if depth == len(order):
            answer = tuple(assignment)
            if reference.support_family_valid(answer):
                yield answer
            return
        variable = order[depth]
        for support in domains[variable]:
            if any((support & prior).bit_count() == 1 for prior in selected):
                continue
            assignment[variable] = support
            chosen = selected + (support,)
            if reference.matching_size(chosen) != len(chosen):
                continue
            yield from visit(depth + 1, chosen)
        assignment[variable] = 0

    yield from visit(0, ())


def propagate_one_mask(mask: int, fixed_supports: Sequence[int]) -> tuple[int, int]:
    """Return the canonical singleton-overlap closure and deletion count."""

    current = int(mask)
    deleted = 0
    while True:
        forced_absent = 0
        for support in fixed_supports:
            intersection = current & support
            if intersection and not (intersection & (intersection - 1)):
                forced_absent |= intersection
        if not forced_absent:
            return current, deleted
        current &= ~forced_absent
        deleted += forced_absent.bit_count()


def propagate_n_masks(
    masks: Sequence[int], fixed_supports: Sequence[int]
) -> tuple[tuple[int, ...], int]:
    """Propagate every N mask independently against the fixed Z supports."""

    output = []
    deletions = 0
    for mask in masks:
        closed, count = propagate_one_mask(mask, fixed_supports)
        output.append(closed)
        deletions += count
    return tuple(output), deletions


@dataclass(frozen=True)
class AssignmentAnalysis:
    failure: str | None
    propagated_masks: tuple[int, ...]
    deletions: int
    n_term_rank: int
    total_term_rank: int
    k_rank_upper: int
    b_rank_upper: int


def analyze_support_assignment(
    graph_n: Sequence[int],
    z_supports: Sequence[int],
    n_allowed_masks: Sequence[int],
    zero_forcing: reference.ZeroForcingSolver,
    clique_solver: reference.CliqueStructureSolver,
) -> AssignmentAnalysis:
    """Apply propagation and every frozen exact rank/Schur screen."""

    if not reference.support_family_valid(z_supports):
        raise ValueError("invalid fixed zero-factor support family")
    if len(graph_n) != len(n_allowed_masks):
        raise ValueError("N graph and support-mask order disagree")
    propagated, deletions = propagate_n_masks(n_allowed_masks, z_supports)

    n_term_rank = reference.matching_size(propagated)
    total_term_rank = reference.matching_size((*z_supports, *propagated))
    k_upper, b_upper = reference.rank_upper_bounds(
        len(graph_n), n_term_rank, total_term_rank, len(z_supports)
    )

    def answer(failure: str | None) -> AssignmentAnalysis:
        return AssignmentAnalysis(
            failure=failure,
            propagated_masks=propagated,
            deletions=deletions,
            n_term_rank=n_term_rank,
            total_term_rank=total_term_rank,
            k_rank_upper=k_upper,
            b_rank_upper=b_upper,
        )

    if any(not mask for mask in propagated):
        return answer("empty_propagated_mask")
    for first, second in combinations(range(len(graph_n)), 2):
        if (
            graph_n[first] & (1 << second)
            and not (propagated[first] & propagated[second])
        ):
            return answer("disjoint_required_edge")

    structure = clique_solver.solve(graph_n, k_upper)
    if structure.clique_number > k_upper:
        return answer("clique")
    if graph_n and structure.f_maximum_degree > k_upper - 1:
        return answer("degree")
    if structure.basis_kernel_failures:
        return answer("basis")
    if structure.saturating_mask_failures:
        return answer("mask")

    k_zero_forcing = zero_forcing.solve(graph_n).number
    if len(graph_n) - k_zero_forcing > k_upper:
        return answer("subspace_K")
    fgraph = reference.complement_graph(graph_n)
    component_cap = reference.component_inertia_nullity_cap(
        fgraph, zero_forcing
    )
    if len(graph_n) - component_cap > b_upper:
        return answer("component_B")
    return answer(None)


@dataclass(frozen=True)
class CoverRefinement:
    passes: bool
    families_checked: int
    deletions: int
    failures: tuple[tuple[str, int], ...]
    first_failure: AssignmentAnalysis | None
    passing_assignment: tuple[int, ...] | None
    passing_analysis: AssignmentAnalysis | None


def refine_baseline_passing_cover(
    graph_n: Sequence[int],
    z_allowed_masks: Sequence[int],
    n_allowed_masks: Sequence[int],
    zero_forcing: reference.ZeroForcingSolver,
    clique_solver: reference.CliqueStructureSolver,
) -> CoverRefinement:
    """Existentially quantify labeled supports for one baseline-passing cover."""

    failures: Counter[str] = Counter()
    families = 0
    deletions = 0
    first_failure = None
    for supports in labeled_support_families(z_allowed_masks):
        families += 1
        analysis = analyze_support_assignment(
            graph_n,
            supports,
            n_allowed_masks,
            zero_forcing,
            clique_solver,
        )
        deletions += analysis.deletions
        if analysis.failure is None:
            return CoverRefinement(
                passes=True,
                families_checked=families,
                deletions=deletions,
                failures=tuple(sorted(failures.items())),
                first_failure=first_failure,
                passing_assignment=supports,
                passing_analysis=analysis,
            )
        failures[analysis.failure] += 1
        if first_failure is None:
            first_failure = analysis
    return CoverRefinement(
        passes=False,
        families_checked=families,
        deletions=deletions,
        failures=tuple(sorted(failures.items())),
        first_failure=first_failure,
        passing_assignment=None,
        passing_analysis=None,
    )


_support_solver: reference.SupportSolver | None = None
_zero_forcing: reference.ZeroForcingSolver | None = None
_clique_solver: reference.CliqueStructureSolver | None = None
_refinement_cache: dict[
    tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]], CoverRefinement
] | None = None


def initialize_worker() -> None:
    global _support_solver, _zero_forcing, _clique_solver, _refinement_cache
    set_worker_qos()
    _support_solver = reference.SupportSolver()
    _zero_forcing = reference.ZeroForcingSolver()
    _clique_solver = reference.CliqueStructureSolver()
    _refinement_cache = {}


def trim_worker_caches() -> None:
    assert _support_solver is not None
    assert _zero_forcing is not None
    assert _clique_solver is not None
    assert _refinement_cache is not None
    for cache in (
        _support_solver.cache,
        _zero_forcing.cache,
        _clique_solver.cache,
        _refinement_cache,
    ):
        if len(cache) > CACHE_LIMIT:
            cache.clear()


def evaluate_graph(graph: dict) -> tuple:
    """Return one deterministic TSV decision row for an input graph."""

    assert _support_solver is not None
    assert _zero_forcing is not None
    assert _clique_solver is not None
    assert _refinement_cache is not None

    adj = tuple(int(row) for row in graph["adjacency"])
    index = int(graph["index"])
    reference.validate_graph(adj)
    seed_masks = tuple(reference.clique_masks(adj, 7))
    if not seed_masks:
        return (
            index, 0, "NOT_APPLICABLE", "NOT_APPLICABLE", 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, *(0 for _ in REASONS),
        )

    covers_count = 0
    baseline_passing_covers = 0
    refined_passing_covers = 0
    newly_failed_covers = 0
    families_checked = 0
    deletions = 0
    cache_hits = 0
    failure_counts: Counter[str] = Counter()
    first_baseline_failing_seed = 0
    first_refined_failing_seed = 0

    for seed_mask in seed_masks:
        _, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        baseline_seed_passes = False
        refined_seed_passes = False
        for zmask in reference.eligible_covers(ladj, eligible, cap=3):
            covers_count += 1
            baseline = reference.analyze_cover(
                adj,
                outside,
                defects,
                zmask,
                _support_solver,
                _zero_forcing,
                total_term_rank,
                _clique_solver,
            )
            if baseline.enhanced_joint_failed:
                continue
            baseline_seed_passes = True
            baseline_passing_covers += 1

            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex for vertex in range(len(outside))
                if not (zmask & (1 << vertex))
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nvertices]
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            key = (graph_n, z_allowed, n_allowed)
            refinement = _refinement_cache.get(key)
            if refinement is None:
                refinement = refine_baseline_passing_cover(
                    graph_n,
                    z_allowed,
                    n_allowed,
                    _zero_forcing,
                    _clique_solver,
                )
                _refinement_cache[key] = refinement
            else:
                cache_hits += 1
            families_checked += refinement.families_checked
            deletions += refinement.deletions
            failure_counts.update(dict(refinement.failures))
            if refinement.passes:
                refined_seed_passes = True
                refined_passing_covers += 1
            else:
                newly_failed_covers += 1

        if not baseline_seed_passes and not first_baseline_failing_seed:
            first_baseline_failing_seed = seed_mask
        if not refined_seed_passes and not first_refined_failing_seed:
            first_refined_failing_seed = seed_mask

    trim_worker_caches()
    baseline_decision = (
        "REJECTED" if first_baseline_failing_seed else "SURVIVOR"
    )
    refined_decision = (
        "REJECTED" if first_refined_failing_seed else "SURVIVOR"
    )
    if baseline_decision == "REJECTED" and refined_decision != "REJECTED":
        raise AssertionError("refinement lost a baseline rejection")
    return (
        index,
        1,
        baseline_decision,
        refined_decision,
        first_baseline_failing_seed,
        first_refined_failing_seed,
        len(seed_masks),
        covers_count,
        baseline_passing_covers,
        refined_passing_covers,
        newly_failed_covers,
        families_checked,
        deletions,
        cache_hits,
        *(failure_counts[reason] for reason in REASONS),
    )


def atomic_decisions(
    path: Path, rows: Sequence[Sequence[object]]
) -> None:
    """Atomically replace a decision file with one validated row prefix."""

    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(DECISION_FIELDS)
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def load_completed(
    path: Path,
    expected: Sequence[int],
    committed_rows: int,
) -> list[tuple[str, ...]]:
    """Load exactly the manifest-committed prefix and discard a torn tail."""

    if committed_rows < 0 or committed_rows > len(expected):
        raise ValueError("checkpoint committed-row count is out of range")
    if not path.exists():
        if committed_rows:
            raise ValueError("checkpoint manifest names rows but TSV is absent")
        return []
    with path.open(newline="", encoding="ascii") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, None)
        if header != list(DECISION_FIELDS):
            raise ValueError(f"bad partial decision header in {path}")
        rows = [tuple(row) for row in reader]
    if len(rows) < committed_rows:
        raise ValueError(
            "partial decision file is shorter than its committed prefix"
        )
    completed = rows[:committed_rows]
    for position, row in enumerate(completed):
        if len(row) != len(DECISION_FIELDS):
            raise ValueError(
                f"committed decision row {position} has {len(row)} fields"
            )
        if int(row[0]) != expected[position]:
            raise ValueError(
                f"partial decision index {row[0]} at position {position}"
            )
        if row[2] not in {"REJECTED", "SURVIVOR", "NOT_APPLICABLE"}:
            raise ValueError(f"bad baseline decision in row {position}")
        if row[3] not in {"REJECTED", "SURVIVOR", "NOT_APPLICABLE"}:
            raise ValueError(f"bad refined decision in row {position}")
        for field in (*row[:2], *row[4:]):
            int(field)
    if len(rows) != committed_rows:
        # The manifest is updated only after fsync.  Rows beyond its count may
        # be a complete-but-uncommitted buffer or a torn final write; neither
        # belongs to the durable checkpoint.
        atomic_decisions(path, completed)
    return completed


def batches(values: Sequence[dict], size: int) -> Iterator[Sequence[dict]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def summarize_rows(rows: Sequence[Sequence[object]]) -> dict:
    baseline_rejected = [
        int(row[0]) for row in rows if row[2] == "REJECTED"
    ]
    refined_rejected = [
        int(row[0]) for row in rows if row[3] == "REJECTED"
    ]
    baseline_set = set(baseline_rejected)
    refined_set = set(refined_rejected)
    if not baseline_set <= refined_set:
        raise AssertionError("refined rejection set is not monotone")
    numeric_positions = {
        name: DECISION_FIELDS.index(name)
        for name in DECISION_FIELDS[6:]
    }
    return {
        "graphs": len(rows),
        "applicable_K7_graphs": sum(int(row[1]) for row in rows),
        "baseline_rejected": len(baseline_rejected),
        "baseline_survivors": len(rows) - len(baseline_rejected),
        "refined_rejected": len(refined_rejected),
        "refined_survivors": len(rows) - len(refined_rejected),
        "incremental_rejected": len(refined_set - baseline_set),
        "baseline_rejected_indices": baseline_rejected,
        "refined_rejected_indices": refined_rejected,
        "incremental_rejected_indices": sorted(refined_set - baseline_set),
        "totals": {
            name: sum(int(row[position]) for row in rows)
            for name, position in numeric_positions.items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k7_rank_sample.json"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
    )
    parser.add_argument(
        "--decisions", type=Path,
        default=Path("/tmp/d6_k7_support_propagation_decisions.tsv"),
    )
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_support_propagation_report.json",
    )
    parser.add_argument(
        "--pid-file", type=Path,
        default=Path("/tmp/d6_k7_support_propagation.pid.json"),
    )
    parser.add_argument("--oracle-report", type=Path)
    parser.add_argument(
        "--comparison-strict-h-report",
        type=Path,
        help=(
            "optional strict-H sample report; record set overlap only, "
            "without affecting support-propagation decisions"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--map-chunksize", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=32)
    args = parser.parse_args()

    reference_hash = sha256(REFERENCE_PATH)
    if reference_hash != EXPECTED_REFERENCE_SHA256:
        raise SystemExit(
            f"frozen reference hash is {reference_hash}, expected "
            f"{EXPECTED_REFERENCE_SHA256}"
        )
    input_hash = sha256(args.input)
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    graphs = payload["graphs"]
    indices = [int(graph["index"]) for graph in graphs]
    if len(indices) != len(set(indices)):
        raise SystemExit("input contains duplicate graph indices")

    source_hash = sha256(Path(__file__))
    command = [sys.executable, *sys.argv]
    compatibility = {
        "schema": 1,
        "input_sha256": input_hash,
        "source_sha256": source_hash,
        "reference_sha256": reference_hash,
        "graph_count": len(graphs),
        "decision_fields": list(DECISION_FIELDS),
    }
    partial = args.decisions.with_name(args.decisions.name + ".partial")
    checkpoint = args.decisions.with_name(args.decisions.name + ".checkpoint.json")
    checkpoint_previous = None
    if args.resume:
        if not checkpoint.exists():
            raise SystemExit(f"resume checkpoint is absent: {checkpoint}")
        checkpoint_previous = json.loads(
            checkpoint.read_text(encoding="utf-8")
        )
        if checkpoint_previous.get("compatibility") != compatibility:
            raise SystemExit(
                "resume checkpoint does not match the current input/source/"
                "reference/decision schema"
            )
        committed_rows = int(checkpoint_previous["committed_rows"])
        completed = load_completed(partial, indices, committed_rows)
    else:
        completed = []
        for stale in (partial, checkpoint):
            if stale.exists():
                stale.unlink()
        atomic_decisions(partial, ())

    initial_started_unix = (
        float(checkpoint_previous["initial_started_unix"])
        if checkpoint_previous is not None
        else time.time()
    )
    initial_command = (
        checkpoint_previous["initial_command"]
        if checkpoint_previous is not None
        else command
    )
    resume_count = (
        int(checkpoint_previous.get("resume_count", 0)) + 1
        if checkpoint_previous is not None
        else 0
    )

    def write_checkpoint(
        committed_rows: int,
        status: str,
        decisions_hash: str | None = None,
    ) -> None:
        value = {
            "schema": 1,
            "status": status,
            "compatibility": compatibility,
            "initial_started_unix": initial_started_unix,
            "initial_command": initial_command,
            "last_session_started_unix": session_started_unix,
            "last_command": command,
            "resume_count": resume_count,
            "committed_rows": committed_rows,
            "partial_decisions": str(partial),
            "final_decisions": str(args.decisions),
            "workers": args.workers,
            "batch_size": args.batch_size,
            "map_chunksize": args.map_chunksize,
            "progress_every": args.progress_every,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        }
        if decisions_hash is not None:
            value["decisions_sha256"] = decisions_hash
        atomic_json(checkpoint, value)

    started = time.monotonic()
    session_started_unix = time.time()
    write_checkpoint(len(completed), "RUNNING")
    atomic_json(
        args.pid_file,
        {
            "pid": os.getpid(),
            "started_unix": session_started_unix,
            "command": command,
            "input": str(args.input),
            "input_sha256": input_hash,
            "workers": args.workers,
            "completed_at_resume": len(completed),
            "source_sha256": source_hash,
            "reference_sha256": reference_hash,
            "checkpoint": str(checkpoint),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
    )
    print(
        f"PID {os.getpid()}; {len(graphs)} graphs; {args.workers} workers; "
        f"resume position {len(completed)}",
        file=sys.stderr,
        flush=True,
    )

    all_rows: list[Sequence[object]] = list(completed)
    remaining = graphs[len(completed):]
    with partial.open("a", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=initialize_worker,
        ) as executor:
            for group in batches(remaining, args.batch_size):
                for row in executor.map(
                    evaluate_graph, group, chunksize=args.map_chunksize
                ):
                    writer.writerow(row)
                    all_rows.append(row)
                    done = len(all_rows)
                    if done % args.progress_every == 0 or done == len(graphs):
                        stream.flush()
                        os.fsync(stream.fileno())
                        write_checkpoint(done, "RUNNING")
                        elapsed = time.monotonic() - started
                        rate = (done - len(completed)) / max(elapsed, 1e-9)
                        print(
                            f"progress {done}/{len(graphs)} "
                            f"({done/len(graphs):.2%}); {rate:.2f} graph/s; "
                            f"elapsed {elapsed:.1f}s",
                            file=sys.stderr,
                            flush=True,
                        )
    if len(all_rows) != len(graphs):
        raise SystemExit(f"wrote {len(all_rows)} decisions for {len(graphs)} graphs")
    os.replace(partial, args.decisions)
    decisions_hash = sha256(args.decisions)
    write_checkpoint(len(all_rows), "COMPLETE", decisions_hash)

    summary = summarize_rows(all_rows)
    oracle_status = "not requested"
    oracle_provenance = None
    if args.oracle_report is not None:
        oracle = json.loads(args.oracle_report.read_text(encoding="utf-8"))
        expected = set(
            oracle["individual_graph_decisions"][
                "enhanced_joint_existential_rejected"
            ]
        ) & set(indices)
        observed = set(summary["baseline_rejected_indices"])
        if observed != expected:
            raise SystemExit(
                "baseline oracle mismatch: "
                f"profiler-only={sorted(observed-expected)[:12]}, "
                f"oracle-only={sorted(expected-observed)[:12]}"
            )
        oracle_status = "PASS"
        oracle_provenance = {
            "path": str(args.oracle_report),
            "sha256": sha256(args.oracle_report),
        }

    comparisons = {}
    if args.comparison_strict_h_report is not None:
        strict_h = json.loads(
            args.comparison_strict_h_report.read_text(encoding="utf-8")
        )
        strict_h_rejected = set(
            strict_h["graph_decisions"]["marginal_h_rejected"]
        ) & set(indices)
        support_incremental = set(summary["incremental_rejected_indices"])
        comparisons["strict_H_sample"] = {
            "report": str(args.comparison_strict_h_report),
            "report_sha256": sha256(args.comparison_strict_h_report),
            "strict_H_incremental_rejected": sorted(strict_h_rejected),
            "support_propagation_incremental_rejected": sorted(
                support_incremental
            ),
            "intersection": sorted(strict_h_rejected & support_incremental),
            "support_only": sorted(support_incremental - strict_h_rejected),
            "strict_H_only": sorted(strict_h_rejected - support_incremental),
            "union": sorted(strict_h_rejected | support_incremental),
        }

    elapsed = time.monotonic() - started
    report = {
        "schema": 1,
        "method": "exact_K7_labeled_support_singleton_propagation",
        "command": command,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "input": str(args.input),
        "input_sha256": input_hash,
        "source_sha256": source_hash,
        "reference_sha256": reference_hash,
        "workers": args.workers,
        "wall_seconds": elapsed,
        "decisions": str(args.decisions),
        "decisions_sha256": decisions_hash,
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": sha256(checkpoint),
        "baseline_oracle_graph_set_agreement": oracle_status,
        "baseline_oracle": oracle_provenance,
        "comparisons": comparisons,
        **summary,
    }
    atomic_json(args.report, report)
    print(
        f"complete: baseline {summary['baseline_rejected']}, refined "
        f"{summary['refined_rejected']}, incremental "
        f"{summary['incremental_rejected']}; wall {elapsed:.1f}s; "
        f"decisions {report['decisions_sha256']}",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()
