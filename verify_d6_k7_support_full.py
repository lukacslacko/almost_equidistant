#!/usr/bin/env python3
"""Independently verify the archived full K7 support-propagation rejections.

The production support-propagation module is deliberately neither imported
nor called.  This checker uses the frozen rank reference only to identify the
baseline-passing covers.  It independently enumerates every *labeled* actual
support family for those covers and independently computes singleton-overlap
closure.  A recorded rejection is accepted only when every support family of
every baseline-passing cover at its recorded failing seed ends in either an
empty propagated mask or a disjoint required edge.

The checker also binds the compressed archive to its uncompressed decision
hash, the complete 17,764-graph input population, the production report, the
frozen reference, the production source, and (when present) its checkpoint.
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
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

import d6_k7_rank_reference as rank_reference


ROOT = Path(__file__).resolve().parent
EXPECTED_REFERENCE_SHA256 = (
    "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
)
EXPECTED_UNCOMPRESSED_DECISIONS_SHA256 = (
    "117e057a79135fd48cac316ee4f1e08e11a2be84e00d9d10393bc17d4a518a31"
)
EXPECTED_GRAPH_COUNT = 17_764
EXPECTED_REJECTED_COUNT = 1_536
QOS_CLASS_USER_INITIATED = 0x19

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
    "assignment_fail_empty_propagated_mask",
    "assignment_fail_disjoint_required_edge",
    "assignment_fail_clique",
    "assignment_fail_degree",
    "assignment_fail_basis",
    "assignment_fail_mask",
    "assignment_fail_subspace_K",
    "assignment_fail_component_B",
)
TOTAL_FIELDS = DECISION_FIELDS[6:]
SIMPLE_FAILURES = ("empty_propagated_mask", "disjoint_required_edge")


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
    """Keep campaign workers out of macOS background QoS."""

    if sys.platform != "darwin":
        return
    libc = ctypes.CDLL(None)
    setter = libc.pthread_set_qos_class_self_np
    setter.argtypes = (ctypes.c_uint, ctypes.c_int)
    setter.restype = ctypes.c_int
    error = int(setter(QOS_CLASS_USER_INITIATED, 0))
    if error:
        if os.environ.get("D6_ALLOW_BACKGROUND_TEST_ONLY") == "1":
            return
        raise OSError(error, "pthread_set_qos_class_self_np failed")


def bit_positions(mask: int) -> Iterator[int]:
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def has_injective_representatives(supports: Sequence[int]) -> bool:
    """Check Hall feasibility via augmenting paths, independently of term rank."""

    coordinate_owner = [-1] * 7

    def augment(column: int, seen: list[bool]) -> bool:
        for coordinate in bit_positions(supports[column] & 0x7F):
            if seen[coordinate]:
                continue
            seen[coordinate] = True
            owner = coordinate_owner[coordinate]
            if owner < 0 or augment(owner, seen):
                coordinate_owner[coordinate] = column
                return True
        return False

    return all(augment(column, [False] * 7) for column in range(len(supports)))


def independently_valid_support_family(supports: Sequence[int]) -> bool:
    """Reimplement the exact support rules needed for cap-three covers."""

    if len(supports) > 3:
        raise ValueError("this checker is specialized to cap-three covers")
    if any(support <= 0 or support & ~0x7F or support.bit_count() < 3
           for support in supports):
        return False
    if any((first & second).bit_count() == 1
           for first, second in combinations(supports, 2)):
        return False
    return has_injective_representatives(supports)


def independent_support_domains(allowed: int) -> tuple[int, ...]:
    if allowed & ~0x7F:
        raise ValueError("allowed support mask exceeds seven coordinates")
    return tuple(
        support
        for support in range(1, 128)
        if not (support & ~allowed) and support.bit_count() >= 3
    )


def independent_labeled_support_families(
    allowed_masks: Sequence[int],
) -> Iterator[tuple[int, ...]]:
    """Enumerate all labeled families; equal domains are not quotiented."""

    if len(allowed_masks) > 3:
        raise ValueError("only cap-three covers are in scope")
    domains = tuple(independent_support_domains(mask) for mask in allowed_masks)
    if any(not domain for domain in domains):
        return
    assignment = [0] * len(domains)
    order = tuple(sorted(range(len(domains)), key=lambda i: (len(domains[i]), i)))

    def visit(depth: int) -> Iterator[tuple[int, ...]]:
        if depth == len(order):
            family = tuple(assignment)
            if independently_valid_support_family(family):
                yield family
            return
        variable = order[depth]
        for support in domains[variable]:
            if any(
                assignment[prior]
                and (support & assignment[prior]).bit_count() == 1
                for prior in order[:depth]
            ):
                continue
            assignment[variable] = support
            chosen = tuple(assignment[i] for i in order[: depth + 1])
            if has_injective_representatives(chosen):
                yield from visit(depth + 1)
        assignment[variable] = 0

    yield from visit(0)


def independent_singleton_closure(
    allowed: int, fixed_supports: Sequence[int]
) -> tuple[int, int]:
    """Delete singleton overlaps to a fixed point using a fresh work-list."""

    current = int(allowed)
    deletions = 0
    while True:
        deletion = 0
        for support in fixed_supports:
            overlap = current & support
            if overlap.bit_count() == 1:
                deletion |= overlap
        if not deletion:
            return current, deletions
        current &= ~deletion
        deletions += deletion.bit_count()


def independently_propagate_family(
    graph_n: Sequence[int],
    n_allowed_masks: Sequence[int],
    fixed_supports: Sequence[int],
) -> tuple[str | None, int]:
    """Return the simple exact failure, or None if this family survives it."""

    propagated: list[int] = []
    deletions = 0
    for allowed in n_allowed_masks:
        closed, count = independent_singleton_closure(allowed, fixed_supports)
        propagated.append(closed)
        deletions += count
    if any(mask == 0 for mask in propagated):
        return "empty_propagated_mask", deletions
    for first, second in combinations(range(len(propagated)), 2):
        if (
            graph_n[first] & (1 << second)
            and not (propagated[first] & propagated[second])
        ):
            return "disjoint_required_edge", deletions
    return None, deletions


def independent_seed_instance(
    adj: Sequence[int], seed_mask: int
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...], int]:
    """Construct outside labels, defect masks, L, and eligible Z vertices."""

    seed = tuple(bit_positions(seed_mask))
    if len(seed) != 7:
        raise ValueError("recorded seed does not have seven vertices")
    if any(not (adj[first] & (1 << second))
           for first, second in combinations(seed, 2)):
        raise ValueError("recorded seed is not a required K7")
    outside = tuple(v for v in range(len(adj)) if not (seed_mask & (1 << v)))
    defects = tuple(
        sum(1 << coordinate for coordinate, q in enumerate(seed)
            if not (adj[x] & (1 << q)))
        for x in outside
    )
    ladj = [0] * len(outside)
    for first, second in combinations(range(len(outside)), 2):
        if (
            adj[outside[first]] & (1 << outside[second])
            and not (defects[first] & defects[second])
        ):
            ladj[first] |= 1 << second
            ladj[second] |= 1 << first
    eligible = sum(
        1 << vertex
        for vertex, allowed in enumerate(defects)
        if allowed.bit_count() >= 3
    )
    return outside, defects, tuple(ladj), eligible


def independent_eligible_covers(
    ladj: Sequence[int], eligible: int
) -> Iterator[int]:
    """Enumerate all vertex covers of L contained in eligible, of size <=3."""

    full = (1 << len(ladj)) - 1
    for zmask in range(1 << len(ladj)):
        if zmask & ~eligible or zmask.bit_count() > 3:
            continue
        remaining = full & ~zmask
        if all(not (ladj[vertex] & remaining)
               for vertex in bit_positions(remaining)):
            yield zmask


def independent_induced_graph(
    adj: Sequence[int], selected: Sequence[int]
) -> tuple[int, ...]:
    position = {vertex: i for i, vertex in enumerate(selected)}
    return tuple(
        sum(1 << position[other] for other in selected
            if adj[vertex] & (1 << other))
        for vertex in selected
    )


_support_solver: rank_reference.SupportSolver | None = None
_zero_forcing: rank_reference.ZeroForcingSolver | None = None
_clique_solver: rank_reference.CliqueStructureSolver | None = None


def initialize_worker() -> None:
    global _support_solver, _zero_forcing, _clique_solver
    set_worker_qos()
    _support_solver = rank_reference.SupportSolver()
    _zero_forcing = rank_reference.ZeroForcingSolver()
    _clique_solver = rank_reference.CliqueStructureSolver()


def verify_rejected_graph(task: tuple[dict, int]) -> dict:
    """Exhaust the independently reconstructed failing seed of one graph."""

    assert _support_solver is not None
    assert _zero_forcing is not None
    assert _clique_solver is not None
    graph, seed_mask = task
    index = int(graph["index"])
    adj = tuple(int(row) for row in graph["adjacency"])
    rank_reference.validate_graph(adj)
    outside, defects, ladj, eligible = independent_seed_instance(adj, seed_mask)

    # Cross-check only the seed construction at the frozen baseline boundary.
    _, reference_outside, reference_defects, reference_ladj, reference_eligible = (
        rank_reference.seed_instance(adj, seed_mask)
    )
    if (
        tuple(reference_outside) != outside
        or tuple(reference_defects) != defects
        or tuple(reference_ladj) != ladj
        or reference_eligible != eligible
    ):
        raise AssertionError(f"independent seed construction mismatch at {index}")
    covers = tuple(independent_eligible_covers(ladj, eligible))
    if covers != tuple(rank_reference.eligible_covers(ladj, eligible, cap=3)):
        raise AssertionError(f"independent cover enumeration mismatch at {index}")

    total_term_rank = rank_reference.matching_size(defects)
    baseline_passing = 0
    families = 0
    deletions = 0
    failures: Counter[str] = Counter()
    cover_summaries: list[dict] = []
    for zmask in covers:
        baseline = rank_reference.analyze_cover(
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
        baseline_passing += 1
        zvertices = tuple(bit_positions(zmask))
        nvertices = tuple(
            vertex for vertex in range(len(outside))
            if not (zmask & (1 << vertex))
        )
        z_allowed = tuple(defects[vertex] for vertex in zvertices)
        n_allowed = tuple(defects[vertex] for vertex in nvertices)
        graph_n = independent_induced_graph(
            adj, tuple(outside[vertex] for vertex in nvertices)
        )
        cover_families = 0
        cover_failures: Counter[str] = Counter()
        cover_deletions = 0
        for fixed_supports in independent_labeled_support_families(z_allowed):
            cover_families += 1
            failure, count = independently_propagate_family(
                graph_n, n_allowed, fixed_supports
            )
            cover_deletions += count
            if failure is None:
                raise AssertionError(
                    f"graph {index}, seed {seed_mask}, cover {zmask} has "
                    f"an independently surviving support family {fixed_supports}"
                )
            cover_failures[failure] += 1
        if cover_families == 0:
            raise AssertionError(
                f"baseline-passing cover {zmask} at graph {index} has no family"
            )
        families += cover_families
        deletions += cover_deletions
        failures.update(cover_failures)
        cover_summaries.append(
            {
                "zmask": zmask,
                "z_size": zmask.bit_count(),
                "families": cover_families,
                "failures": dict(sorted(cover_failures.items())),
                "coordinate_deletions": cover_deletions,
            }
        )
    if baseline_passing == 0:
        raise AssertionError(
            f"recorded refined-only failure at graph {index} has no baseline pass"
        )
    if families != sum(failures.values()):
        raise AssertionError(f"family population mismatch at graph {index}")
    return {
        "index": index,
        "recorded_failing_seed": seed_mask,
        "eligible_covers": len(covers),
        "baseline_passing_covers": baseline_passing,
        "support_families_exhausted": families,
        "failure_counts": dict(sorted(failures.items())),
        "coordinate_deletions": deletions,
        "covers": cover_summaries,
    }


def load_and_validate_archive(
    archive: Path,
    production_report_path: Path,
    input_path: Path,
    reference_path: Path,
    production_source: Path,
) -> tuple[list[dict], list[dict], dict]:
    """Validate all immutable bindings before running the decision replay."""

    report = json.loads(production_report_path.read_text(encoding="utf-8"))
    if sha256(reference_path) != EXPECTED_REFERENCE_SHA256:
        raise ValueError("frozen reference hash mismatch")
    if report.get("reference_sha256") != EXPECTED_REFERENCE_SHA256:
        raise ValueError("production report reference binding mismatch")
    input_hash = sha256(input_path)
    if report.get("input_sha256") != input_hash:
        raise ValueError("production report input binding mismatch")
    source_hash = sha256(production_source)
    if report.get("source_sha256") != source_hash:
        raise ValueError("production report source binding mismatch")

    compressed_hash = sha256(archive)
    with gzip.open(archive, "rb") as stream:
        uncompressed = stream.read()
    uncompressed_hash = hashlib.sha256(uncompressed).hexdigest()
    if uncompressed_hash != EXPECTED_UNCOMPRESSED_DECISIONS_SHA256:
        raise ValueError("archived uncompressed decision hash mismatch")
    if report.get("decisions_sha256") != uncompressed_hash:
        raise ValueError("production report decision binding mismatch")
    reader = csv.DictReader(
        io.StringIO(uncompressed.decode("ascii")), delimiter="\t"
    )
    if tuple(reader.fieldnames or ()) != DECISION_FIELDS:
        raise ValueError("decision archive header mismatch")
    rows: list[dict] = []
    for position, raw in enumerate(reader, start=1):
        row = {
            field: (
                raw[field]
                if field in ("baseline_decision", "refined_decision")
                else int(raw[field])
            )
            for field in DECISION_FIELDS
        }
        if row["baseline_decision"] not in {"SURVIVOR", "REJECTED"}:
            raise ValueError(f"bad baseline decision at row {position}")
        if row["refined_decision"] not in {"SURVIVOR", "REJECTED"}:
            raise ValueError(f"bad refined decision at row {position}")
        rows.append(row)

    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list):
        raise ValueError("input has no graph list")
    input_indices = [int(graph["index"]) for graph in graphs]
    archive_indices = [row["index"] for row in rows]
    if archive_indices != input_indices or len(set(input_indices)) != len(input_indices):
        raise ValueError("archive/input graph population or order mismatch")
    if len(graphs) != EXPECTED_GRAPH_COUNT:
        raise ValueError("unexpected input population")
    if any(row["applicable"] != 1 for row in rows):
        raise ValueError("rank-survivor population contains a non-K7 row")

    baseline_rejected = [
        row["index"] for row in rows if row["baseline_decision"] == "REJECTED"
    ]
    refined_rejected = [
        row["index"] for row in rows if row["refined_decision"] == "REJECTED"
    ]
    if baseline_rejected:
        raise ValueError("rank-survivor input has a baseline rejection")
    if len(refined_rejected) != EXPECTED_REJECTED_COUNT:
        raise ValueError("unexpected refined rejection count")
    if any(
        bool(row["first_refined_failing_seed"])
        != (row["refined_decision"] == "REJECTED")
        for row in rows
    ):
        raise ValueError("recorded failing seed/decision mismatch")

    totals = {
        field: sum(int(row[field]) for row in rows)
        for field in TOTAL_FIELDS
    }
    expected_scalars = {
        "graphs": len(rows),
        "applicable_K7_graphs": len(rows),
        "baseline_rejected": len(baseline_rejected),
        "baseline_survivors": len(rows) - len(baseline_rejected),
        "refined_rejected": len(refined_rejected),
        "refined_survivors": len(rows) - len(refined_rejected),
        "incremental_rejected": len(set(refined_rejected) - set(baseline_rejected)),
    }
    for field, value in expected_scalars.items():
        if report.get(field) != value:
            raise ValueError(f"report scalar mismatch for {field}")
    for field, observed in totals.items():
        if report.get("totals", {}).get(field) != observed:
            raise ValueError(f"report total mismatch for {field}")
    if report.get("baseline_rejected_indices") != baseline_rejected:
        raise ValueError("report baseline rejection list mismatch")
    if report.get("refined_rejected_indices") != refined_rejected:
        raise ValueError("report refined rejection list mismatch")
    if report.get("incremental_rejected_indices") != refined_rejected:
        raise ValueError("report incremental rejection list mismatch")
    if any(totals[f"assignment_fail_{reason}"] <= 0 for reason in SIMPLE_FAILURES):
        raise ValueError("production archive lacks an expected simple failure class")
    for field in TOTAL_FIELDS:
        if field.startswith("assignment_fail_") and field not in {
            "assignment_fail_empty_propagated_mask",
            "assignment_fail_disjoint_required_edge",
        } and totals[field] != 0:
            raise ValueError("archive reports a non-simple propagation failure")

    checkpoint_status: dict | None = None
    checkpoint = ROOT / str(report.get("checkpoint", ""))
    if checkpoint.is_file():
        checkpoint_hash = sha256(checkpoint)
        if checkpoint_hash != report.get("checkpoint_sha256"):
            raise ValueError("checkpoint hash mismatch")
        checkpoint_status = {
            "path": str(checkpoint),
            "sha256": checkpoint_hash,
            "validated": True,
        }

    sources_status: dict | None = None
    sources = payload.get("sources", {})
    source_archive_name = sources.get("decision_archive")
    if source_archive_name:
        source_archive = ROOT / str(source_archive_name)
        if not source_archive.is_file():
            raise ValueError("input population source archive is absent")
        observed = sha256(source_archive)
        if observed != sources.get("decision_archive_sha256"):
            raise ValueError("input population source archive hash mismatch")
        sources_status = {
            "decision_archive": str(source_archive),
            "decision_archive_sha256": observed,
            "declared_decision_rows": sources.get("decision_rows"),
        }

    bindings = {
        "production_report": str(production_report_path),
        "production_report_sha256": sha256(production_report_path),
        "archive": str(archive),
        "archive_compressed_sha256": compressed_hash,
        "archive_uncompressed_sha256": uncompressed_hash,
        "archive_uncompressed_bytes": len(uncompressed),
        "archive_rows": len(rows),
        "input": str(input_path),
        "input_sha256": input_hash,
        "input_graphs": len(graphs),
        "reference": str(reference_path),
        "reference_sha256": EXPECTED_REFERENCE_SHA256,
        "production_source": str(production_source),
        "production_source_sha256": source_hash,
        "checkpoint": checkpoint_status,
        "population_source": sources_status,
        "reported_totals_recomputed": totals,
        "rejected_count": len(refined_rejected),
        "rejected_indices_sha256": hashlib.sha256(
            ("\n".join(map(str, refined_rejected)) + "\n").encode("ascii")
        ).hexdigest(),
    }
    return graphs, rows, bindings


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / ".runs/d6_k7_rank_survivors.json"
    )
    parser.add_argument(
        "--archive", type=Path,
        default=ROOT / "d6_k7_support_rank_survivors_decisions.tsv.gz",
    )
    parser.add_argument(
        "--production-report", type=Path,
        default=ROOT / "d6_k7_support_rank_survivors_report.json",
    )
    parser.add_argument(
        "--reference", type=Path, default=ROOT / "d6_k7_rank_reference.py"
    )
    parser.add_argument(
        "--production-source", type=Path,
        default=ROOT / "d6_k7_support_propagation.py",
    )
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_support_full_verification_report.json",
    )
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument(
        "--limit", type=int,
        help="verify only the first N recorded rejections (bounded control)",
    )
    parser.add_argument("--progress-every", type=int, default=64)
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("workers must be positive")
    if args.limit is not None and args.limit < 1:
        raise SystemExit("limit must be positive")

    started_wall = time.time()
    started = time.monotonic()
    graphs, rows, bindings = load_and_validate_archive(
        args.archive,
        args.production_report,
        args.input,
        args.reference,
        args.production_source,
    )
    graph_by_index = {int(graph["index"]): graph for graph in graphs}
    rejected_rows = [
        row for row in rows if row["refined_decision"] == "REJECTED"
    ]
    selected_rows = (
        rejected_rows if args.limit is None else rejected_rows[: args.limit]
    )
    tasks = [
        (graph_by_index[row["index"]], row["first_refined_failing_seed"])
        for row in selected_rows
    ]
    command = [sys.executable, *sys.argv]
    results: list[dict] = []
    if args.workers == 1:
        initialize_worker()
        iterator = map(verify_rejected_graph, tasks)
        for position, result in enumerate(iterator, start=1):
            results.append(result)
            if position % args.progress_every == 0 or position == len(tasks):
                elapsed = time.monotonic() - started
                print(
                    f"verified {position}/{len(tasks)}; "
                    f"{position/max(elapsed, 1e-9):.2f} graph/s; "
                    f"elapsed {elapsed:.1f}s",
                    file=sys.stderr,
                    flush=True,
                )
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=initialize_worker
        ) as executor:
            for position, result in enumerate(
                executor.map(verify_rejected_graph, tasks, chunksize=1), start=1
            ):
                results.append(result)
                if position % args.progress_every == 0 or position == len(tasks):
                    elapsed = time.monotonic() - started
                    print(
                        f"verified {position}/{len(tasks)}; "
                        f"{position/max(elapsed, 1e-9):.2f} graph/s; "
                        f"elapsed {elapsed:.1f}s",
                        file=sys.stderr,
                        flush=True,
                    )

    if [item["index"] for item in results] != [row["index"] for row in selected_rows]:
        raise SystemExit("parallel replay changed result order or population")
    total_families = sum(item["support_families_exhausted"] for item in results)
    total_failures: Counter[str] = Counter()
    for item in results:
        total_failures.update(item["failure_counts"])
    if total_families != sum(total_failures.values()):
        raise SystemExit("global family/failure count mismatch")
    elapsed = time.monotonic() - started
    complete = len(results) == EXPECTED_REJECTED_COUNT
    report = {
        "schema": 1,
        "method": "independent_K7_labeled_support_singleton_rejection_replay",
        "status": "PASS" if complete else "BOUNDED_CONTROL_PASS",
        "complete_rejection_population_verified": complete,
        "command": command,
        "started_unix": started_wall,
        "wall_seconds": elapsed,
        "workers": args.workers,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "verifier_source_sha256": sha256(Path(__file__)),
        "production_module_imported": False,
        "bindings": bindings,
        "recorded_rejections": len(rejected_rows),
        "selected_rejections": len(selected_rows),
        "verified_rejections": len(results),
        "baseline_passing_covers_exhausted": sum(
            item["baseline_passing_covers"] for item in results
        ),
        "labeled_support_families_exhausted": total_families,
        "independent_failure_counts": dict(sorted(total_failures.items())),
        "coordinate_deletions": sum(
            item["coordinate_deletions"] for item in results
        ),
        "results": results,
        "trust_scope": (
            "Exact integer/bit-mask replay. The frozen rank reference is "
            "trusted only for baseline-cover classification; labeled support "
            "enumeration and singleton propagation are independently coded."
        ),
    }
    atomic_json(args.report, report)
    print(
        f"{report['status']}: {len(results)} rejections, "
        f"{report['baseline_passing_covers_exhausted']} covers, "
        f"{total_families} labeled families; wall {elapsed:.1f}s; "
        f"report {sha256(args.report)}",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()
