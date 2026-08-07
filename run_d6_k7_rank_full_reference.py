#!/usr/bin/env python3
"""Fast independent full-residue decision wrapper for the frozen K7 reference.

The mathematical primitives remain in ``d6_k7_rank_reference.py``.  This
wrapper changes only quantifier scheduling: covers stop after the first exact
pass, and a graph stops after the first K7 seed for which every cap-three
cover fails.  It writes ordered restartable decisions and a compact report.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

import d6_k7_rank_reference as reference


ROOT = Path(__file__).resolve().parent
REFERENCE_PATH = ROOT / "d6_k7_rank_reference.py"
EXPECTED_REFERENCE_SHA256 = (
    "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
)
EXPECTED_FULL_INPUT_SHA256 = (
    "7097ddd5f333326bd4a4ab1c30379b054182542ca4cd3820f3cc3bc4308ec6da"
)
EXPECTED_FULL_GRAPHS = 113_136
QOS_CLASS_USER_INITIATED = 0x19
CACHE_LIMIT = 50_000
DECISION_FIELDS = (
    "index",
    "decision",
    "first_failing_seed",
    "seeds_checked",
    "covers_checked",
    "clique_failures",
    "degree_failures",
    "mask_failures",
    "basis_failures",
    "support_failures",
    "subspace_K_failures",
    "component_B_failures",
)


_support_solver: reference.SupportSolver | None = None
_zero_forcing: reference.ZeroForcingSolver | None = None
_structure_cache: dict[tuple[tuple[int, ...], int], str] | None = None


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


def initialize_worker() -> None:
    global _support_solver, _zero_forcing, _structure_cache
    set_worker_qos()
    _support_solver = reference.SupportSolver()
    _zero_forcing = reference.ZeroForcingSolver()
    _structure_cache = {}


def trim_worker_caches() -> None:
    assert _support_solver is not None
    assert _zero_forcing is not None
    assert _structure_cache is not None
    if len(_support_solver.cache) > CACHE_LIMIT:
        _support_solver.cache.clear()
    if len(_zero_forcing.cache) > CACHE_LIMIT:
        _zero_forcing.cache.clear()
    if len(_structure_cache) > CACHE_LIMIT:
        _structure_cache.clear()


def structure_failure(graph: Sequence[int], rank_upper: int) -> str | None:
    """Return the first exact equality-structure failure, with memoization."""

    assert _structure_cache is not None
    key = (tuple(graph), rank_upper)
    cached = _structure_cache.get(key)
    if cached is not None:
        return cached or None
    order = len(graph)
    answer: str | None = None
    if rank_upper < order and next(
        reference.clique_masks(graph, rank_upper + 1), None
    ) is not None:
        answer = "clique"
    else:
        complement = reference.complement_graph(graph)
        maximum_degree = max((row.bit_count() for row in complement), default=0)
        if order and maximum_degree > rank_upper - 1:
            answer = "degree"
        elif 0 < rank_upper <= order:
            for clique in reference.clique_masks(graph, rank_upper):
                compatible, _ = reference.saturating_clique_mask_compatibility(
                    graph, clique
                )
                if not compatible:
                    answer = "mask"
                    break
                compatible, _, _ = reference.basis_kernel_compatibility(
                    graph, clique
                )
                if not compatible:
                    answer = "basis"
                    break
    _structure_cache[key] = answer or ""
    return answer


def cover_failure(
    adj: Sequence[int],
    outside: Sequence[int],
    defects: Sequence[int],
    total_term_rank: int,
    zmask: int,
) -> str | None:
    """Return one exact rejection reason, or None when this cover survives."""

    assert _support_solver is not None
    assert _zero_forcing is not None
    zvertices = list(reference.bits(zmask))
    nvertices = [
        vertex for vertex in range(len(outside))
        if not (zmask & (1 << vertex))
    ]
    ndefects = [defects[vertex] for vertex in nvertices]
    term_rank = reference.matching_size(ndefects)
    k_upper, b_upper = reference.rank_upper_bounds(
        len(nvertices), term_rank, total_term_rank, len(zvertices)
    )
    global_n = [outside[vertex] for vertex in nvertices]
    graph_n = reference.induced_graph(adj, global_n)

    for first, second in combinations(range(len(nvertices)), 2):
        intersects = bool(ndefects[first] & ndefects[second])
        required = bool(graph_n[first] & (1 << second))
        if intersects != required:
            raise AssertionError("cover/intersection graph invariant failed")

    failed = structure_failure(graph_n, k_upper)
    if failed is not None:
        return failed
    if _support_solver.solve(defects[vertex] for vertex in zvertices) is None:
        return "support"
    k_zero_forcing = _zero_forcing.solve(graph_n).number
    if len(nvertices) - k_zero_forcing > k_upper:
        return "subspace_K"
    fgraph = reference.complement_graph(graph_n)
    nullity_cap = reference.component_inertia_nullity_cap(
        fgraph, _zero_forcing
    )
    if len(nvertices) - nullity_cap > b_upper:
        return "component_B"
    return None


def evaluate_graph(graph: dict) -> tuple:
    adj = tuple(int(row) for row in graph["adjacency"])
    index = int(graph["index"])
    reference.validate_graph(adj)
    seeds_checked = 0
    covers_checked = 0
    failures: Counter[str] = Counter()
    any_seed = False
    for seed_mask in reference.clique_masks(adj, 7):
        any_seed = True
        seeds_checked += 1
        _, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        seed_passes = False
        for cover in reference.eligible_covers(ladj, eligible, cap=3):
            covers_checked += 1
            failure = cover_failure(
                adj, outside, defects, total_term_rank, cover
            )
            if failure is None:
                seed_passes = True
                break
            failures[failure] += 1
        if not seed_passes:
            trim_worker_caches()
            return (
                index,
                "REJECTED",
                seed_mask,
                seeds_checked,
                covers_checked,
                *(failures[name] for name in (
                    "clique",
                    "degree",
                    "mask",
                    "basis",
                    "support",
                    "subspace_K",
                    "component_B",
                )),
            )
    if not any_seed:
        raise AssertionError(f"selected graph {index} has no K7")
    trim_worker_caches()
    return (
        index,
        "SURVIVOR",
        0,
        seeds_checked,
        covers_checked,
        *(failures[name] for name in (
            "clique",
            "degree",
            "mask",
            "basis",
            "support",
            "subspace_K",
            "component_B",
        )),
    )


def load_completed(path: Path, expected: Sequence[int]) -> list[tuple[str, ...]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="ascii") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, None)
        if header != list(DECISION_FIELDS):
            raise ValueError(f"bad partial decision header in {path}")
        rows = [tuple(row) for row in reader]
    if len(rows) > len(expected):
        raise ValueError("partial decision file is longer than the input")
    for position, row in enumerate(rows):
        if int(row[0]) != expected[position]:
            raise ValueError(
                f"partial decision index {row[0]} at position {position}"
            )
    return rows


def batches(values: Sequence[dict], size: int) -> Iterator[Sequence[dict]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def summarize_rows(rows: Sequence[Sequence[object]]) -> dict:
    rejected = [int(row[0]) for row in rows if row[1] == "REJECTED"]
    survivors = [int(row[0]) for row in rows if row[1] == "SURVIVOR"]
    return {
        "graphs": len(rows),
        "rejected": len(rejected),
        "survivors": len(survivors),
        "rejected_indices": rejected,
        "survivor_indices": survivors,
        "seeds_checked": sum(int(row[3]) for row in rows),
        "covers_checked": sum(int(row[4]) for row in rows),
        "checked_cover_failure_reasons": {
            name: sum(int(row[position]) for row in rows)
            for position, name in enumerate(
                (
                    "clique",
                    "degree",
                    "mask",
                    "basis",
                    "support",
                    "subspace_K",
                    "component_B",
                ),
                start=5,
            )
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path("/tmp/d6_k7_rank_full_residue.json")
    )
    parser.add_argument("--workers", type=int, default=min(24, os.cpu_count() or 1))
    parser.add_argument(
        "--decisions", type=Path,
        default=Path("/tmp/d6_k7_rank_python_full_decisions.tsv"),
    )
    parser.add_argument(
        "--report", type=Path,
        default=Path("/tmp/d6_k7_rank_python_full_report.json"),
    )
    parser.add_argument(
        "--pid-file", type=Path,
        default=Path("/tmp/d6_k7_rank_python_full.pid.json"),
    )
    parser.add_argument("--oracle-report", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--batch-size", type=int, default=384)
    parser.add_argument("--map-chunksize", type=int, default=4)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()

    actual_reference = sha256(REFERENCE_PATH)
    if actual_reference != EXPECTED_REFERENCE_SHA256:
        raise SystemExit(
            f"frozen reference hash is {actual_reference}, expected "
            f"{EXPECTED_REFERENCE_SHA256}"
        )
    input_hash = sha256(args.input)
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    graphs = payload["graphs"]
    if args.input.name == "d6_k7_rank_full_residue.json":
        if input_hash != EXPECTED_FULL_INPUT_SHA256:
            raise SystemExit(f"full input hash is {input_hash}")
        if len(graphs) != EXPECTED_FULL_GRAPHS:
            raise SystemExit(f"full input has {len(graphs)} graphs")
    indices = [int(graph["index"]) for graph in graphs]
    if len(indices) != len(set(indices)):
        raise SystemExit("input contains duplicate graph indices")

    partial = args.decisions.with_name(args.decisions.name + ".partial")
    completed = load_completed(partial, indices) if args.resume else []
    if partial.exists() and not args.resume:
        partial.unlink()
    mode = "a" if completed else "w"
    started = time.monotonic()
    atomic_json(
        args.pid_file,
        {
            "pid": os.getpid(),
            "started_unix": time.time(),
            "input": str(args.input),
            "input_sha256": input_hash,
            "workers": args.workers,
            "completed_at_resume": len(completed),
        },
    )
    print(
        f"PID {os.getpid()}; {len(graphs)} graphs; {args.workers} workers; "
        f"resume position {len(completed)}",
        file=sys.stderr,
        flush=True,
    )

    all_rows: list[Sequence[object]] = list(completed)
    remaining = graphs[len(completed) :]
    with partial.open(mode, newline="", encoding="ascii") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        if not completed:
            writer.writerow(DECISION_FIELDS)
            stream.flush()
            os.fsync(stream.fileno())
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
                        elapsed = time.monotonic() - started
                        rate = (done - len(completed)) / max(elapsed, 1e-9)
                        remaining_seconds = (len(graphs) - done) / max(rate, 1e-9)
                        print(
                            f"progress {done}/{len(graphs)} "
                            f"({done/len(graphs):.2%}); {rate:.2f} graph/s; "
                            f"elapsed {elapsed:.1f}s; projected remaining "
                            f"{remaining_seconds:.1f}s",
                            file=sys.stderr,
                            flush=True,
                        )
    if len(all_rows) != len(graphs):
        raise SystemExit(f"wrote {len(all_rows)} decisions for {len(graphs)} graphs")
    os.replace(partial, args.decisions)

    summary = summarize_rows(all_rows)
    if args.oracle_report is not None:
        oracle = json.loads(args.oracle_report.read_text(encoding="utf-8"))
        expected = set(
            oracle["individual_graph_decisions"][
                "enhanced_joint_existential_rejected"
            ]
        )
        observed = set(summary["rejected_indices"])
        if observed != expected:
            raise SystemExit(
                "oracle mismatch: "
                f"wrapper-only={sorted(observed-expected)[:12]}, "
                f"oracle-only={sorted(expected-observed)[:12]}"
            )
        oracle_status = "PASS"
    else:
        oracle_status = "not requested"
    elapsed = time.monotonic() - started
    report = {
        "schema": 1,
        "method": "frozen_Python_K7_reference_cap3_short_circuit",
        "input": str(args.input),
        "input_sha256": input_hash,
        "reference_sha256": actual_reference,
        "workers": args.workers,
        "wall_seconds": elapsed,
        "decisions": str(args.decisions),
        "decisions_sha256": sha256(args.decisions),
        "oracle_graph_set_agreement": oracle_status,
        **summary,
    }
    atomic_json(args.report, report)
    print(
        f"complete: rejected {summary['rejected']}, survivors "
        f"{summary['survivors']}; wall {elapsed:.1f}s; "
        f"decisions {report['decisions_sha256']}",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()
