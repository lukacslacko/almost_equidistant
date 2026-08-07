#!/usr/bin/env python3
"""Restartable certified interval search on the exact dimension-six residue.

For one placement order, ``KILLED`` is returned only when every requested
theta slice is certified empty by the outward-rounded C kernel.  One such
order certifies the graph.  ``ABORT`` (node cap), ``UNRESOLVED`` (survivor
boxes or no eligible order), and ``INFRA_ERROR`` are recorded separately and
make no mathematical claim.

Every graph result is atomically checkpointed as soon as it finishes.  A
matching checkpoint directory resumes automatically; its campaign hash pins
all mathematical parameters, selected indices, sources, and input artifacts.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import fcntl
import gzip
import hashlib
import json
import math
import os
import platform
import random
import shlex
import subprocess
import sys
import time
import traceback
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from threading import Lock
from typing import Iterable, Sequence, TextIO

import cdriver6


ROOT = Path(__file__).resolve().parent
N = 19
EXPECTED_CORPUS_GRAPHS = 3_971_787
EXPECTED_CORPUS_SHA256 = (
    "12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5"
)
EXPECTED_K7_INPUT_SHA256 = {
    # Exhaustive C TSV and its tracked deterministic Python decision archive.
    "eaa4d5061b44cc86e6ab98ac6315b29a509539e5ac22ab33d1b7a03a34358391",
    "d331016042c14ba412a29e42b1f2ee06ee101a7d4376067b64e99e249677f10f",
}
EXPECTED_K6_SUPPORT_REPORT_SHA256 = (
    "8f344566b9fa84d042a05712bf48892f013510740141471ac16e2f186fb95c14"
)
EXPECTED_K7_STRICT_H_REPORT_SHA256 = (
    "59be9c6a2d4cd2e8e8214d28d3198487846e58697046415ae1e15f4f69dc3add"
)
EXPECTED_K7_SUPPORT_REPORT_SHA256 = {
    # Audited 512-graph sample.  Add the frozen full report in the commit that
    # authorizes it for production selection.
    "d0743b4e0977fa9333e152d6c290b1ca2e5d121ac0c9d829e092093e547c58f0",
    # Exact 17,764-graph rank-survivor pass; independent rejection replay is
    # archived alongside this runner before production use.
    "5138d207543967b393baae9b14b47b49a276273cbc4468ee843856891917691e",
}
EXPECTED_K6_BIPARTITE_REPORT_SHA256 = (
    "ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc"
)
EXPECTED_BASE_COUNTS = {"K7": 17_764, "K6": 1_098}
RESULT_STATUSES = ("KILLED", "ABORT", "UNRESOLVED", "INFRA_ERROR")
KERNEL_STATUSES = ("KILLED", "SURVIVORS", "ABORT")
QOS_CLASS_USER_INITIATED = 0x19
_KERNEL_LOCK = Lock()


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


def git_provenance() -> dict:
    """Return the exact commit/branch and a hashed porcelain worktree state."""

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


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    try:
        directory = os.open(path.parent, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def open_ascii(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open(encoding="ascii", newline="")


def load_k7_survivors(path: Path) -> tuple[set[int], dict]:
    actual = sha256(path)
    if actual not in EXPECTED_K7_INPUT_SHA256:
        raise ValueError(f"unexpected K7 decision SHA-256 {actual}")
    survivors: set[int] = set()
    rows = 0
    with open_ascii(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if "index" not in fields:
            raise ValueError("K7 decisions have no index column")
        if "joint_rejected" in fields:
            format_name = "exhaustive_C_TSV"

            def survives(row: dict[str, str]) -> bool:
                return int(row["joint_rejected"]) == 0

        elif "decision" in fields:
            format_name = "short_circuit_Python_TSV"

            def survives(row: dict[str, str]) -> bool:
                if row["decision"] not in ("REJECTED", "SURVIVOR"):
                    raise ValueError(f"bad K7 decision {row['decision']!r}")
                return row["decision"] == "SURVIVOR"

        else:
            raise ValueError("unrecognized K7 decision columns")
        seen: set[int] = set()
        for row in reader:
            index = int(row["index"])
            if index in seen:
                raise ValueError(f"duplicate K7 decision index {index}")
            seen.add(index)
            rows += 1
            if survives(row):
                survivors.add(index)
    if rows != 113_136 or len(survivors) != EXPECTED_BASE_COUNTS["K7"]:
        raise ValueError(
            f"unexpected K7 decisions rows={rows}, survivors={len(survivors)}"
        )
    return survivors, {
        "kind": "K7_enhanced_rank_survivors",
        "path": str(path),
        "sha256": actual,
        "format": format_name,
        "input_rows": rows,
        "survivors": len(survivors),
    }


def load_k6_survivors(path: Path) -> tuple[set[int], dict]:
    actual = sha256(path)
    if actual != EXPECTED_K6_SUPPORT_REPORT_SHA256:
        raise ValueError(f"unexpected K6 support report SHA-256 {actual}")
    report = json.loads(path.read_text(encoding="utf-8"))
    results = report.get("graph_results")
    if not isinstance(results, list):
        raise ValueError("K6 support report has no graph_results")
    survivors: set[int] = set()
    seen: set[int] = set()
    for result in results:
        index = int(result["index"])
        if index in seen:
            raise ValueError(f"duplicate K6 report index {index}")
        seen.add(index)
        decision = result.get("decision", {})
        if not bool(decision.get("rejected")):
            survivors.add(index)
    if len(survivors) != EXPECTED_BASE_COUNTS["K6"]:
        raise ValueError(f"unexpected K6 support survivors {len(survivors)}")
    return survivors, {
        "kind": "K6_support_survivors",
        "path": str(path),
        "sha256": actual,
        "input_rows": len(results),
        "survivors": len(survivors),
    }


def apply_k7_strict_h(
    selected: dict[int, str], k7_universe: set[int], path: Path
) -> dict:
    actual = sha256(path)
    if actual != EXPECTED_K7_STRICT_H_REPORT_SHA256:
        raise ValueError(f"unapproved strict-H report SHA-256 {actual}")
    report = json.loads(path.read_text(encoding="utf-8"))
    decisions = report.get("decisions", {})
    rejected = set(map(int, decisions.get("strict_H_rejected", ())))
    survivors = set(map(int, decisions.get("strict_H_survivors", ())))
    if rejected & survivors or rejected | survivors != k7_universe:
        raise ValueError("strict-H report is not a partition of the K7 base")
    before = {index for index, population in selected.items() if population == "K7"}
    if before != k7_universe:
        raise ValueError("strict-H must be the first optional K7 selection layer")
    for index in rejected:
        selected.pop(index)
    return {
        "kind": "K7_strict_H_survivors",
        "path": str(path),
        "sha256": actual,
        "before": len(before),
        "removed": len(rejected),
        "after": len(survivors),
        "removed_indices_sha256": stable_hash(sorted(rejected)),
    }


def apply_k7_support_report(
    selected: dict[int, str], k7_universe: set[int], path: Path
) -> dict:
    """Remove exact K7 support-propagation rejections from any prior layer.

    A pilot report may cover only a sample, so this hook consumes its explicit
    incremental rejection list rather than assuming unlisted graphs passed.
    """

    actual = sha256(path)
    if actual not in EXPECTED_K7_SUPPORT_REPORT_SHA256:
        raise ValueError(f"unapproved K7 support report SHA-256 {actual}")
    report = json.loads(path.read_text(encoding="utf-8"))
    for key in (
        "incremental_rejected_indices",
        "rejected_indices",
        "refined_rejected_indices",
    ):
        if key in report:
            rejected = set(map(int, report[key]))
            break
    else:
        raise ValueError("K7 support report has no explicit rejection indices")
    if not rejected <= k7_universe:
        raise ValueError("K7 support report contains a non-K7-base index")
    if len(rejected) != len(report[key]):
        raise ValueError("K7 support report repeats a rejection index")
    before = {index for index, population in selected.items() if population == "K7"}
    removed = before & rejected
    for index in removed:
        selected.pop(index)
    return {
        "kind": "K7_support_propagation_rejections",
        "path": str(path),
        "sha256": actual,
        "decision_field": key,
        "report_graphs": report.get("graphs"),
        "before": len(before),
        "report_rejections": len(rejected),
        "removed": len(removed),
        "after": len(before - removed),
        "removed_indices_sha256": stable_hash(sorted(removed)),
    }


def apply_k6_bipartite_report(
    selected: dict[int, str], k6_universe: set[int], path: Path
) -> dict:
    actual = sha256(path)
    if actual != EXPECTED_K6_BIPARTITE_REPORT_SHA256:
        raise ValueError(f"unapproved K6 bipartite report SHA-256 {actual}")
    report = json.loads(path.read_text(encoding="utf-8"))
    rejected_values = report.get("rejected_indices")
    if not isinstance(rejected_values, list):
        raise ValueError("K6 bipartite report has no rejected_indices")
    rejected = set(map(int, rejected_values))
    if len(rejected) != len(rejected_values) or not rejected <= k6_universe:
        raise ValueError("invalid K6 bipartite rejection set")
    if report.get("input_graphs") != len(k6_universe):
        raise ValueError("K6 bipartite report does not cover the K6 base")
    before = {index for index, population in selected.items() if population == "K6"}
    removed = before & rejected
    for index in removed:
        selected.pop(index)
    return {
        "kind": "K6_bipartite_rank_rejections",
        "path": str(path),
        "sha256": actual,
        "before": len(before),
        "removed": len(removed),
        "after": len(before - removed),
        "removed_indices_sha256": stable_hash(sorted(removed)),
    }


def exact_residue_indices(
    k7_decisions: Path,
    k6_report: Path,
    k7_strict_h_report: Path | None,
    k7_support_report: Path | None,
    k6_bipartite_report: Path | None,
) -> tuple[dict[int, str], list[dict]]:
    k7, k7_provenance = load_k7_survivors(k7_decisions)
    k6, k6_provenance = load_k6_survivors(k6_report)
    if k7 & k6:
        raise ValueError("K6 and K7 base survivor strata overlap")
    selected = {index: "K7" for index in k7}
    selected.update((index, "K6") for index in k6)
    provenance = [k7_provenance, k6_provenance]
    if k7_strict_h_report is not None:
        provenance.append(apply_k7_strict_h(selected, k7, k7_strict_h_report))
    if k7_support_report is not None:
        provenance.append(apply_k7_support_report(selected, k7, k7_support_report))
    if k6_bipartite_report is not None:
        provenance.append(
            apply_k6_bipartite_report(selected, k6, k6_bipartite_report)
        )
    return selected, provenance


def validate_adjacency(adj: Sequence[int]) -> None:
    if len(adj) != N:
        raise ValueError("wrong adjacency order")
    full = (1 << N) - 1
    for vertex, row in enumerate(adj):
        if row & ~full or row & (1 << vertex):
            raise ValueError("invalid adjacency mask")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(adj[other] & (1 << vertex)):
                raise ValueError("asymmetric adjacency")


def load_graphs(corpus: Path, selected: dict[int, str]) -> tuple[list[dict], str]:
    graphs: list[dict] = []
    remaining = set(selected)
    digest = hashlib.sha256()
    lines = 0
    with corpus.open("rb") as stream:
        for index, raw in enumerate(stream):
            digest.update(raw)
            lines += 1
            if index not in remaining:
                continue
            fields = raw.split()
            if len(fields) != N + 1 or fields[0] != str(N).encode("ascii"):
                raise ValueError(f"malformed corpus line {index + 1}")
            adjacency = tuple(map(int, fields[1:]))
            validate_adjacency(adjacency)
            graphs.append({
                "index": index,
                "population": selected[index],
                "adjacency": adjacency,
            })
            remaining.remove(index)
    actual = digest.hexdigest()
    if lines != EXPECTED_CORPUS_GRAPHS or actual != EXPECTED_CORPUS_SHA256:
        raise ValueError(f"unexpected corpus lines/hash: {lines}, {actual}")
    if remaining:
        raise ValueError(f"corpus omitted {len(remaining)} selected indices")
    return graphs, actual


def ncircle(adj: Sequence[int], seed: Sequence[int], order: Sequence[int]) -> int:
    placed = set(seed)
    answer = 0
    for vertex in order:
        neighbours = sum(bool(adj[vertex] & (1 << other)) for other in placed)
        if neighbours == 5:
            answer += 1
        placed.add(vertex)
    return answer


@dataclass(frozen=True)
class Result:
    ordinal: int
    index: int
    population: str
    status: str
    winning_order: int
    winning_circles: int
    orders_available: int
    orders_tried: int
    kernel_calls: int
    kernel_killed: int
    kernel_survivors: int
    kernel_aborts: int
    total_nodes: int
    unresolved_cells: int
    elapsed_seconds: float
    error_type: str | None = None
    error_message: str | None = None
    error_traceback: str | None = None
    winning_seed: tuple[int, ...] | None = None
    winning_placement_order: tuple[int, ...] | None = None
    winning_slice_records: tuple[dict, ...] | None = None


def analyze_graph(
    graph: dict,
    max_orders: int,
    max_nodes: int,
    zero_circle_only: bool,
    include_bulk_order: bool,
    slices: int,
) -> Result:
    """Return KILLED iff one order is KILLED on every theta slice."""

    started = time.perf_counter()
    adj = graph["adjacency"]
    orders = cdriver6.gen_orders(adj, N, kmax=max_orders)
    indexed = [
        (position, seed, order, ncircle(adj, seed, order))
        for position, (seed, order) in enumerate(orders)
    ]
    candidates = [entry for entry in indexed if include_bulk_order or entry[0] != 0]
    if zero_circle_only:
        candidates = [entry for entry in candidates if entry[3] == 0]
    counts: Counter[str] = Counter()
    total_nodes = 0
    unresolved_cells = 0
    for tried, (position, seed, order, circles) in enumerate(candidates, 1):
        order_killed = True
        slice_records = []
        for part in range(slices):
            lo = 2.0 * math.pi * part / slices
            hi = 2.0 * math.pi * (part + 1) / slices
            status, nodes, cells = cdriver6.decide6(
                adj,
                N,
                seed=seed,
                order=order,
                th0=(lo, hi),
                max_nodes=max_nodes,
            )
            if status not in KERNEL_STATUSES:
                raise RuntimeError(f"unknown kernel status {status!r}")
            if nodes < 0 or cells < 0:
                raise RuntimeError("kernel returned a negative work counter")
            counts[status] += 1
            total_nodes += nodes
            unresolved_cells += cells
            slice_records.append({
                "part": part,
                "lo": lo,
                "hi": hi,
                "status": status,
                "nodes": nodes,
                "unresolved_cells": cells,
            })
            if status != "KILLED":
                order_killed = False
                break
        if order_killed:
            return Result(
                graph["ordinal"], graph["index"], graph["population"],
                "KILLED", position, circles, len(orders), tried,
                sum(counts.values()), counts["KILLED"], counts["SURVIVORS"],
                counts["ABORT"], total_nodes, unresolved_cells,
                time.perf_counter() - started,
                winning_seed=tuple(seed),
                winning_placement_order=tuple(order),
                winning_slice_records=tuple(slice_records),
            )
    final_status = "ABORT" if counts["ABORT"] else "UNRESOLVED"
    return Result(
        graph["ordinal"], graph["index"], graph["population"], final_status,
        -1, -1, len(orders), len(candidates), sum(counts.values()),
        counts["KILLED"], counts["SURVIVORS"], counts["ABORT"], total_nodes,
        unresolved_cells, time.perf_counter() - started,
    )


def analyze_graph_safe(graph: dict, parameters: dict) -> Result:
    started = time.perf_counter()
    try:
        return analyze_graph(
            graph,
            parameters["orders"],
            parameters["cap"],
            parameters["zero_circle_only"],
            parameters["include_bulk_order"],
            parameters["slices"],
        )
    except Exception as error:
        return Result(
            graph["ordinal"], graph["index"], graph["population"],
            "INFRA_ERROR", -1, -1, 0, 0, 0, 0, 0, 0, 0, 0,
            time.perf_counter() - started,
            type(error).__name__, str(error), traceback.format_exc(),
        )


def initialize_native_thread() -> None:
    if sys.platform != "darwin":
        return
    libc = ctypes.CDLL(None)
    setter = libc.pthread_set_qos_class_self_np
    setter.argtypes = (ctypes.c_uint, ctypes.c_int)
    setter.restype = ctypes.c_int
    error = int(setter(QOS_CLASS_USER_INITIATED, 0))
    if error:
        raise OSError(error, "pthread_set_qos_class_self_np failed")


def control_graph(negative: bool) -> dict:
    """Return a K8 negative control or a circle-family positive control."""

    adjacency = [0] * N
    clique_order = 8 if negative else 7
    for vertex in range(clique_order):
        for other in range(vertex):
            adjacency[vertex] |= 1 << other
            adjacency[other] |= 1 << vertex
    # In the positive control, all twelve outside points can be chosen as
    # distinct generic points on the common-unit circle of seed vertices 0..4.
    # The same edges also make every vertex placeable by the interval driver.
    for vertex in range(clique_order, N):
        for other in range(5):
            adjacency[vertex] |= 1 << other
            adjacency[other] |= 1 << vertex
    validate_adjacency(adjacency)
    return {
        "ordinal": 0,
        "index": -2 if negative else -1,
        "population": "NEGATIVE_CONTROL" if negative else "POSITIVE_CONTROL",
        "adjacency": tuple(adjacency),
    }


def control_slices(graph: dict, slices: int) -> dict:
    """Run one fixed eligible order over the exact production slice tiling."""

    orders = cdriver6.gen_orders(graph["adjacency"], N, kmax=12)
    if not orders:
        raise RuntimeError("kernel control has no eligible placement order")
    seed, order = orders[0]
    records = []
    for part in range(slices):
        lo = 2.0 * math.pi * part / slices
        hi = 2.0 * math.pi * (part + 1) / slices
        status, nodes, cells = cdriver6.decide6(
            graph["adjacency"],
            N,
            seed=seed,
            order=order,
            th0=(lo, hi),
            max_nodes=1_000,
        )
        if status not in KERNEL_STATUSES:
            raise RuntimeError(f"unknown control kernel status {status!r}")
        records.append({
            "part": part,
            "lo": lo,
            "hi": hi,
            "status": status,
            "nodes": nodes,
            "unresolved_cells": cells,
        })
    return {
        "seed": list(seed),
        "placement_order": list(order),
        "slices": slices,
        "status_counts": dict(Counter(record["status"] for record in records)),
        "records": records,
    }


def run_kernel_controls(slices: int) -> dict:
    positive = control_slices(control_graph(False), slices)
    negative = control_slices(control_graph(True), slices)
    if any(record["status"] == "KILLED" for record in positive["records"]):
        raise RuntimeError("interval kernel falsely killed a positive-control slice")
    if any(record["status"] != "KILLED" for record in negative["records"]):
        raise RuntimeError("interval kernel did not kill every K8 control slice")
    return {
        "status": "PASS",
        "production_slices_exercised": slices,
        "positive": {
            "description": (
                "K7 plus twelve distinct generic points on the common-unit "
                "circle of five seed vertices"
            ),
            "required": "every production slice is not KILLED",
            **positive,
        },
        "negative": {
            "description": "K8 subgraph, impossible in R6",
            "required": "every production slice is KILLED",
            **negative,
        },
    }


def result_path(directory: Path, graph: dict) -> Path:
    return directory / "results" / (
        f"result_{graph['ordinal']:06d}_{graph['index']:07d}.json"
    )


def validate_checkpoint(
    value: dict,
    config_hash: str,
    graph: dict,
    expected_slices: int,
) -> dict:
    if value.get("schema") != 1 or value.get("config_sha256") != config_hash:
        raise ValueError("result checkpoint schema/config mismatch")
    result = value.get("result")
    if not isinstance(result, dict):
        raise ValueError("result checkpoint has no result")
    if (
        result.get("ordinal") != graph["ordinal"]
        or result.get("index") != graph["index"]
        or result.get("population") != graph["population"]
        or result.get("status") not in RESULT_STATUSES
    ):
        raise ValueError("result checkpoint identity/status mismatch")
    integer_fields = (
        "winning_order", "winning_circles", "orders_available", "orders_tried",
        "kernel_calls", "kernel_killed", "kernel_survivors", "kernel_aborts",
        "total_nodes", "unresolved_cells",
    )
    if any(not isinstance(result.get(name), int) for name in integer_fields):
        raise ValueError("result checkpoint has a noninteger work field")
    nonnegative = integer_fields[2:]
    if any(result[name] < 0 for name in nonnegative):
        raise ValueError("result checkpoint has a negative work field")
    if result["kernel_calls"] != sum(
        result[name]
        for name in ("kernel_killed", "kernel_survivors", "kernel_aborts")
    ):
        raise ValueError("result checkpoint kernel counts do not add up")
    status = result["status"]
    slices = value.get("search_slices")
    if slices != expected_slices:
        raise ValueError("result checkpoint slice count/config mismatch")
    if status == "KILLED" and (
        result["winning_order"] < 0 or result["kernel_killed"] < slices
    ):
        raise ValueError("invalid KILLED checkpoint witness counters")
    if status == "KILLED":
        seed = result.get("winning_seed")
        order = result.get("winning_placement_order")
        records = result.get("winning_slice_records")
        if (
            not isinstance(seed, (list, tuple))
            or len(seed) not in (6, 7)
            or not isinstance(order, (list, tuple))
            or sorted(list(seed) + list(order)) != list(range(N))
            or not isinstance(records, (list, tuple))
            or len(records) != slices
        ):
            raise ValueError("KILLED checkpoint omits its exact winning order")
        for part, record in enumerate(records):
            if (
                not isinstance(record, dict)
                or record.get("part") != part
                or record.get("status") != "KILLED"
                or not isinstance(record.get("nodes"), int)
                or record["nodes"] < 0
                or not isinstance(record.get("unresolved_cells"), int)
                or record["unresolved_cells"] < 0
            ):
                raise ValueError("invalid KILLED per-slice checkpoint witness")
            expected_lo = 2.0 * math.pi * part / slices
            expected_hi = 2.0 * math.pi * (part + 1) / slices
            if record.get("lo") != expected_lo or record.get("hi") != expected_hi:
                raise ValueError("KILLED checkpoint has a wrong slice interval")
    if status == "ABORT" and result["kernel_aborts"] == 0:
        raise ValueError("ABORT checkpoint has no kernel abort")
    if status == "UNRESOLVED" and result["kernel_aborts"] != 0:
        raise ValueError("UNRESOLVED checkpoint contains a kernel abort")
    if status == "INFRA_ERROR" and not result.get("error_type"):
        raise ValueError("INFRA_ERROR checkpoint has no exception type")
    return result


def status_summary(results: Iterable[dict], total: int) -> dict:
    materialized = list(results)
    counts = Counter(result["status"] for result in materialized)
    return {
        "completed_graphs": len(materialized),
        "total_graphs": total,
        "status_counts": {status: counts[status] for status in RESULT_STATUSES},
        "certified_killed": counts["KILLED"],
        "unresolved_total": len(materialized) - counts["KILLED"],
    }


def decision_text(results: Sequence[dict]) -> str:
    fields = (
        "ordinal", "index", "population", "status", "winning_order",
        "winning_circles", "orders_available", "orders_tried", "kernel_calls",
        "kernel_killed", "kernel_survivors", "kernel_aborts", "total_nodes",
        "unresolved_cells", "elapsed_seconds", "error_type", "error_message",
    )
    lines = ["\t".join(fields) + "\n"]
    for result in results:
        values = []
        for field in fields:
            value = result.get(field)
            if field == "elapsed_seconds":
                value = f"{float(value):.9f}"
            if value is None:
                value = ""
            values.append(str(value).replace("\t", " ").replace("\n", " "))
        lines.append("\t".join(values) + "\n")
    return "".join(lines)


@contextmanager
def exclusive_campaign_lock(directory: Path):
    """Prevent two processes from mutating one checkpoint directory."""

    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "campaign.lock"
    with path.open("a+", encoding="ascii") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                f"checkpoint directory is locked by another process: {directory}"
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


def _run_campaign_unlocked(
    graphs: list[dict],
    *,
    configuration: dict,
    workers: int,
    checkpoint_every: int,
    progress_every: int,
    retry_infra_errors: bool,
    checkpoint_dir: Path,
    report_path: Path,
    decisions_path: Path,
) -> tuple[dict, bool]:
    config_hash = stable_hash(configuration)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    campaign_path = checkpoint_dir / "campaign.json"
    campaign = {
        "schema": 1,
        "configuration": configuration,
        "config_sha256": config_hash,
    }
    if campaign_path.exists():
        existing = json.loads(campaign_path.read_text(encoding="utf-8"))
        if existing.get("configuration") != configuration:
            raise RuntimeError(
                f"checkpoint directory belongs to another campaign: {checkpoint_dir}"
            )
    else:
        campaign.update(
            created_utc=utc_now(),
            creator_pid=os.getpid(),
            argv=[sys.executable, *sys.argv],
            command=shlex.join([sys.executable, *sys.argv]),
        )
        atomic_json(campaign_path, campaign)

    session = {
        "schema": 1,
        "config_sha256": config_hash,
        "pid": os.getpid(),
        "started_utc": utc_now(),
        "argv": [sys.executable, *sys.argv],
        "command": shlex.join([sys.executable, *sys.argv]),
        "workers": workers,
    }
    session_path = checkpoint_dir / "sessions" / (
        f"session_{time.time_ns()}_{os.getpid()}.json"
    )
    atomic_json(session_path, session)

    completed: dict[int, dict] = {}
    pending: deque[dict] = deque()
    for graph in graphs:
        path = result_path(checkpoint_dir, graph)
        if not path.exists():
            pending.append(graph)
            continue
        value = json.loads(path.read_text(encoding="utf-8"))
        result = validate_checkpoint(
            value, config_hash, graph, configuration["search"]["slices"]
        )
        if retry_infra_errors and result["status"] == "INFRA_ERROR":
            atomic_json(
                checkpoint_dir / "failures" / (
                    f"retry_{graph['ordinal']:06d}_{graph['index']:07d}_"
                    f"{time.time_ns()}.json"
                ),
                value,
            )
            pending.append(graph)
        else:
            completed[graph["ordinal"]] = result

    started = time.perf_counter()
    newly_completed = 0
    last_print = len(completed)

    def ordered_results() -> list[dict]:
        return [completed[position] for position in sorted(completed)]

    def write_progress(force: bool = False) -> None:
        nonlocal last_print
        summary = status_summary(completed.values(), len(graphs))
        atomic_json(
            checkpoint_dir / "progress.json",
            {
                "schema": 1,
                "config_sha256": config_hash,
                "session": session_path.name,
                "pid": os.getpid(),
                "updated_utc": utc_now(),
                **summary,
            },
        )
        done = summary["completed_graphs"]
        if force or done == len(graphs) or done - last_print >= progress_every:
            counts = summary["status_counts"]
            print(
                f"progress {done}/{len(graphs)}; KILLED {counts['KILLED']}; "
                f"ABORT {counts['ABORT']}; UNRESOLVED {counts['UNRESOLVED']}; "
                f"INFRA_ERROR {counts['INFRA_ERROR']}",
                file=sys.stderr,
                flush=True,
            )
            last_print = done

    def save_result(graph: dict, result: dict) -> None:
        nonlocal newly_completed
        value = {
            "schema": 1,
            "config_sha256": config_hash,
            "search_slices": configuration["search"]["slices"],
            "completed_utc": utc_now(),
            "result": result,
        }
        validate_checkpoint(
            value, config_hash, graph, configuration["search"]["slices"]
        )
        if result["status"] == "INFRA_ERROR":
            atomic_json(
                checkpoint_dir / "failures" / (
                    f"infra_{graph['ordinal']:06d}_{graph['index']:07d}_"
                    f"{time.time_ns()}.json"
                ),
                value,
            )
        atomic_json(result_path(checkpoint_dir, graph), value)
        completed[graph["ordinal"]] = result
        newly_completed += 1
        if newly_completed % checkpoint_every == 0:
            write_progress()

    atomic_json(
        checkpoint_dir / "run_state.json",
        {
            "schema": 1,
            "config_sha256": config_hash,
            "status": "RUNNING",
            "session": session_path.name,
            "pid": os.getpid(),
            "updated_utc": utc_now(),
        },
    )
    write_progress(force=True)
    executor = ThreadPoolExecutor(
        max_workers=workers, initializer=initialize_native_thread
    )
    active: dict = {}
    try:
        while pending or active:
            while pending and len(active) < 2 * workers:
                graph = pending.popleft()
                try:
                    future = executor.submit(
                        analyze_graph_safe, graph, configuration["search"]
                    )
                except Exception as error:
                    save_result(graph, asdict(Result(
                        graph["ordinal"], graph["index"], graph["population"],
                        "INFRA_ERROR", -1, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                        0.0, type(error).__name__, str(error),
                        traceback.format_exc(),
                    )))
                    continue
                active[future] = graph
            if not active:
                continue
            finished, _ = wait(active, return_when=FIRST_COMPLETED)
            for future in finished:
                graph = active.pop(future)
                try:
                    result = asdict(future.result())
                except Exception as error:
                    result = asdict(Result(
                        graph["ordinal"], graph["index"], graph["population"],
                        "INFRA_ERROR", -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0.0,
                        type(error).__name__, str(error), traceback.format_exc(),
                    ))
                save_result(graph, result)
        executor.shutdown(wait=True)
    except BaseException:
        for future in active:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)
        atomic_json(
            checkpoint_dir / "run_state.json",
            {
                "schema": 1,
                "config_sha256": config_hash,
                "status": "INTERRUPTED",
                "session": session_path.name,
                "pid": os.getpid(),
                "completed_graphs": len(completed),
                "updated_utc": utc_now(),
            },
        )
        write_progress(force=True)
        raise

    write_progress(force=True)
    results = ordered_results()
    if len(results) != len(graphs):
        raise AssertionError("campaign ended without every graph result")
    if [result["ordinal"] for result in results] != list(range(len(graphs))):
        raise AssertionError("final results are not in deterministic order")
    atomic_text(decisions_path, decision_text(results))
    summary = status_summary(results, len(graphs))
    checkpoint_entries = [
        [result_path(checkpoint_dir, graph).name,
         sha256(result_path(checkpoint_dir, graph))]
        for graph in graphs
    ]
    report = {
        "schema": 2,
        "claim": (
            "Only KILLED is a certified non-realizability result. ABORT, "
            "UNRESOLVED, and INFRA_ERROR make no mathematical claim."
        ),
        "KILLED_semantics": (
            "There exists a tried placement order for which every theta "
            "slice returned KILLED from the outward-rounded interval kernel."
        ),
        "configuration": configuration,
        "config_sha256": config_hash,
        "checkpoint_directory": str(checkpoint_dir),
        "checkpoint_result_index_sha256": stable_hash(checkpoint_entries),
        "session": session,
        **summary,
        "killed_by_population": {
            population: sum(
                result["status"] == "KILLED"
                and result["population"] == population
                for result in results
            )
            for population in ("K7", "K6")
        },
        "kernel_status_totals": {
            "KILLED": sum(result["kernel_killed"] for result in results),
            "SURVIVORS": sum(result["kernel_survivors"] for result in results),
            "ABORT": sum(result["kernel_aborts"] for result in results),
        },
        "runtime": {
            "resume_wall_seconds": time.perf_counter() - started,
            "newly_completed_graphs": newly_completed,
            "sum_graph_seconds": sum(result["elapsed_seconds"] for result in results),
            "maximum_graph_seconds": max(
                (result["elapsed_seconds"] for result in results), default=0.0
            ),
        },
        "decisions": {
            "path": str(decisions_path),
            "sha256": sha256(decisions_path),
            "rows_excluding_header": len(results),
        },
        "results": results,
    }
    atomic_json(report_path, report)
    has_infra = summary["status_counts"]["INFRA_ERROR"] != 0
    final_state = "COMPLETE_WITH_INFRA_ERRORS" if has_infra else "COMPLETE"
    atomic_json(
        checkpoint_dir / "run_state.json",
        {
            "schema": 1,
            "config_sha256": config_hash,
            "status": final_state,
            "session": session_path.name,
            "pid": os.getpid(),
            "completed_graphs": len(results),
            "status_counts": summary["status_counts"],
            "report": str(report_path),
            "report_sha256": sha256(report_path),
            "updated_utc": utc_now(),
        },
    )
    session.update(
        completed_utc=utc_now(),
        final_status=final_state,
        status_counts=summary["status_counts"],
        report=str(report_path),
    )
    atomic_json(session_path, session)
    return report, has_infra


def run_campaign(
    graphs: list[dict],
    *,
    configuration: dict,
    workers: int,
    checkpoint_every: int,
    progress_every: int,
    retry_infra_errors: bool,
    checkpoint_dir: Path,
    report_path: Path,
    decisions_path: Path,
) -> tuple[dict, bool]:
    with exclusive_campaign_lock(checkpoint_dir):
        return _run_campaign_unlocked(
            graphs,
            configuration=configuration,
            workers=workers,
            checkpoint_every=checkpoint_every,
            progress_every=progress_every,
            retry_infra_errors=retry_infra_errors,
            checkpoint_dir=checkpoint_dir,
            report_path=report_path,
            decisions_path=decisions_path,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("aeq_d6_n19.txt"))
    parser.add_argument(
        "--k7-decisions", type=Path,
        default=Path("d6_k7_rank_python_full_decisions.tsv.gz"),
    )
    parser.add_argument(
        "--k6-report", type=Path, default=Path("d6_k6_support_report.json")
    )
    parser.add_argument("--k7-strict-h-report", type=Path)
    parser.add_argument("--k7-support-report", type=Path)
    parser.add_argument("--k6-bipartite-report", type=Path)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--orders", type=int, default=12)
    parser.add_argument("--cap", type=int, default=1_000_000)
    parser.add_argument(
        "--slices", type=int, default=1,
        help="tile the first circle into this many certified intervals",
    )
    parser.add_argument("--sample", type=int)
    parser.add_argument("--sample-seed", type=int, default=6)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--include-circle-orders", action="store_true")
    parser.add_argument("--include-bulk-order", action="store_true")
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=Path(".runs/d6_interval_residue_checkpoints"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path(".runs/d6_interval_residue_report.json"),
    )
    parser.add_argument(
        "--decisions", type=Path,
        default=Path(".runs/d6_interval_residue_decisions.tsv"),
    )
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=50)
    parser.add_argument("--retry-infra-errors", action="store_true")
    parser.add_argument(
        "--outer-launch-command",
        help=(
            "exact user-supplied outer command, including wrappers such as "
            "caffeinate/taskpolicy, recorded verbatim in campaign provenance"
        ),
    )
    parser.add_argument(
        "--selection-only", action="store_true",
        help="validate and print the selected population without kernel calls",
    )
    args = parser.parse_args()
    if (
        args.workers < 1 or args.orders < 1 or args.cap < 1 or args.slices < 1
        or args.checkpoint_every < 1 or args.progress_every < 1
    ):
        parser.error("workers, orders, cap, slices, and checkpoint periods must be positive")
    if args.shards < 1 or not (0 <= args.shard < args.shards):
        parser.error("require shards >= 1 and 0 <= shard < shards")
    if args.sample is not None and args.sample < 1:
        parser.error("sample must be positive")

    selected, selection_provenance = exact_residue_indices(
        args.k7_decisions,
        args.k6_report,
        args.k7_strict_h_report,
        args.k7_support_report,
        args.k6_bipartite_report,
    )
    graphs, corpus_hash = load_graphs(args.corpus, selected)
    graphs = [
        graph for position, graph in enumerate(graphs)
        if position % args.shards == args.shard
    ]
    if args.sample is not None and args.sample < len(graphs):
        random.Random(args.sample_seed).shuffle(graphs)
        graphs = graphs[:args.sample]
    graphs.sort(key=lambda graph: graph["index"])
    for ordinal, graph in enumerate(graphs):
        graph["ordinal"] = ordinal

    selected_counts = Counter(graph["population"] for graph in graphs)
    selected_summary = {
        "graphs": len(graphs),
        "population_counts": dict(selected_counts),
        "indices_sha256": stable_hash([graph["index"] for graph in graphs]),
        "selection_layers": selection_provenance,
    }
    if args.selection_only:
        print(json.dumps(selected_summary, indent=2, sort_keys=True))
        return

    with _KERNEL_LOCK:
        cdriver6._kernel()  # type: ignore[attr-defined]
    kernel_binary = Path(cdriver6.HERE) / cdriver6.LIBNAME
    controls = run_kernel_controls(args.slices)
    try:
        compiler = subprocess.run(
            ["cc", "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        ).stdout.splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError) as error:
        compiler = f"unavailable: {type(error).__name__}: {error}"
    source_hashes = {
        "runner": sha256(Path(__file__).resolve()),
        "cdriver6.py": sha256(ROOT / "cdriver6.py"),
        "ckernel6.c": sha256(ROOT / "ckernel6.c"),
        "ival.py": sha256(ROOT / "ival.py"),
        cdriver6.LIBNAME: sha256(kernel_binary),
    }
    search = {
        "orders": args.orders,
        "cap": args.cap,
        "slices": args.slices,
        "zero_circle_only": not args.include_circle_orders,
        "include_bulk_order": args.include_bulk_order,
    }
    configuration = {
        "schema": 2,
        "sources": source_hashes,
        "inputs": {
            "corpus": str(args.corpus),
            "corpus_sha256": corpus_hash,
        },
        "selection": selected_summary,
        "search": search,
        "sampling": {
            "sample": args.sample,
            "sample_seed": args.sample_seed,
            "shards": args.shards,
            "shard": args.shard,
        },
        "trust_assumptions": {
            "basic_arithmetic": (
                "IEEE-754 binary64 basic operations and sqrt are correctly rounded; "
                "each interval endpoint is expanded with nextafter."
            ),
            "transcendentals": (
                "macOS libm cos endpoint values are assumed within 8 ulps; the "
                "kernel pads both directions by 8 nextafter steps and includes all "
                "interior extrema."
            ),
            "candidate_nonedges": "unconstrained and allowed to be unit distance",
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "compiler": compiler,
            "kernel_compile_command": (
                f"cc -O2 -shared -o {cdriver6.LIBNAME} ckernel6.c -lm"
            ),
        },
        "git": git_provenance(),
        "launch": {
            "python_argv": [sys.executable, *sys.argv],
            "python_command": shlex.join([sys.executable, *sys.argv]),
            "outer_launch_command": args.outer_launch_command,
            "outer_launch_command_user_supplied": (
                args.outer_launch_command is not None
            ),
        },
        "kernel_controls": controls,
    }
    print(
        f"interval residue: {len(graphs)} graphs, {args.workers} native threads, "
        f"up to {args.orders} orders, cap {args.cap}; checkpoint "
        f"{args.checkpoint_dir}",
        file=sys.stderr,
        flush=True,
    )
    report, has_infra = run_campaign(
        graphs,
        configuration=configuration,
        workers=args.workers,
        checkpoint_every=args.checkpoint_every,
        progress_every=args.progress_every,
        retry_infra_errors=args.retry_infra_errors,
        checkpoint_dir=args.checkpoint_dir,
        report_path=args.output,
        decisions_path=args.decisions,
    )
    counts = report["status_counts"]
    print(
        f"complete: KILLED {counts['KILLED']}; ABORT {counts['ABORT']}; "
        f"UNRESOLVED {counts['UNRESOLVED']}; INFRA_ERROR "
        f"{counts['INFRA_ERROR']}; wrote {args.output}",
        file=sys.stderr,
        flush=True,
    )
    if has_infra:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
