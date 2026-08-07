#!/usr/bin/env python3
"""Extract the hash-pinned pattern-954 algebraic core and current targets.

The source graph is line/index 954 (zero based) of ``aeq_d6_n14.txt``.
The certified obstruction discards source vertex 1 and two unneeded required
edges, leaving a 13-vertex, 54-edge unit-edge graph.  Candidate nonedges are
not constraints.

The optional target TSV is the exact 12,839-graph K7 residue after the full
degree-one positive-polynomial pass.  It is a transient containment input;
the tracked JSON records its content hash and the complete selection
provenance.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

import cdriver6


ROOT = Path(__file__).resolve().parent
SOURCE_INDEX = 954
SOURCE_CORPUS_SHA256 = (
    "0e3d74c081b272731848d655da0c68cfba09435e39a2fd2107c32c2bd3b378f0"
)
SOURCE_CORPUS_GRAPHS = 1_052
EXPECTED_SOURCE_ADJACENCY = (
    5694, 365, 13555, 6995, 16077, 10119, 14750,
    16244, 10986, 15801, 15029, 14296, 11997, 8180,
)
RANK_SURVIVORS_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
PARENT_SELECTION_SHA256 = (
    "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0"
)
FULL_DUAL_REPORT_SHA256 = (
    "c2de7e06bbe4f3d163a4f747d42721f7f60da985b867ce74669664431738e345"
)
FULL_DUAL_DECISIONS_SHA256 = (
    "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
)
CURRENT_SELECTION_SIZE = 12_839
CURRENT_SELECTION_INDICES_SHA256 = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)

# Seed order defines defect-coordinate positions 0,...,6.
SEED = (4, 7, 9, 10, 11, 12, 13)
ROLES = {
    "A": 0,
    "C": 5,
    "E": 2,
    "B": 3,
    "D": 8,
    "F": 6,
}
DEFECT_BOUNDS = {
    "A": (1, 4, 6),
    "C": (0, 4, 5),
    "E": (2, 4),
    "B": (1, 3, 6),
    "D": (0, 3, 5),
    "F": (2, 3),
}
TRIANGLES = (("A", "C", "E"), ("B", "D", "F"))
BRIDGE = ("E", "F")
UNUSED_SOURCE_VERTEX = 1
UNUSED_SOURCE_EDGES = ((0, 3), (5, 8))


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


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="ascii", newline="") as stream:
        stream.write(value)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def validate_graph(adjacency: Sequence[int]) -> None:
    count = len(adjacency)
    full = (1 << count) - 1
    for vertex, row in enumerate(adjacency):
        if not isinstance(row, int) or isinstance(row, bool):
            raise ValueError("noninteger adjacency")
        if row & ~full or row & (1 << vertex):
            raise ValueError("invalid adjacency bit or loop")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(
                adjacency[other] & (1 << vertex)
            ):
                raise ValueError("asymmetric adjacency")


def bits(mask: int) -> Iterator[int]:
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def source_pattern(path: Path) -> tuple[int, ...]:
    if sha256(path) != SOURCE_CORPUS_SHA256:
        raise ValueError("n=14 corpus hash mismatch")
    rows = path.read_text(encoding="ascii").splitlines()
    if len(rows) != SOURCE_CORPUS_GRAPHS:
        raise ValueError("n=14 corpus population mismatch")
    fields = tuple(map(int, rows[SOURCE_INDEX].split()))
    if fields[0] != 14 or len(fields) != 15:
        raise ValueError("pattern-954 source row is malformed")
    adjacency = fields[1:]
    validate_graph(adjacency)
    if adjacency != EXPECTED_SOURCE_ADJACENCY:
        raise ValueError("pattern-954 adjacency changed")
    return adjacency


def core_graph(source: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    original_vertices = tuple(
        vertex for vertex in range(14) if vertex != UNUSED_SOURCE_VERTEX
    )
    position = {vertex: index for index, vertex in enumerate(original_vertices)}
    removed = {frozenset(edge) for edge in UNUSED_SOURCE_EDGES}
    adjacency = [0] * len(original_vertices)
    for first, second in combinations(original_vertices, 2):
        if not (source[first] & (1 << second)):
            continue
        if frozenset((first, second)) in removed:
            continue
        left, right = position[first], position[second]
        adjacency[left] |= 1 << right
        adjacency[right] |= 1 << left
    validate_graph(adjacency)
    return original_vertices, tuple(adjacency)


def defect_mask(adjacency: Sequence[int], seed: Sequence[int], vertex: int) -> int:
    return sum(
        1 << coordinate
        for coordinate, seed_vertex in enumerate(seed)
        if not (adjacency[vertex] & (1 << seed_vertex))
    )


def circle_stages(
    adjacency: Sequence[int], seed: Sequence[int], order: Sequence[int]
) -> list[dict]:
    placed = set(seed)
    answer = []
    for position, vertex in enumerate(order):
        neighbours = sum(
            bool(adjacency[vertex] & (1 << other)) for other in placed
        )
        if neighbours == 5:
            answer.append({"position": position, "vertex": vertex})
        placed.add(vertex)
    return answer


def current_selection(
    parent_selection_path: Path,
    dual_report_path: Path,
    decisions_path: Path,
) -> tuple[list[int], dict]:
    if sha256(parent_selection_path) != PARENT_SELECTION_SHA256:
        raise ValueError("parent K7 selection hash mismatch")
    if sha256(dual_report_path) != FULL_DUAL_REPORT_SHA256:
        raise ValueError("full-dual report hash mismatch")
    if sha256(decisions_path) != FULL_DUAL_DECISIONS_SHA256:
        raise ValueError("full-dual decision archive hash mismatch")
    parent = json.loads(parent_selection_path.read_text(encoding="utf-8"))
    parent_indices = parent.get("selected_indices")
    if (
        not isinstance(parent_indices, list)
        or parent.get("selected") != 12_941
        or len(parent_indices) != 12_941
        or len(set(parent_indices)) != len(parent_indices)
    ):
        raise ValueError("parent K7 selection is malformed")
    report = json.loads(dual_report_path.read_text(encoding="utf-8"))
    summary = report.get("summary", {})
    if (
        summary.get("graphs") != 12_941
        or summary.get("complete") != 12_941
        or summary.get("marginal_dual_rejected") != 102
        or summary.get("survivors") != CURRENT_SELECTION_SIZE
        or summary.get("infra_errors") != 0
    ):
        raise ValueError("full-dual report accounting mismatch")
    reported_rejected = summary.get("marginal_dual_rejected_indices")
    if not isinstance(reported_rejected, list) or len(set(reported_rejected)) != 102:
        raise ValueError("full-dual report rejection list is malformed")

    with gzip.open(decisions_path, "rt", encoding="ascii", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != len(parent_indices):
        raise ValueError("full-dual decisions do not cover the parent selection")
    observed = [int(row["index"]) for row in rows]
    if observed != parent_indices:
        raise ValueError("full-dual decision order differs from parent selection")
    if any(row["status"] not in ("REJECTED", "SURVIVOR") for row in rows):
        raise ValueError("full-dual decision archive contains a bad status")
    rejected = [int(row["index"]) for row in rows if row["status"] == "REJECTED"]
    if rejected != reported_rejected:
        raise ValueError("full-dual report/archive rejection lists differ")
    selected = [
        int(row["index"]) for row in rows if row["status"] == "SURVIVOR"
    ]
    if len(selected) != CURRENT_SELECTION_SIZE:
        raise ValueError("wrong current K7 survivor count")
    if stable_hash(selected) != CURRENT_SELECTION_INDICES_SHA256:
        raise ValueError("current K7 survivor index hash mismatch")
    return selected, {
        "population": len(selected),
        "indices_sha256": stable_hash(selected),
        "parent_selection": parent_selection_path.name,
        "parent_selection_sha256": PARENT_SELECTION_SHA256,
        "full_dual_report": dual_report_path.name,
        "full_dual_report_sha256": FULL_DUAL_REPORT_SHA256,
        "full_dual_decisions": decisions_path.name,
        "full_dual_decisions_sha256": FULL_DUAL_DECISIONS_SHA256,
        "exact_removed_by_full_dual": len(rejected),
    }


def target_rows(rank_path: Path, selected: Sequence[int]) -> list[dict]:
    if sha256(rank_path) != RANK_SURVIVORS_SHA256:
        raise ValueError("rank-survivor graph input hash mismatch")
    payload = json.loads(rank_path.read_text(encoding="utf-8"))
    rows = payload.get("graphs")
    if not isinstance(rows, list) or len(rows) != 17_764:
        raise ValueError("rank-survivor graph population mismatch")
    by_index = {int(row["index"]): row for row in rows}
    if len(by_index) != len(rows):
        raise ValueError("rank-survivor graph input repeats an index")
    answer = []
    for index in selected:
        graph = by_index.get(index)
        if graph is None:
            raise ValueError(f"current selection index {index} is absent")
        adjacency = tuple(map(int, graph["adjacency"]))
        if len(adjacency) != 19:
            raise ValueError("current target has wrong order")
        validate_graph(adjacency)
        answer.append({"index": index, "adjacency": adjacency})
    return answer


def build_report(args: argparse.Namespace) -> dict:
    source = source_pattern(args.corpus)
    original_vertices, core = core_graph(source)
    selected, selection_provenance = current_selection(
        args.parent_selection, args.dual_report, args.dual_decisions
    )
    targets = target_rows(args.rank_survivors, selected)
    target_text = "".join(
        f"{row['index']} " + " ".join(map(str, row["adjacency"])) + "\n"
        for row in targets
    )
    atomic_text(args.targets_output, target_text)
    target_hash = sha256(args.targets_output)

    engine_orders = []
    for seed, order in cdriver6.gen_orders(source, 14, kmax=12):
        engine_orders.append({
            "seed": list(seed),
            "placement_order": list(order),
            "circle_stages": circle_stages(source, seed, order),
        })
    source_position = {vertex: index for index, vertex in enumerate(original_vertices)}
    core_seed = [source_position[vertex] for vertex in SEED]
    core_roles = {name: source_position[vertex] for name, vertex in ROLES.items()}
    defects = {
        name: list(bits(defect_mask(core, core_seed, vertex)))
        for name, vertex in core_roles.items()
    }
    if defects != {name: list(value) for name, value in DEFECT_BOUNDS.items()}:
        raise AssertionError("extracted core defect masks changed")

    return {
        "schema": 1,
        "kind": "d6_n14_pattern_954_algebraic_core",
        "claim": (
            "The 13-vertex core is a required-unit-edge obstruction in R6; "
            "source/candidate nonedges are unconstrained."
        ),
        "source_pattern": {
            "corpus": args.corpus.name,
            "corpus_sha256": SOURCE_CORPUS_SHA256,
            "corpus_graphs": SOURCE_CORPUS_GRAPHS,
            "zero_based_index": SOURCE_INDEX,
            "adjacency": list(source),
            "adjacency_sha256": stable_hash(list(source)),
            "edges": sum(row.bit_count() for row in source) // 2,
            "unique_K7_seed": list(SEED),
        },
        "obstruction_core": {
            "order": len(core),
            "edges": sum(row.bit_count() for row in core) // 2,
            "original_vertex_labels": list(original_vertices),
            "adjacency": list(core),
            "adjacency_sha256": stable_hash(list(core)),
            "seed": core_seed,
            "roles": core_roles,
            "defect_coordinate_upper_bounds": defects,
            "required_triangles": [list(triangle) for triangle in TRIANGLES],
            "required_bridge": list(BRIDGE),
            "discarded_source_vertex": UNUSED_SOURCE_VERTEX,
            "discarded_source_edges": [list(edge) for edge in UNUSED_SOURCE_EDGES],
            "nonedges_used_as_distance_constraints": False,
        },
        "algebraic_certificate": {
            "simplex_radicand": 7,
            "first_triangle_shared_coordinate": 4,
            "second_triangle_shared_coordinate": 3,
            "bridge_shared_coordinate": 2,
            "sign_cases": [
                {"sigma": sigma, "tau": tau}
                for sigma in (-1, 1) for tau in (-1, 1)
            ],
            "terminal_equation": "3 + sigma*tau = sqrt(7)*(sigma + tau)",
        },
        "interval_engine_inspection": {
            "engine_usable": bool(engine_orders),
            "orders": engine_orders,
            "minimum_circle_stages": min(
                len(order["circle_stages"]) for order in engine_orders
            ),
            "cdriver6_sha256": sha256(ROOT / "cdriver6.py"),
            "ckernel6_sha256": sha256(ROOT / "ckernel6.c"),
            "ival_sha256": sha256(ROOT / "ival.py"),
            "conclusion": (
                "Two orders exist and each has two circle stages. The exact "
                "algebraic certificate supersedes an interval launch."
            ),
        },
        "current_K7_targets": {
            **selection_provenance,
            "rank_survivors": str(args.rank_survivors.relative_to(ROOT)),
            "rank_survivors_sha256": RANK_SURVIVORS_SHA256,
            "transient_tsv": str(args.targets_output.relative_to(ROOT)),
            "transient_tsv_sha256": target_hash,
            "transient_tsv_bytes": args.targets_output.stat().st_size,
            "row_format": "corpus_index followed by 19 adjacency masks",
        },
        "extractor": {
            "file": Path(__file__).name,
            "sha256": sha256(Path(__file__)),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=ROOT / "aeq_d6_n14.txt")
    parser.add_argument(
        "--rank-survivors",
        type=Path,
        default=ROOT / ".runs/d6_k7_rank_survivors.json",
    )
    parser.add_argument(
        "--parent-selection",
        type=Path,
        default=ROOT / "d6_k7_positive_dual_selection.json",
    )
    parser.add_argument(
        "--dual-report",
        type=Path,
        default=ROOT / "d6_k7_positive_dual_full_report.json",
    )
    parser.add_argument(
        "--dual-decisions",
        type=Path,
        default=ROOT / "d6_k7_positive_dual_full_decisions.tsv.gz",
    )
    parser.add_argument(
        "--targets-output",
        type=Path,
        default=ROOT / ".runs/d6_n14_pattern_954_targets_12839.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_n14_pattern_954_input.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = build_report(args)
    atomic_json(args.output, report)
    print(
        f"wrote {args.output}: pattern {report['source_pattern']['edges']} edges, "
        f"core {report['obstruction_core']['order']} vertices/"
        f"{report['obstruction_core']['edges']} edges, "
        f"targets {report['current_K7_targets']['population']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
