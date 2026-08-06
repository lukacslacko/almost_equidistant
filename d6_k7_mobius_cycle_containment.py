#!/usr/bin/env python3
"""Exact required-edge containment for K7 Mobius-cycle obstructions.

The two patterns are a K7 seed plus one outside role for each edge of a
simple 3- or 4-cycle on seed coordinates.  A role of type ``{i,j}`` has
required unit edges to the other five seed vertices, and consecutive roles
have a required unit edge.  Pattern nonedges impose no condition.

The matcher exploits only this required-edge structure.  It exhaustively
enumerates every target K7, every coordinate triangle/quadrilateral, and
every injective role assignment.  It has no heuristic or node cap.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence


ROOT = Path(__file__).resolve().parent
TARGETS = ROOT / ".runs/d6_n14_pattern_954_targets_12839.tsv"
SELECTION = ROOT / "d6_k7_rankone_tetrad_full_selection.json"
PROOF = ROOT / "d6_k7_mobius_cycle_obstructions.md"
VERIFIER = ROOT / "verify_d6_k7_mobius_cycle_containment.py"
TEST = ROOT / "test_d6_k7_mobius_cycle_containment.py"
TARGETS_SHA256 = (
    "9355854172c324f9d93cc4085020974fb9622a010b2e8fe19f5a954d17a7277d"
)
SELECTION_SHA256 = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
INDICES_SHA256 = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)
TARGET_COUNT = 12_839
TARGET_ORDER = 19


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def validate_graph(adjacency: Sequence[int], order: int) -> tuple[int, ...]:
    if len(adjacency) != order:
        raise ValueError(f"expected order {order}, got {len(adjacency)}")
    full = (1 << order) - 1
    answer = tuple(adjacency)
    for vertex, row in enumerate(answer):
        if type(row) is not int or row & ~full or row & (1 << vertex):
            raise ValueError("invalid adjacency row")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(
                answer[other] & (1 << vertex)
            ):
                raise ValueError("asymmetric adjacency")
    return answer


def add_edge(adjacency: list[int], first: int, second: int) -> None:
    if first == second:
        raise ValueError("loop in pattern")
    adjacency[first] |= 1 << second
    adjacency[second] |= 1 << first


def cycle_pattern(length: int) -> dict:
    if length not in (3, 4):
        raise ValueError("only the forbidden 3- and 4-cycles are frozen")
    order = 7 + length
    adjacency = [0] * order
    for first, second in combinations(range(7), 2):
        add_edge(adjacency, first, second)
    types = [(coordinate, (coordinate + 1) % length)
             for coordinate in range(length)]
    for role, defect_type in enumerate(types, 7):
        for seed in range(7):
            if seed not in defect_type:
                add_edge(adjacency, role, seed)
    for role in range(length):
        add_edge(adjacency, 7 + role, 7 + ((role + 1) % length))
    frozen = validate_graph(adjacency, order)
    return {
        "schema": 1,
        "kind": f"d6_k7_two_defect_{length}_cycle_obstruction",
        "order": order,
        "edges": sum(row.bit_count() for row in frozen) // 2,
        "adjacency": list(frozen),
        "adjacency_sha256": stable_hash(list(frozen)),
        "seed_vertices": list(range(7)),
        "role_vertices": list(range(7, order)),
        "role_defect_upper_bounds": [list(pair) for pair in types],
        "required_role_cycle": [
            [7 + role, 7 + ((role + 1) % length)]
            for role in range(length)
        ],
        "nonedges_used_as_distance_constraints": False,
    }


PATTERNS = {length: cycle_pattern(length) for length in (3, 4)}


def clique_masks(adjacency: Sequence[int], size: int) -> Iterator[int]:
    """Enumerate all required cliques exactly once."""

    def visit(candidates: int, need: int, chosen: int) -> Iterator[int]:
        if need == 0:
            yield chosen
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            yield from visit(
                candidates & adjacency[vertex], need - 1, chosen | bit
            )

    yield from visit((1 << len(adjacency)) - 1, size, 0)


def bits(mask: int) -> list[int]:
    result = []
    while mask:
        bit = mask & -mask
        mask ^= bit
        result.append(bit.bit_length() - 1)
    return result


def coordinate_cycles(length: int) -> Iterator[tuple[int, ...]]:
    if length == 3:
        # Every permutation is an automorphism of the triangle core.
        yield from combinations(range(7), 3)
        return
    if length != 4:
        raise ValueError("unsupported cycle length")
    # A four-set has exactly three unoriented Hamilton cycles.  Fixing its
    # least coordinate first gives these representatives modulo D_8.
    for first, second, third, fourth in combinations(range(7), 4):
        yield (first, second, third, fourth)
        yield (first, second, fourth, third)
        yield (first, third, second, fourth)


@dataclass
class MatchStats:
    k7_seeds: int = 0
    coordinate_cycles: int = 0
    role_nodes: int = 0


def role_domains(
    adjacency: Sequence[int],
    seed_vertices: Sequence[int],
    cycle: Sequence[int],
    outside_mask: int,
) -> list[int]:
    seed_mask = sum(1 << vertex for vertex in seed_vertices)
    domains = []
    for position in range(len(cycle)):
        first = seed_vertices[cycle[position]]
        second = seed_vertices[cycle[(position + 1) % len(cycle)]]
        required = seed_mask & ~(1 << first) & ~(1 << second)
        domain = 0
        candidates = outside_mask
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            if adjacency[vertex] & required == required:
                domain |= bit
        domains.append(domain)
    return domains


def assign_roles(
    adjacency: Sequence[int], domains: Sequence[int], stats: MatchStats
) -> list[int] | None:
    length = len(domains)
    assignment = [-1] * length

    def visit(position: int, used: int) -> bool:
        stats.role_nodes += 1
        if position == length:
            return bool(
                adjacency[assignment[-1]] & (1 << assignment[0])
            )
        candidates = domains[position] & ~used
        if position:
            candidates &= adjacency[assignment[position - 1]]
        if position == length - 1 and assignment[0] >= 0:
            candidates &= adjacency[assignment[0]]
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            assignment[position] = bit.bit_length() - 1
            if visit(position + 1, used | bit):
                return True
        assignment[position] = -1
        return False

    return list(assignment) if visit(0, 0) else None


def verify_mapping(
    pattern: Sequence[int], target: Sequence[int], mapping: Sequence[int]
) -> None:
    if len(mapping) != len(pattern) or len(set(mapping)) != len(mapping):
        raise ValueError("mapping is not an injection of the pattern")
    if any(type(vertex) is not int or not 0 <= vertex < len(target)
           for vertex in mapping):
        raise ValueError("mapping has an invalid target vertex")
    for first, neighbours in enumerate(pattern):
        todo = neighbours
        while todo:
            bit = todo & -todo
            todo ^= bit
            second = bit.bit_length() - 1
            if not (target[mapping[first]] & (1 << mapping[second])):
                raise ValueError(f"mapping loses required edge {first}-{second}")


def find_embedding(
    target: Sequence[int], length: int, stats: MatchStats | None = None
) -> list[int] | None:
    target = validate_graph(target, TARGET_ORDER)
    pattern = PATTERNS[length]
    stats = stats if stats is not None else MatchStats()
    full = (1 << len(target)) - 1
    for seed_mask in clique_masks(target, 7):
        stats.k7_seeds += 1
        seed = bits(seed_mask)
        outside = full ^ seed_mask
        for cycle in coordinate_cycles(length):
            stats.coordinate_cycles += 1
            domains = role_domains(target, seed, cycle, outside)
            if any(not domain for domain in domains):
                continue
            roles = assign_roles(target, domains, stats)
            if roles is None:
                continue
            remaining = [coordinate for coordinate in range(7)
                         if coordinate not in cycle]
            seed_mapping = [seed[coordinate] for coordinate in cycle]
            seed_mapping.extend(seed[coordinate] for coordinate in remaining)
            mapping = seed_mapping + roles
            verify_mapping(pattern["adjacency"], target, mapping)
            return mapping
    return None


def load_targets() -> list[dict]:
    if sha256(TARGETS) != TARGETS_SHA256:
        raise ValueError("target TSV hash mismatch")
    if sha256(SELECTION) != SELECTION_SHA256:
        raise ValueError("selection hash mismatch")
    selection = json.loads(SELECTION.read_text(encoding="utf-8"))
    expected_indices = selection["selected_indices"]
    if (
        len(expected_indices) != TARGET_COUNT
        or stable_hash(expected_indices) != INDICES_SHA256
    ):
        raise ValueError("selection index list mismatch")
    rows = []
    with TARGETS.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = [int(value) for value in line.split()]
            if len(fields) != TARGET_ORDER + 1:
                raise ValueError(f"malformed target row {line_number}")
            rows.append({
                "index": fields[0],
                "adjacency": validate_graph(fields[1:], TARGET_ORDER),
            })
    if [row["index"] for row in rows] != expected_indices:
        raise ValueError("target TSV order differs from selection")
    return rows


def embedded_target(pattern: Sequence[int]) -> tuple[int, ...]:
    return validate_graph(
        tuple(pattern) + (0,) * (TARGET_ORDER - len(pattern)), TARGET_ORDER
    )


def complete_target() -> tuple[int, ...]:
    full = (1 << TARGET_ORDER) - 1
    return tuple(full ^ (1 << vertex) for vertex in range(TARGET_ORDER))


def controls() -> dict:
    cases = {}
    for length, pattern in PATTERNS.items():
        adjacency = pattern["adjacency"]
        embedded_mapping = find_embedding(embedded_target(adjacency), length)
        complete_mapping = find_embedding(complete_target(), length)
        empty_mapping = find_embedding((0,) * TARGET_ORDER, length)
        if embedded_mapping is None or complete_mapping is None:
            raise AssertionError("positive containment control failed")
        if empty_mapping is not None:
            raise AssertionError("empty containment control failed")
        verify_mapping(adjacency, embedded_target(adjacency), embedded_mapping)
        verify_mapping(adjacency, complete_target(), complete_mapping)
        cases[str(length)] = {
            "embedded_target": {"status": "HIT", "mapping": embedded_mapping},
            "complete_target": {"status": "HIT", "mapping": complete_mapping},
            "empty_target": {"status": "NO_HIT"},
        }
    return {"status": "PASS", "cases": cases}


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            capture_output=True, text=True,
        ).stdout.strip()

    try:
        porcelain = git("status", "--porcelain=v1", "--untracked-files=all")
        lines = porcelain.splitlines() if porcelain else []
        return {
            "available": True,
            "commit": git("rev-parse", "HEAD"),
            "branch": git("branch", "--show-current"),
            "dirty": bool(lines),
            "porcelain_lines": lines,
            "porcelain_sha256": stable_hash(lines),
        }
    except (OSError, subprocess.CalledProcessError) as error:
        return {"available": False, "error": str(error)}


def screen(lengths: Sequence[int], limit: int | None) -> dict:
    started = time.perf_counter()
    targets = load_targets()
    if limit is not None:
        targets = targets[:limit]
    hits = {length: [] for length in lengths}
    totals = {length: MatchStats() for length in lengths}
    for ordinal, target in enumerate(targets):
        for length in lengths:
            local = MatchStats()
            mapping = find_embedding(target["adjacency"], length, local)
            totals[length].k7_seeds += local.k7_seeds
            totals[length].coordinate_cycles += local.coordinate_cycles
            totals[length].role_nodes += local.role_nodes
            if mapping is not None:
                hits[length].append({
                    "ordinal": ordinal,
                    "index": target["index"],
                    "mapping": mapping,
                })
    combined_hits = {
        item["index"]
        for length in lengths
        for item in hits[length]
    }
    target_indices = [target["index"] for target in targets]
    residue = [index for index in target_indices if index not in combined_hits]
    return {
        "schema": 1,
        "kind": "d6_k7_mobius_cycle_containment_screen",
        "status": "COMPLETE",
        "required_edge_only": True,
        "targets": len(targets),
        "target_indices_sha256": stable_hash(
            target_indices
        ),
        "patterns": {
            str(length): {
                "pattern": PATTERNS[length],
                "hits": hits[length],
                "hit_count": len(hits[length]),
                "hit_indices_sha256": stable_hash(
                    [item["index"] for item in hits[length]]
                ),
                "stats": asdict(totals[length]),
            }
            for length in lengths
        },
        "combined_hit_indices": [
            index for index in target_indices if index in combined_hits
        ],
        "combined_hit_indices_sha256": stable_hash(
            [index for index in target_indices if index in combined_hits]
        ),
        "combined_residue_indices": residue,
        "combined_residue_indices_sha256": stable_hash(residue),
        "runtime": {
            "wall_seconds": time.perf_counter() - started,
            "finished_utc": utc_now(),
        },
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lengths", type=int, nargs="+", choices=(3, 4), default=(3, 4)
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = screen(tuple(dict.fromkeys(args.lengths)), args.limit)
    result["controls"] = controls()
    result["provenance"] = {
        "command": shlex.join(sys.argv),
        "source_hashes": {
            Path(__file__).name: sha256(Path(__file__)),
            PROOF.name: sha256(PROOF),
            VERIFIER.name: sha256(VERIFIER),
            TEST.name: sha256(TEST),
            TARGETS.relative_to(ROOT).as_posix(): sha256(TARGETS),
            SELECTION.name: sha256(SELECTION),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "pid": os.getpid(),
        },
        "git": git_provenance(),
    }
    if args.output is not None:
        atomic_json(args.output, result)
    print(json.dumps({
        "targets": result["targets"],
        "patterns": {
            key: {
                "hit_count": value["hit_count"],
                "stats": value["stats"],
            }
            for key, value in result["patterns"].items()
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()
