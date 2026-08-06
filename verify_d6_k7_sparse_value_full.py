#!/usr/bin/env python3
"""Independent, restartable verifier for the full K7 sparse-value archive.

Neither ``d6_k7_small_support_value`` nor its profiler/production runner is
imported.  The checker reads those files only to validate the hashes recorded
by the production report.  For every one of the 3,195 archived rejections it
reconstructs the recorded failing K7 seed, all baseline-passing covers, every
labeled zero-factor support family, singleton propagation, all actual
one/two-coordinate support assignments, and the exact finite value rules.

The frozen rank reference is used only for the already independently checked
baseline-cover classification.  The independently committed support verifier
supplies the separate labeled-support and graph-construction kernel.
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
import shutil
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

import d6_k7_rank_reference as rank_reference
import verify_d6_k7_support_full as support_verifier


ROOT = Path(__file__).resolve().parent
EXPECTED = {
    "input": "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9",
    "support_archive": "3c35b228c4768b06f881dc9046329196b683d390f4763f83816c6a3464ecf938",
    "support_uncompressed": "117e057a79135fd48cac316ee4f1e08e11a2be84e00d9d10393bc17d4a518a31",
    "support_report": "5138d207543967b393baae9b14b47b49a276273cbc4468ee843856891917691e",
    "support_checkpoint": "daeee836f1f2c632656bbae812ef884cb5f2eca882fb2a5d5131f1c18edbdb5b",
    "support_verifier": "c6b9dce47618bd9a4247eb5a9f093682265447729e1868c4dc0677a24f51c628",
    "rank_reference": "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0",
    "value_archive": "9f70e0e267fe97bc2f6a2890cae47d8b1c882974d5f054f16b4673bfb45ee003",
    "value_uncompressed": "0ce9d2d2aeb18403fce17a0611594843946d4c07764fbb24946069f0997aab4c",
    "value_report": "cea64f3dde804e0766c5b77e373aa50338d5bee6e5e54f44749f4e387cf529aa",
    "value_checkpoint": "ae91cc60463f728697d0cc8fb807ca430c0842ba2e56b205b37d9ba877ddad9e",
    "value_source": "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f",
    "value_proof": "e8ec6a012f9ea824e9cddb76510840537dba0e603e207cd3220c1b06866cb42a",
    "value_profiler": "c8231ff47263b348e3daba2b32cd3eb693ff70f3db5b10aa808b1b5296589573",
    "value_runner": "af7dd2f0173c8593044bf6a32d53fc0ab7a3cc4ee45c9b9af21809b3a729be07",
}
EXPECTED_RANK_SURVIVORS = 17_764
EXPECTED_SUPPORT_SURVIVORS = 16_228
EXPECTED_VALUE_REJECTED = 3_195
QOS_CLASS_USER_INITIATED = 0x19

BASE_COUNT_FIELDS = (
    "seeds",
    "covers",
    "baseline_passing_covers",
    "labeled_z_families",
    "propagation_passing_families",
    "small_support_intersection_failures",
    "small_support_intersection_passes",
    "sparse_value_failures",
    "sparse_value_passes",
    "value_passing_covers",
)
PROPAGATION_REASONS = (
    "empty_propagated_mask",
    "disjoint_required_edge",
    "clique",
    "degree",
    "basis",
    "mask",
    "subspace_K",
    "component_B",
)
VALUE_REASONS = (
    "disjoint_required_small_support",
    "three_one_defects",
    "duplicate_one_defect",
    "three_same_two_defects",
    "two_defects_repeat_at_one_defect",
    "forbidden_two_defect_cycle",
    "parallel_two_defect_type_touches_another",
    "one_defects_joined_by_two_defects",
    "one_defect_component_branches",
    "incompatible_two_defect_branch_distance",
    "inconsistent_two_defect_branch_signs",
)
VALUE_DECISION_FIELDS = (
    "index",
    "decision",
    "first_failing_seed",
    *BASE_COUNT_FIELDS,
    *(f"propagation_failure_{reason}" for reason in PROPAGATION_REASONS),
    *(f"sparse_value_branch_{reason}" for reason in VALUE_REASONS),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def decompressed_bytes_and_hash(path: Path) -> tuple[bytes, str]:
    with gzip.open(path, "rb") as stream:
        content = stream.read()
    return content, hashlib.sha256(content).hexdigest()


def stable_hash(value: object) -> str:
    content = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(content).hexdigest()


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def deterministic_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with source.open("rb") as input_stream, temporary.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0
        ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream, length=1 << 20)
        raw.flush()
        os.fsync(raw.fileno())
    os.replace(temporary, destination)


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
            return
        raise OSError(error, "pthread_set_qos_class_self_np failed")


# Exact Q(sqrt(7)) arithmetic used as a mandatory proof control.
Quadratic = tuple[int, int]
Matrix = tuple[tuple[Quadratic, Quadratic], tuple[Quadratic, Quadratic]]


def qadd(first: Quadratic, second: Quadratic) -> Quadratic:
    return first[0] + second[0], first[1] + second[1]


def qscale(scale: int, value: Quadratic) -> Quadratic:
    return scale * value[0], scale * value[1]


def qmul(first: Quadratic, second: Quadratic) -> Quadratic:
    return (
        first[0] * second[0] + 7 * first[1] * second[1],
        first[0] * second[1] + first[1] * second[0],
    )


def matrix_multiply(first: Matrix, second: Matrix) -> Matrix:
    return tuple(  # type: ignore[return-value]
        tuple(
            qadd(qmul(first[row][0], second[0][column]),
                 qmul(first[row][1], second[1][column]))
            for column in range(2)
        )
        for row in range(2)
    )


def validate_mobius_lemmas() -> dict:
    """Independently confirm M^6 and all five real fixed-point obstructions."""

    zero = (0, 0)
    one = (1, 0)
    identity: Matrix = ((one, zero), (zero, one))
    matrix: Matrix = ((one, (0, -1)), ((0, 1), (-4, 0)))
    power = identity
    discriminants = []
    powers = []
    for exponent in range(1, 7):
        power = matrix_multiply(power, matrix)
        powers.append(power)
        a, b = power[0]
        c, d = power[1]
        d_minus_a = qadd(d, qscale(-1, a))
        discriminant = qadd(
            qmul(d_minus_a, d_minus_a), qscale(4, qmul(b, c))
        )
        discriminants.append(discriminant)
    expected_discriminants = [(-3, 0), (-27, 0), (-108, 0), (-243, 0), (-243, 0)]
    if discriminants[:5] != expected_discriminants:
        raise AssertionError("independent Mobius fixed-point discriminants disagree")
    expected_sixth: Matrix = (((-27, 0), zero), (zero, (-27, 0)))
    if powers[5] != expected_sixth:
        raise AssertionError("independent Mobius order-six identity disagrees")

    # For x=c*sqrt(7), T sends c to (c-1)/(7c-4).
    orbit = [Fraction(1, 4)]
    for _ in range(5):
        coefficient = orbit[-1]
        orbit.append((coefficient - 1) / (7 * coefficient - 4))
    expected_orbit = [
        Fraction(1, 4), Fraction(1, 3), Fraction(2, 5),
        Fraction(1, 2), Fraction(1), Fraction(0),
    ]
    if orbit != expected_orbit:
        raise AssertionError("independent one-defect orbit disagrees")
    return {
        "M6": [[list(value) for value in row] for row in powers[5]],
        "fixed_point_discriminants": [list(value) for value in discriminants[:5]],
        "one_defect_orbit_sqrt7_coefficients": [str(value) for value in orbit],
    }


def bit_positions(mask: int) -> Iterator[int]:
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def independent_cycle_lengths(edges: frozenset[tuple[int, int]]) -> frozenset[int]:
    """Enumerate simple-cycle lengths by a rooted DFS, not cyclic permutations."""

    adjacency = [set() for _ in range(7)]
    for first, second in edges:
        if not 0 <= first < second < 7:
            raise ValueError("bad two-support edge")
        adjacency[first].add(second)
        adjacency[second].add(first)
    lengths: set[int] = set()
    for root in range(7):
        path = [root]
        used = {root}

        def visit(vertex: int) -> None:
            for neighbour in adjacency[vertex]:
                if neighbour == root:
                    if len(path) >= 3:
                        lengths.add(len(path))
                    continue
                # Root must be the least cycle vertex; this removes rotations.
                if neighbour < root or neighbour in used:
                    continue
                used.add(neighbour)
                path.append(neighbour)
                visit(neighbour)
                path.pop()
                used.remove(neighbour)

        visit(root)
    return frozenset(lengths)


def edge_components(edges: frozenset[tuple[int, int]]) -> tuple[frozenset[int], ...]:
    adjacency = [set() for _ in range(7)]
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = {vertex for vertex in range(7) if adjacency[vertex]}
    output = []
    while unseen:
        root = min(unseen)
        component = {root}
        stack = [root]
        unseen.remove(root)
        while stack:
            vertex = stack.pop()
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        output.append(frozenset(component))
    return tuple(output)


class ParityUnionFind:
    """Union-find with xor labels for independent branch-sign consistency."""

    def __init__(self, size: int = 7) -> None:
        self.parent = list(range(size))
        self.parity = [0] * size

    def find(self, vertex: int) -> tuple[int, int]:
        if self.parent[vertex] == vertex:
            return vertex, 0
        root, above = self.find(self.parent[vertex])
        self.parity[vertex] ^= above
        self.parent[vertex] = root
        return root, self.parity[vertex]

    def constrain(self, first: int, second: int, xor: int) -> bool:
        root_first, parity_first = self.find(first)
        root_second, parity_second = self.find(second)
        if root_first == root_second:
            return (parity_first ^ parity_second) == xor
        self.parent[root_second] = root_first
        self.parity[root_second] = parity_first ^ parity_second ^ xor
        return True


@lru_cache(maxsize=200_000)
def independent_value_failure(sorted_supports: tuple[int, ...]) -> str | None:
    """Apply the exact one/two-support rules to a partial assignment."""

    supports = sorted_supports
    if any(mask <= 0 or mask >= 128 or mask.bit_count() not in (1, 2)
           for mask in supports):
        raise ValueError("value rules require nonempty one/two-bit supports")
    singletons = [next(bit_positions(mask)) for mask in supports
                  if mask.bit_count() == 1]
    if len(singletons) > 2:
        return "three_one_defects"
    if len(singletons) != len(set(singletons)):
        return "duplicate_one_defect"

    multiplicities: Counter[tuple[int, int]] = Counter(
        tuple(bit_positions(mask)) for mask in supports if mask.bit_count() == 2
    )
    if any(number > 2 for number in multiplicities.values()):
        return "three_same_two_defects"
    singleton_set = set(singletons)
    for coordinate in singleton_set:
        incident = sum(
            number for edge, number in multiplicities.items()
            if coordinate in edge
        )
        if incident > 1:
            return "two_defects_repeat_at_one_defect"

    edges = frozenset(multiplicities)
    if any(length % 6 for length in independent_cycle_lengths(edges)):
        return "forbidden_two_defect_cycle"
    degrees = [0] * 7
    adjacency = [set() for _ in range(7)]
    for first, second in edges:
        degrees[first] += 1
        degrees[second] += 1
        adjacency[first].add(second)
        adjacency[second].add(first)
    if any(
        number == 2 and (degrees[first] >= 2 or degrees[second] >= 2)
        for (first, second), number in multiplicities.items()
    ):
        return "parallel_two_defect_type_touches_another"

    for component in edge_components(edges):
        component_singletons = component & singleton_set
        if len(component_singletons) > 1:
            return "one_defects_joined_by_two_defects"
        if component_singletons and any(degrees[v] >= 3 for v in component):
            return "one_defect_component_branches"

    branches = {vertex for vertex in range(7) if degrees[vertex] >= 3}
    constraints: set[tuple[int, int, int]] = set()
    for start in branches:
        for neighbour in adjacency[start]:
            previous, current, length = start, neighbour, 1
            while current not in branches and degrees[current] == 2:
                following = adjacency[current] - {previous}
                if len(following) != 1:
                    raise AssertionError("bad degree-two branch path")
                previous, current = current, next(iter(following))
                length += 1
            if current not in branches or current == start:
                continue
            if length % 3:
                return "incompatible_two_defect_branch_distance"
            constraints.add((min(start, current), max(start, current),
                             (length // 3) & 1))
    signs = ParityUnionFind()
    for first, second, parity in sorted(constraints):
        if not signs.constrain(first, second, parity):
            return "inconsistent_two_defect_branch_signs"
    return None


@dataclass(frozen=True)
class SmallSupportCheck:
    feasible: bool
    search_nodes: int
    complete_assignments: int
    prunes: tuple[tuple[str, int], ...]
    witness: tuple[tuple[int, int], ...] | None


def independent_small_support_check(
    graph_n: Sequence[int], propagated_masks: Sequence[int]
) -> SmallSupportCheck:
    """Exhaust actual supports for all vertices forced to size one or two."""

    selected = tuple(
        vertex for vertex, mask in enumerate(propagated_masks)
        if mask.bit_count() <= 2
    )
    if not selected:
        return SmallSupportCheck(True, 1, 1, (), ())
    selected_set = set(selected)
    domains = {}
    for vertex in selected:
        mask = propagated_masks[vertex]
        domains[vertex] = tuple(
            support for support in range(1, 128)
            if not (support & ~mask)
        )
    order = tuple(sorted(
        selected,
        key=lambda vertex: (
            len(domains[vertex]),
            -sum(1 for neighbour in bit_positions(graph_n[vertex])
                 if neighbour in selected_set),
            vertex,
        ),
    ))
    assignment: dict[int, int] = {}
    prunes: Counter[str] = Counter()
    nodes = 0
    complete = 0

    def visit(depth: int) -> tuple[tuple[int, int], ...] | None:
        nonlocal nodes, complete
        nodes += 1
        if depth == len(order):
            complete += 1
            return tuple(sorted(assignment.items()))
        vertex = order[depth]
        for support in domains[vertex]:
            if any(
                graph_n[vertex] & (1 << prior)
                and not (support & prior_support)
                for prior, prior_support in assignment.items()
            ):
                prunes["disjoint_required_small_support"] += 1
                continue
            assignment[vertex] = support
            failure = independent_value_failure(
                tuple(sorted(assignment.values()))
            )
            if failure is None:
                witness = visit(depth + 1)
                if witness is not None:
                    return witness
            else:
                prunes[failure] += 1
            del assignment[vertex]
        return None

    witness = visit(0)
    return SmallSupportCheck(
        feasible=witness is not None,
        search_nodes=nodes,
        complete_assignments=complete,
        prunes=tuple(sorted(prunes.items())),
        witness=witness,
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


def propagated_masks(
    n_allowed: Sequence[int], fixed_supports: Sequence[int]
) -> tuple[tuple[int, ...], int]:
    output = []
    deletions = 0
    for allowed in n_allowed:
        closed, count = support_verifier.independent_singleton_closure(
            allowed, fixed_supports
        )
        output.append(closed)
        deletions += count
    return tuple(output), deletions


def simple_propagation_failure(
    graph_n: Sequence[int], masks: Sequence[int]
) -> str | None:
    if any(mask == 0 for mask in masks):
        return "empty_propagated_mask"
    for first, second in combinations(range(len(masks)), 2):
        if graph_n[first] & (1 << second) and not (masks[first] & masks[second]):
            return "disjoint_required_edge"
    return None


def verify_rejected_graph(task: tuple[dict, int]) -> dict:
    """Replay the recorded sparse-value failing seed of one graph."""

    assert _support_solver is not None
    assert _zero_forcing is not None
    assert _clique_solver is not None
    graph, seed_mask = task
    index = int(graph["index"])
    adj = tuple(int(row) for row in graph["adjacency"])
    rank_reference.validate_graph(adj)
    outside, defects, ladj, eligible = support_verifier.independent_seed_instance(
        adj, seed_mask
    )
    covers = tuple(support_verifier.independent_eligible_covers(ladj, eligible))
    total_term_rank = rank_reference.matching_size(defects)
    counts: Counter[str] = Counter()
    branch_counts: Counter[str] = Counter()
    value_search_nodes = 0
    complete_assignments = 0
    coordinate_deletions = 0
    baseline_passing = 0

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
        zvertices = tuple(support_verifier.bit_positions(zmask))
        nvertices = tuple(
            vertex for vertex in range(len(outside))
            if not (zmask & (1 << vertex))
        )
        z_allowed = tuple(defects[vertex] for vertex in zvertices)
        n_allowed = tuple(defects[vertex] for vertex in nvertices)
        graph_n = support_verifier.independent_induced_graph(
            adj, tuple(outside[vertex] for vertex in nvertices)
        )
        for fixed_supports in support_verifier.independent_labeled_support_families(
            z_allowed
        ):
            counts["labeled_z_families"] += 1
            masks, deletions = propagated_masks(n_allowed, fixed_supports)
            coordinate_deletions += deletions
            propagation_failure = simple_propagation_failure(graph_n, masks)
            if propagation_failure is not None:
                counts[f"propagation_failure:{propagation_failure}"] += 1
                continue
            counts["propagation_passing_families"] += 1
            check = independent_small_support_check(graph_n, masks)
            value_search_nodes += check.search_nodes
            complete_assignments += check.complete_assignments
            branch_counts.update(dict(check.prunes))
            if check.feasible:
                return {
                    "verified": False,
                    "index": index,
                    "recorded_failing_seed": seed_mask,
                    "surviving_cover": zmask,
                    "surviving_z_supports": list(fixed_supports),
                    "surviving_small_support_assignment": list(check.witness or ()),
                }
            counts["sparse_value_failed_families"] += 1
    if baseline_passing == 0:
        raise AssertionError(
            f"value-only failing seed at graph {index} has no baseline pass"
        )
    if counts["propagation_passing_families"] != counts["sparse_value_failed_families"]:
        raise AssertionError(f"value family population mismatch at graph {index}")
    return {
        "verified": True,
        "index": index,
        "recorded_failing_seed": seed_mask,
        "eligible_covers": len(covers),
        "baseline_passing_covers": baseline_passing,
        "labeled_z_families": counts["labeled_z_families"],
        "propagation_failure_counts": {
            reason: counts[f"propagation_failure:{reason}"]
            for reason in ("empty_propagated_mask", "disjoint_required_edge")
        },
        "propagation_passing_families": counts["propagation_passing_families"],
        "sparse_value_failed_families": counts["sparse_value_failed_families"],
        "value_search_nodes": value_search_nodes,
        "complete_small_support_assignments": complete_assignments,
        "value_branch_prunes": dict(sorted(branch_counts.items())),
        "coordinate_deletions": coordinate_deletions,
    }


def parse_tsv_archive(path: Path, fields: Sequence[str]) -> tuple[list[dict], str, int]:
    content, uncompressed_hash = decompressed_bytes_and_hash(path)
    reader = csv.DictReader(io.StringIO(content.decode("ascii")), delimiter="\t")
    if tuple(reader.fieldnames or ()) != tuple(fields):
        raise ValueError(f"archive schema mismatch: {path}")
    rows = list(reader)
    return rows, uncompressed_hash, len(content)


def validate_bindings(args: argparse.Namespace) -> tuple[list[dict], list[dict], dict]:
    """Validate both archived campaigns and reconstruct the selected population."""

    paths_and_hashes = {
        "input": (args.input, EXPECTED["input"]),
        "support_archive": (args.support_archive, EXPECTED["support_archive"]),
        "support_report": (args.support_report, EXPECTED["support_report"]),
        "support_checkpoint": (args.support_checkpoint, EXPECTED["support_checkpoint"]),
        "support_verifier": (ROOT / "verify_d6_k7_support_full.py", EXPECTED["support_verifier"]),
        "rank_reference": (ROOT / "d6_k7_rank_reference.py", EXPECTED["rank_reference"]),
        "value_archive": (args.value_archive, EXPECTED["value_archive"]),
        "value_report": (args.value_report, EXPECTED["value_report"]),
        "value_checkpoint": (args.value_checkpoint, EXPECTED["value_checkpoint"]),
        "value_source": (ROOT / "d6_k7_small_support_value.py", EXPECTED["value_source"]),
        "value_proof": (ROOT / "d6_k7_small_support_value.md", EXPECTED["value_proof"]),
        "value_profiler": (ROOT / "profile_d6_k7_small_support_value.py", EXPECTED["value_profiler"]),
        "value_runner": (ROOT / "run_d6_k7_small_support_value_full.py", EXPECTED["value_runner"]),
    }
    observed_hashes = {}
    for name, (path, expected) in paths_and_hashes.items():
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"{name} hash {observed} != expected {expected}")
        observed_hashes[name] = observed

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != EXPECTED_RANK_SURVIVORS:
        raise ValueError("unexpected rank-residue population")
    graph_indices = [int(graph["index"]) for graph in graphs]
    if len(set(graph_indices)) != len(graph_indices):
        raise ValueError("duplicate graph index in rank residue")

    support_rows, support_uncompressed, support_bytes = parse_tsv_archive(
        args.support_archive, support_verifier.DECISION_FIELDS
    )
    if support_uncompressed != EXPECTED["support_uncompressed"]:
        raise ValueError("support archive uncompressed hash mismatch")
    if len(support_rows) != len(graphs):
        raise ValueError("support archive population mismatch")
    selected_graphs = []
    support_rejected = []
    for graph, row in zip(graphs, support_rows):
        if int(row["index"]) != int(graph["index"]):
            raise ValueError("support archive/input order mismatch")
        if row["refined_decision"] == "SURVIVOR":
            selected_graphs.append(graph)
        elif row["refined_decision"] == "REJECTED":
            support_rejected.append(int(row["index"]))
        else:
            raise ValueError("bad support archive decision")
    if len(selected_graphs) != EXPECTED_SUPPORT_SURVIVORS:
        raise ValueError("unexpected support-survivor population")

    value_rows, value_uncompressed, value_bytes = parse_tsv_archive(
        args.value_archive, VALUE_DECISION_FIELDS
    )
    if value_uncompressed != EXPECTED["value_uncompressed"]:
        raise ValueError("value archive uncompressed hash mismatch")
    selected_indices = [int(graph["index"]) for graph in selected_graphs]
    value_indices = [int(row["index"]) for row in value_rows]
    if value_indices != selected_indices:
        raise ValueError("value archive/support-survivor population mismatch")
    if any(row["decision"] not in {"REJECTED", "SURVIVOR"}
           for row in value_rows):
        raise ValueError("bad value archive decision")
    if any(int(row[field]) < 0 for row in value_rows
           for field in VALUE_DECISION_FIELDS if field != "decision"):
        raise ValueError("negative integer in value archive")
    rejected_rows = [row for row in value_rows if row["decision"] == "REJECTED"]
    rejected_indices = [int(row["index"]) for row in rejected_rows]
    if len(rejected_rows) != EXPECTED_VALUE_REJECTED:
        raise ValueError("unexpected value rejection count")
    if any(
        (row["decision"] == "REJECTED") != (int(row["first_failing_seed"]) != 0)
        for row in value_rows
    ):
        raise ValueError("value decision/failing-seed mismatch")

    value_report = json.loads(args.value_report.read_text(encoding="utf-8"))
    decisions = value_report.get("decisions", {})
    if (
        decisions.get("graphs") != len(value_rows)
        or decisions.get("rejected") != len(rejected_rows)
        or decisions.get("survivors") != len(value_rows) - len(rejected_rows)
        or decisions.get("rejected_indices") != rejected_indices
    ):
        raise ValueError("value report decision summary mismatch")
    totals = {
        field: sum(int(row[field]) for row in value_rows)
        for field in VALUE_DECISION_FIELDS[3:]
    }
    if decisions.get("totals") != totals:
        raise ValueError("value report column totals mismatch")
    if decisions.get("survivor_indices_sha256") != stable_hash(
        [int(row["index"]) for row in value_rows if row["decision"] == "SURVIVOR"]
    ):
        raise ValueError("value survivor population hash mismatch")
    artifacts = value_report.get("artifacts", {})
    if (
        artifacts.get("archive_sha256") != EXPECTED["value_archive"]
        or artifacts.get("decisions_sha256") != EXPECTED["value_uncompressed"]
        or artifacts.get("checkpoint_sha256") != EXPECTED["value_checkpoint"]
    ):
        raise ValueError("value report artifact binding mismatch")
    configuration = value_report.get("configuration", {})
    if stable_hash(configuration) != value_report.get("configuration_sha256"):
        raise ValueError("value report configuration hash mismatch")
    if (
        configuration.get("run_graphs") != EXPECTED_SUPPORT_SURVIVORS
        or configuration.get("run_indices_sha256") != stable_hash(selected_indices)
        or configuration.get("runner_source_sha256") != EXPECTED["value_runner"]
        or configuration.get("decision_fields") != list(VALUE_DECISION_FIELDS)
    ):
        raise ValueError("value report configuration mismatch")
    recorded_dependencies = configuration.get("dependencies_sha256", {})
    for name, expected in {
        "d6_k7_rank_reference.py": EXPECTED["rank_reference"],
        "d6_k7_small_support_value.py": EXPECTED["value_source"],
        "d6_k7_small_support_value.md": EXPECTED["value_proof"],
        "profile_d6_k7_small_support_value.py": EXPECTED["value_profiler"],
    }.items():
        if recorded_dependencies.get(name) != expected:
            raise ValueError(f"value dependency report mismatch for {name}")

    checkpoint = json.loads(args.value_checkpoint.read_text(encoding="utf-8"))
    if (
        checkpoint.get("status") != "COMPLETE"
        or checkpoint.get("committed_rows") != EXPECTED_SUPPORT_SURVIVORS
        or checkpoint.get("decisions_sha256") != EXPECTED["value_uncompressed"]
        or checkpoint.get("archive_sha256") != EXPECTED["value_archive"]
        or checkpoint.get("compatibility") != configuration
    ):
        raise ValueError("value checkpoint/report binding mismatch")

    bindings = {
        "hashes": observed_hashes,
        "support_archive_uncompressed_sha256": support_uncompressed,
        "support_archive_uncompressed_bytes": support_bytes,
        "value_archive_uncompressed_sha256": value_uncompressed,
        "value_archive_uncompressed_bytes": value_bytes,
        "rank_residue_graphs": len(graphs),
        "support_rejected": len(support_rejected),
        "support_survivors": len(selected_graphs),
        "support_survivor_indices_sha256": stable_hash(selected_indices),
        "value_rejected": len(rejected_rows),
        "value_rejected_indices_sha256": stable_hash(rejected_indices),
        "production_totals_recomputed": totals,
        "production_value_modules_imported": False,
    }
    return selected_graphs, rejected_rows, bindings


def load_completed(
    partial: Path, expected_indices: Sequence[int], committed: int
) -> list[dict]:
    if not partial.is_file():
        if committed:
            raise ValueError("checkpoint names results but partial JSONL is absent")
        return []
    lines = partial.read_text(encoding="utf-8").splitlines()
    if len(lines) < committed:
        raise ValueError("partial JSONL is shorter than checkpoint")
    results = [json.loads(line) for line in lines[:committed]]
    if [int(item["index"]) for item in results] != list(expected_indices[:committed]):
        raise ValueError("partial JSONL prefix population mismatch")
    if len(lines) != committed:
        atomic_text(
            partial,
            "".join(json.dumps(item, sort_keys=True) + "\n" for item in results),
        )
    return results


def summarize(results: Sequence[dict]) -> dict:
    propagation = Counter()
    value_prunes = Counter()
    for result in results:
        propagation.update(result["propagation_failure_counts"])
        value_prunes.update(result["value_branch_prunes"])
    return {
        "verified_rejections": len(results),
        "eligible_covers": sum(item["eligible_covers"] for item in results),
        "baseline_passing_covers": sum(
            item["baseline_passing_covers"] for item in results
        ),
        "labeled_z_families": sum(item["labeled_z_families"] for item in results),
        "propagation_failure_counts": dict(sorted(propagation.items())),
        "propagation_passing_families": sum(
            item["propagation_passing_families"] for item in results
        ),
        "sparse_value_failed_families": sum(
            item["sparse_value_failed_families"] for item in results
        ),
        "value_search_nodes": sum(item["value_search_nodes"] for item in results),
        "complete_small_support_assignments": sum(
            item["complete_small_support_assignments"] for item in results
        ),
        "value_branch_prunes": dict(sorted(value_prunes.items())),
        "coordinate_deletions": sum(item["coordinate_deletions"] for item in results),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=ROOT / ".runs/d6_k7_rank_survivors.json")
    parser.add_argument("--support-archive", type=Path,
                        default=ROOT / "d6_k7_support_rank_survivors_decisions.tsv.gz")
    parser.add_argument("--support-report", type=Path,
                        default=ROOT / "d6_k7_support_rank_survivors_report.json")
    parser.add_argument("--support-checkpoint", type=Path,
                        default=ROOT / "d6_k7_support_rank_survivors_checkpoint.json")
    parser.add_argument("--value-archive", type=Path,
                        default=ROOT / "d6_k7_small_support_value_full_decisions.tsv.gz")
    parser.add_argument("--value-report", type=Path,
                        default=ROOT / "d6_k7_small_support_value_full_report.json")
    parser.add_argument("--value-checkpoint", type=Path,
                        default=ROOT / "d6_k7_small_support_value_full_checkpoint.json")
    parser.add_argument("--decisions", type=Path,
                        default=ROOT / ".runs/d6_k7_sparse_value_full_verification.jsonl")
    parser.add_argument("--archive", type=Path,
                        default=ROOT / "d6_k7_sparse_value_full_verification.jsonl.gz")
    parser.add_argument("--checkpoint", type=Path,
                        default=ROOT / ".runs/d6_k7_sparse_value_full_verification_checkpoint.json")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "d6_k7_sparse_value_full_verification_report.json")
    parser.add_argument("--pid-file", type=Path,
                        default=ROOT / ".runs/d6_k7_sparse_value_full_verification.pid.json")
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--checkpoint-every", type=int, default=32)
    parser.add_argument("--progress-every", type=int, default=64)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.checkpoint_every < 1 or args.progress_every < 1:
        raise SystemExit("worker/checkpoint/progress counts must be positive")
    if args.limit is not None and args.limit < 1:
        raise SystemExit("limit must be positive")

    algebra_controls = validate_mobius_lemmas()
    forbidden_imports = {
        "d6_k7_small_support_value",
        "profile_d6_k7_small_support_value",
        "run_d6_k7_small_support_value_full",
    }
    imported = forbidden_imports & set(sys.modules)
    if imported:
        raise SystemExit(
            "production value module unexpectedly imported: "
            + ", ".join(sorted(imported))
        )
    selected_graphs, rejected_rows, bindings = validate_bindings(args)
    graph_by_index = {int(graph["index"]): graph for graph in selected_graphs}
    if args.limit is not None:
        rejected_rows = rejected_rows[: args.limit]
    expected_indices = [int(row["index"]) for row in rejected_rows]
    tasks = [
        (graph_by_index[index], int(row["first_failing_seed"]))
        for index, row in zip(expected_indices, rejected_rows)
    ]
    source_hash = sha256(Path(__file__))
    compatibility = {
        "schema": 1,
        "verifier_source_sha256": source_hash,
        "input_sha256": EXPECTED["input"],
        "value_archive_sha256": EXPECTED["value_archive"],
        "value_archive_uncompressed_sha256": EXPECTED["value_uncompressed"],
        "selected_rejections": len(tasks),
        "selected_indices_sha256": stable_hash(expected_indices),
        "limit": args.limit,
    }
    partial = args.decisions.with_name(args.decisions.name + ".partial")
    checkpoint_previous = None
    if args.resume:
        checkpoint_previous = json.loads(args.checkpoint.read_text(encoding="utf-8"))
        if checkpoint_previous.get("compatibility") != compatibility:
            raise SystemExit("resume checkpoint is incompatible")
        committed = int(checkpoint_previous.get("committed_results", -1))
        # A crash in the narrow finalization window may occur after the
        # completed partial was renamed but before archive/report creation.
        if not partial.exists() and args.decisions.exists():
            os.replace(args.decisions, partial)
        results = load_completed(partial, expected_indices, committed)
    else:
        stale = [path for path in (partial, args.decisions, args.checkpoint,
                                   args.archive, args.report) if path.exists()]
        if stale:
            raise SystemExit(
                "fresh-run output already exists; use new paths or --resume: "
                + ", ".join(map(str, stale))
            )
        atomic_text(partial, "")
        results = []
    started_unix = time.time()
    started = time.monotonic()
    command = [sys.executable, *sys.argv]

    def checkpoint(
        status: str,
        error: str | None = None,
        artifacts: dict | None = None,
    ) -> None:
        value = {
            "schema": 1,
            "status": status,
            "compatibility": compatibility,
            "committed_results": len(results),
            "partial_decisions": str(partial),
            "final_decisions": str(args.decisions),
            "archive": str(args.archive),
            "report": str(args.report),
            "workers": args.workers,
            "command": command,
            "started_unix": started_unix,
            "updated_unix": time.time(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        }
        if error is not None:
            value["error"] = error
        if artifacts is not None:
            value["artifacts"] = artifacts
        atomic_json(args.checkpoint, value)

    atomic_json(args.pid_file, {
        "pid": os.getpid(),
        "command": command,
        "workers": args.workers,
        "started_unix": started_unix,
        "compatibility": compatibility,
        "checkpoint": str(args.checkpoint),
    })
    checkpoint("RUNNING")
    remaining = tasks[len(results):]
    try:
        with partial.open("a", encoding="utf-8", newline="") as stream:
            if args.workers == 1:
                initialize_worker()
                iterator = map(verify_rejected_graph, remaining)
                executor = None
            else:
                executor = ProcessPoolExecutor(
                    max_workers=args.workers, initializer=initialize_worker
                )
                iterator = executor.map(verify_rejected_graph, remaining, chunksize=1)
            try:
                for result in iterator:
                    if not result.get("verified"):
                        raise AssertionError(
                            "archived rejection has an independent survivor: "
                            + json.dumps(result, sort_keys=True)
                        )
                    stream.write(json.dumps(result, sort_keys=True) + "\n")
                    results.append(result)
                    done = len(results)
                    if done % args.checkpoint_every == 0 or done == len(tasks):
                        stream.flush()
                        os.fsync(stream.fileno())
                        checkpoint("RUNNING")
                    if done % args.progress_every == 0 or done == len(tasks):
                        elapsed = time.monotonic() - started
                        print(
                            f"verified {done}/{len(tasks)}; "
                            f"{(done-len(tasks)+len(remaining))/max(elapsed,1e-9):.2f} "
                            f"graph/s; elapsed {elapsed:.1f}s",
                            file=sys.stderr,
                            flush=True,
                        )
            finally:
                if executor is not None:
                    executor.shutdown()
    except KeyboardInterrupt:
        checkpoint("ABORT", "KeyboardInterrupt")
        raise
    except BaseException as error:
        checkpoint("INFRA_ERROR", "".join(traceback.format_exception(error)))
        raise

    if [int(item["index"]) for item in results] != expected_indices:
        checkpoint("INFRA_ERROR", "final population/order mismatch")
        raise SystemExit("final population/order mismatch")
    os.replace(partial, args.decisions)
    deterministic_gzip(args.decisions, args.archive)
    summary = summarize(results)
    elapsed = time.monotonic() - started
    complete = len(results) == EXPECTED_VALUE_REJECTED and args.limit is None
    report = {
        "schema": 1,
        "method": "independent_K7_sparse_one_two_support_rejection_replay",
        "status": "PASS" if complete else "BOUNDED_CONTROL_PASS",
        "complete_rejection_population_verified": complete,
        "command": command,
        "started_unix": started_unix,
        "wall_seconds": elapsed,
        "workers": args.workers,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "verifier_source_sha256": source_hash,
        "production_value_modules_imported": False,
        "bindings": bindings,
        "algebra_controls": algebra_controls,
        "recorded_rejections": EXPECTED_VALUE_REJECTED,
        "selected_rejections": len(tasks),
        "decision_jsonl": str(args.decisions),
        "decision_jsonl_sha256": sha256(args.decisions),
        "decision_archive": str(args.archive),
        "decision_archive_sha256": sha256(args.archive),
        **summary,
        "trust_scope": (
            "Exact integer, rational, and bit-mask replay. The frozen rank "
            "reference is trusted only for baseline-cover classification; "
            "the production sparse-value locator and evaluator are not imported."
        ),
    }
    atomic_json(args.report, report)
    checkpoint("COMPLETE", artifacts={
        "decision_jsonl_sha256": report["decision_jsonl_sha256"],
        "decision_archive_sha256": report["decision_archive_sha256"],
        "report_sha256": sha256(args.report),
    })
    print(
        f"{report['status']}: {len(results)} rejections; "
        f"{summary['labeled_z_families']} labeled Z families; "
        f"{summary['propagation_passing_families']} value families; "
        f"wall {elapsed:.1f}s; report {sha256(args.report)}",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()
