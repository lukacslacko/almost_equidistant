#!/usr/bin/env python3
"""Exact propagation through the K7--K6 overlap complex.

Seven pairwise unit points in R6 form a full-dimensional regular simplex.
If two required K7 cliques share a K6 facet, their two nonshared vertices
must be the two points at unit distance from that facet.  Distinctness thus
forces one apex to be the reflection of the other.  Starting from any K7,
this module propagates those reflections with exact rational affine
coordinates and rejects only one of the following rigorous contradictions:

* one graph vertex receives two different coordinates;
* two different graph vertices receive the same coordinate; or
* a required unit edge receives squared distance different from one.

Candidate nonedges are never inspected as distance constraints.  A graph
which survives this necessary-condition filter is not claimed realizable.
The CLI profiles a deterministic sample of the hash-pinned current K7
residue; use ``--full`` only for an intentional full pass.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter, defaultdict, deque
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import d6_k7_reflection_verify as independent_verifier


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_SELECTION_SHA256 = (
    "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0"
)
EXPECTED_SELECTION_SIZE = 12_941

Coordinate = tuple[Fraction, ...]


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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


def bits(mask: int) -> Iterator[int]:
    """Yield set-bit positions in increasing order."""

    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate_graph(adj: Sequence[int]) -> None:
    """Check only the simple undirected graph representation."""

    n = len(adj)
    full = (1 << n) - 1
    for vertex, row in enumerate(adj):
        if not isinstance(row, int) or isinstance(row, bool):
            raise ValueError(f"adjacency row {vertex} is not an integer")
        if row & ~full:
            raise ValueError(f"adjacency row {vertex} has out-of-range bits")
        if row & (1 << vertex):
            raise ValueError(f"loop at vertex {vertex}")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(adj[other] & (1 << vertex)):
                raise ValueError(f"asymmetric pair {other},{vertex}")


def clique_masks(adj: Sequence[int], size: int = 7) -> Iterator[int]:
    """Enumerate required cliques once, in deterministic mask order."""

    def visit(candidates: int, need: int, chosen: int) -> Iterator[int]:
        if need == 0:
            yield chosen
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            yield from visit(candidates & adj[vertex], need - 1, chosen | bit)

    yield from visit((1 << len(adj)) - 1, size, 0)


def is_clique(adj: Sequence[int], vertices: Sequence[int]) -> bool:
    if len(vertices) != len(set(vertices)):
        return False
    if any(vertex < 0 or vertex >= len(adj) for vertex in vertices):
        return False
    return all(adj[first] & (1 << second) for first, second in combinations(vertices, 2))


def root_coordinates(root: Sequence[int]) -> dict[int, Coordinate]:
    """Give a root K7 its affine-coordinate basis e_0,...,e_6."""

    if len(root) != 7 or len(set(root)) != 7:
        raise ValueError("a root must contain seven distinct vertices")
    return {
        vertex: tuple(Fraction(int(position == coordinate)) for coordinate in range(7))
        for position, vertex in enumerate(root)
    }


def reflect_coordinate(
    coordinates: dict[int, Coordinate], old: int, facet: Sequence[int]
) -> Coordinate:
    """Reflect ``old`` in a regular-simplex K6 facet, exactly.

    The altitude foot is the facet centroid, so the reflected apex is
    ``(1/3) * sum(facet) - old``.
    """

    if len(facet) != 6 or len(set(facet)) != 6 or old in facet:
        raise ValueError("a reflection step requires one apex and a K6 facet")
    if old not in coordinates or any(vertex not in coordinates for vertex in facet):
        raise ValueError("reflection inputs do not yet have coordinates")
    answer = tuple(
        sum((coordinates[vertex][entry] for vertex in facet), Fraction()) / 3
        - coordinates[old][entry]
        for entry in range(7)
    )
    if sum(answer, Fraction()) != 1:
        raise AssertionError("reflection did not preserve affine-coordinate sum")
    return answer


def squared_distance(first: Coordinate, second: Coordinate) -> Fraction:
    """Squared Euclidean distance in a unit-simplex affine basis."""

    if len(first) != 7 or len(second) != 7:
        raise ValueError("K7 affine coordinates have length seven")
    if sum(first, Fraction()) != 1 or sum(second, Fraction()) != 1:
        raise ValueError("affine coordinates must sum to one")
    return sum(
        ((left - right) ** 2 for left, right in zip(first, second)), Fraction()
    ) / 2


def fraction_json(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def step_record(parent: int, child: int, facet: int) -> dict:
    old_mask = parent & ~facet
    new_mask = child & ~facet
    if old_mask.bit_count() != 1 or new_mask.bit_count() != 1:
        raise ValueError("overlap edge does not exchange exactly one apex")
    return {
        "parent_clique": list(bits(parent)),
        "child_clique": list(bits(child)),
        "facet": list(bits(facet)),
        "old_apex": next(bits(old_mask)),
        "new_apex": next(bits(new_mask)),
    }


def replay_step(
    adj: Sequence[int], coordinates: dict[int, Coordinate], step: dict
) -> tuple[int, Coordinate]:
    """Validate one recorded overlap edge and compute its reflected apex."""

    required = {
        "parent_clique", "child_clique", "facet", "old_apex", "new_apex"
    }
    if set(step) != required:
        raise ValueError("reflection step has unexpected fields")
    parent = tuple(step["parent_clique"])
    child = tuple(step["child_clique"])
    facet = tuple(step["facet"])
    old = step["old_apex"]
    new = step["new_apex"]
    if not is_clique(adj, parent) or not is_clique(adj, child):
        raise ValueError("recorded overlap endpoint is not a required K7")
    if len(parent) != 7 or len(child) != 7:
        raise ValueError("recorded overlap endpoints must have size seven")
    if set(parent) & set(child) != set(facet) or len(facet) != 6:
        raise ValueError("recorded facet is not the six-vertex intersection")
    if set(parent) - set(facet) != {old}:
        raise ValueError("recorded old apex is wrong")
    if set(child) - set(facet) != {new}:
        raise ValueError("recorded new apex is wrong")
    return new, reflect_coordinate(coordinates, old, facet)


@dataclass(frozen=True)
class ReflectionResult:
    rejected: bool
    reason: str | None
    k7_cliques: int
    overlap_edges: int
    overlap_components: int
    nontrivial_components: int
    maximum_component_cliques: int
    maximum_component_vertices: int
    certificate: dict | None


def _certificate(
    root: tuple[int, ...], trace: list[dict], terminal: dict
) -> dict:
    return {
        "schema": 1,
        "kind": "K7_K6_reflection_contradiction",
        "root_clique": list(root),
        "trace": trace,
        "terminal": terminal,
    }


def analyze_graph(adj: Sequence[int]) -> ReflectionResult:
    """Apply exact reflection propagation to every overlap component."""

    validate_graph(adj)
    cliques = tuple(clique_masks(adj, 7))
    by_facet: dict[int, list[int]] = defaultdict(list)
    for clique in cliques:
        for vertex in bits(clique):
            by_facet[clique ^ (1 << vertex)].append(clique)

    neighbors: dict[int, list[tuple[int, int]]] = {
        clique: [] for clique in cliques
    }
    overlap_edges = 0
    for facet in sorted(by_facet):
        extensions = sorted(by_facet[facet])
        for first, second in combinations(extensions, 2):
            neighbors[first].append((second, facet))
            neighbors[second].append((first, facet))
            overlap_edges += 1
    for clique in neighbors:
        neighbors[clique].sort()

    seen: set[int] = set()
    components = 0
    nontrivial = 0
    maximum_cliques = 0
    maximum_vertices = 0

    def rejected(reason: str, certificate: dict) -> ReflectionResult:
        return ReflectionResult(
            rejected=True,
            reason=reason,
            k7_cliques=len(cliques),
            overlap_edges=overlap_edges,
            overlap_components=components,
            nontrivial_components=nontrivial,
            maximum_component_cliques=maximum_cliques,
            maximum_component_vertices=maximum_vertices,
            certificate=certificate,
        )

    for root_mask in cliques:
        if root_mask in seen:
            continue
        components += 1
        root = tuple(bits(root_mask))
        coordinates = root_coordinates(root)
        trace: list[dict] = []
        queue = deque([root_mask])
        seen.add(root_mask)
        component_cliques = 0
        component_nontrivial = bool(neighbors[root_mask])
        if component_nontrivial:
            nontrivial += 1

        while queue:
            parent = queue.popleft()
            component_cliques += 1
            for child, facet in neighbors[parent]:
                step = step_record(parent, child, facet)
                new, candidate = replay_step(adj, coordinates, step)
                if new in coordinates:
                    if coordinates[new] != candidate:
                        terminal = {
                            "reason": "coordinate_conflict",
                            "step": step,
                            "vertex": new,
                        }
                        maximum_cliques = max(maximum_cliques, component_cliques)
                        maximum_vertices = max(maximum_vertices, len(coordinates))
                        return rejected(
                            "coordinate_conflict",
                            _certificate(root, trace, terminal),
                        )
                else:
                    collision = next(
                        (
                            vertex
                            for vertex, coordinate in sorted(coordinates.items())
                            if coordinate == candidate
                        ),
                        None,
                    )
                    if collision is not None:
                        terminal = {
                            "reason": "collision",
                            "step": step,
                            "new_vertex": new,
                            "existing_vertex": collision,
                        }
                        maximum_cliques = max(maximum_cliques, component_cliques)
                        maximum_vertices = max(maximum_vertices, len(coordinates))
                        return rejected(
                            "collision", _certificate(root, trace, terminal)
                        )

                    wrong_edge: tuple[int, Fraction] | None = None
                    for other, coordinate in sorted(coordinates.items()):
                        if adj[new] & (1 << other):
                            distance = squared_distance(candidate, coordinate)
                            if distance != 1:
                                wrong_edge = (other, distance)
                                break
                    if wrong_edge is not None:
                        other, distance = wrong_edge
                        terminal = {
                            "reason": "wrong_required_distance",
                            "step": step,
                            "edge": [min(new, other), max(new, other)],
                            "squared_distance": fraction_json(distance),
                        }
                        maximum_cliques = max(maximum_cliques, component_cliques)
                        maximum_vertices = max(maximum_vertices, len(coordinates))
                        return rejected(
                            "wrong_required_distance",
                            _certificate(root, trace, terminal),
                        )
                    coordinates[new] = candidate
                    trace.append(step)

                if child not in seen:
                    seen.add(child)
                    queue.append(child)

        maximum_cliques = max(maximum_cliques, component_cliques)
        maximum_vertices = max(maximum_vertices, len(coordinates))

    return ReflectionResult(
        rejected=False,
        reason=None,
        k7_cliques=len(cliques),
        overlap_edges=overlap_edges,
        overlap_components=components,
        nontrivial_components=nontrivial,
        maximum_component_cliques=maximum_cliques,
        maximum_component_vertices=maximum_vertices,
        certificate=None,
    )


def verify_certificate(adj: Sequence[int], certificate: dict) -> None:
    """Independently replay a compact exact contradiction certificate."""

    validate_graph(adj)
    if certificate.get("schema") != 1:
        raise ValueError("unsupported reflection certificate schema")
    if certificate.get("kind") != "K7_K6_reflection_contradiction":
        raise ValueError("wrong reflection certificate kind")
    if set(certificate) != {"schema", "kind", "root_clique", "trace", "terminal"}:
        raise ValueError("reflection certificate has unexpected fields")
    root = tuple(certificate["root_clique"])
    if len(root) != 7 or not is_clique(adj, root):
        raise ValueError("certificate root is not a required K7")
    coordinates = root_coordinates(root)

    for step in certificate["trace"]:
        new, candidate = replay_step(adj, coordinates, step)
        if new in coordinates:
            raise ValueError("trace tries to reassign an existing vertex")
        if candidate in coordinates.values():
            raise ValueError("trace silently introduces a collision")
        for other, coordinate in coordinates.items():
            if adj[new] & (1 << other) and squared_distance(candidate, coordinate) != 1:
                raise ValueError("trace silently violates a required edge")
        coordinates[new] = candidate

    terminal = certificate["terminal"]
    reason = terminal.get("reason")
    new, candidate = replay_step(adj, coordinates, terminal.get("step", {}))
    if reason == "coordinate_conflict":
        if terminal.get("vertex") != new:
            raise ValueError("coordinate-conflict vertex does not match its step")
        if new not in coordinates or coordinates[new] == candidate:
            raise ValueError("recorded coordinate conflict is not a conflict")
        return
    if reason == "collision":
        other = terminal.get("existing_vertex")
        if terminal.get("new_vertex") != new or new in coordinates:
            raise ValueError("recorded collision has inconsistent labels")
        if other == new or other not in coordinates or coordinates[other] != candidate:
            raise ValueError("recorded collision is not a collision")
        return
    if reason == "wrong_required_distance":
        edge = terminal.get("edge")
        if not isinstance(edge, list) or len(edge) != 2 or new not in edge:
            raise ValueError("recorded wrong-distance edge is malformed")
        other = edge[0] if edge[1] == new else edge[1]
        if other not in coordinates or not (adj[new] & (1 << other)):
            raise ValueError("recorded pair is not a derived required edge")
        distance = squared_distance(candidate, coordinates[other])
        if distance == 1:
            raise ValueError("recorded required edge actually has unit distance")
        expected = terminal.get("squared_distance")
        if expected != fraction_json(distance):
            raise ValueError("recorded squared distance does not match exact replay")
        return
    raise ValueError(f"unknown terminal reflection reason {reason!r}")


def select_graphs(
    input_path: Path,
    selection_path: Path,
    sample_size: int | None,
) -> tuple[list[dict], dict]:
    """Load and deterministically sample the hash-pinned 12,941 residue."""

    input_hash = sha256(input_path)
    selection_hash = sha256(selection_path)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"rank-survivor input hash mismatch: {input_hash}")
    if selection_hash != EXPECTED_SELECTION_SHA256:
        raise ValueError(f"current-residue selection hash mismatch: {selection_hash}")

    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    indices = selection.get("selected_indices")
    if (
        not isinstance(indices, list)
        or len(indices) != EXPECTED_SELECTION_SIZE
        or len(set(indices)) != len(indices)
    ):
        raise ValueError("current-residue selection population is malformed")
    if selection.get("selected") != EXPECTED_SELECTION_SIZE:
        raise ValueError("selection count disagrees with selected_indices")

    payload = json.loads(input_path.read_text(encoding="utf-8"))
    rows = payload.get("graphs")
    if not isinstance(rows, list):
        raise ValueError("rank-survivor input has no graph list")
    by_index = {row.get("index"): row for row in rows}
    if len(by_index) != len(rows):
        raise ValueError("rank-survivor input repeats an index")
    try:
        selected = [by_index[index] for index in indices]
    except KeyError as error:
        raise ValueError(f"selection index absent from rank survivors: {error}") from error

    if sample_size is None:
        chosen = selected
        positions = list(range(len(selected)))
        method = "full selected population in selection order"
    else:
        if sample_size <= 0 or sample_size > len(selected):
            raise ValueError("sample size must lie between one and the residue size")
        positions = [position * len(selected) // sample_size for position in range(sample_size)]
        if len(set(positions)) != len(positions):
            raise AssertionError("evenly spaced sample positions are not distinct")
        chosen = [selected[position] for position in positions]
        method = (
            "deterministic evenly spaced positions floor(i*N/m), "
            "0 <= i < m, in selection order"
        )
    metadata = {
        "input": display_path(input_path),
        "input_sha256": input_hash,
        "selection": display_path(selection_path),
        "selection_sha256": selection_hash,
        "selection_population": len(selected),
        "sample_size": len(chosen),
        "sample_method": method,
        "sample_positions_sha256": stable_hash(positions),
        "sample_indices_sha256": stable_hash([row["index"] for row in chosen]),
    }
    return chosen, metadata


def profile(graphs: Iterable[dict], metadata: dict) -> dict:
    started = time.monotonic()
    results: list[tuple[int, ReflectionResult]] = []
    for graph in graphs:
        adjacency = graph.get("adjacency")
        if not isinstance(adjacency, list):
            raise ValueError(f"graph {graph.get('index')} has no adjacency list")
        result = analyze_graph(adjacency)
        if result.k7_cliques == 0:
            raise ValueError(f"selected graph {graph.get('index')} has no K7")
        if result.certificate is not None:
            verify_certificate(adjacency, result.certificate)
            independent_verifier.verify(adjacency, result.certificate)
        results.append((int(graph["index"]), result))

    reasons = Counter(result.reason for _, result in results if result.rejected)
    overlap_rows = [
        {
            "index": index,
            "k7_cliques": result.k7_cliques,
            "overlap_edges": result.overlap_edges,
            "nontrivial_components": result.nontrivial_components,
            "maximum_component_cliques": result.maximum_component_cliques,
            "maximum_component_vertices": result.maximum_component_vertices,
            "rejected": result.rejected,
            "reason": result.reason,
        }
        for index, result in results
        if result.overlap_edges
    ]
    rejected = [
        {
            "index": index,
            "reason": result.reason,
            "certificate": result.certificate,
        }
        for index, result in results
        if result.rejected
    ]
    return {
        "schema": 1,
        "kind": "d6_k7_reflection_overlap_profile",
        "description": (
            "Exact rational propagation through required K7 cliques sharing "
            "K6 facets. Nonedges remain unconstrained."
        ),
        "implementation": {
            "generator": "d6_k7_reflection_overlap.py",
            "generator_sha256": sha256(ROOT / "d6_k7_reflection_overlap.py"),
            "independent_checker": "d6_k7_reflection_verify.py",
            "independent_checker_sha256": sha256(
                ROOT / "d6_k7_reflection_verify.py"
            ),
            "controls": "d6_k7_reflection_overlap_test.py",
            "controls_sha256": sha256(
                ROOT / "d6_k7_reflection_overlap_test.py"
            ),
        },
        "scope": metadata,
        "counts": {
            "graphs": len(results),
            "rejected": len(rejected),
            "survived": len(results) - len(rejected),
            "graphs_with_multiple_K7": sum(
                result.k7_cliques > 1 for _, result in results
            ),
            "graphs_with_K6_overlap": len(overlap_rows),
            "total_K7_cliques": sum(result.k7_cliques for _, result in results),
            "total_overlap_edges": sum(result.overlap_edges for _, result in results),
            "maximum_K7_cliques": max(
                (result.k7_cliques for _, result in results), default=0
            ),
            "maximum_overlap_edges": max(
                (result.overlap_edges for _, result in results), default=0
            ),
            "maximum_component_cliques": max(
                (result.maximum_component_cliques for _, result in results), default=0
            ),
            "maximum_component_vertices": max(
                (result.maximum_component_vertices for _, result in results), default=0
            ),
        },
        "rejection_reasons": dict(sorted(reasons.items())),
        "overlap_graphs": overlap_rows,
        "rejections": rejected,
        "checks": {
            "all_selected_graphs_contain_K7": True,
            "every_rejection_certificate_replayed_by_generator": True,
            "every_rejection_certificate_replayed_by_independent_checker": True,
            "candidate_nonedges_used_as_constraints": False,
        },
        "runtime": {
            "elapsed_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
            "workers": 1,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / ".runs/d6_k7_rank_survivors.json",
    )
    parser.add_argument(
        "--selection",
        type=Path,
        default=ROOT / "d6_k7_positive_dual_selection.json",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=1_024,
        help="deterministic evenly spaced pilot size (default: 1024)",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="profile all 12,941 selected graphs instead of a sample",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_reflection_overlap_sample.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    graphs, metadata = select_graphs(
        args.input, args.selection, None if args.full else args.sample_size
    )
    report = profile(graphs, metadata)
    atomic_json(args.output, report)
    print(
        f"wrote {args.output}: {report['counts']['graphs']} graphs, "
        f"{report['counts']['graphs_with_K6_overlap']} with K6 overlaps, "
        f"{report['counts']['rejected']} exact rejections in "
        f"{report['runtime']['elapsed_seconds']:.3f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
