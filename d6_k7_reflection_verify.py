#!/usr/bin/env python3
"""Independent exact checker for K7--K6 reflection certificates.

This file deliberately does not import :mod:`d6_k7_reflection_overlap`.
It replays only a compact certificate against a supplied required-edge graph,
using Python ``Fraction`` arithmetic and the unit-simplex affine metric.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Sequence


Coordinate = tuple[Fraction, ...]


def validate_graph(adjacency: Sequence[int]) -> None:
    count = len(adjacency)
    full = (1 << count) - 1
    for vertex, row in enumerate(adjacency):
        if not isinstance(row, int) or isinstance(row, bool):
            raise ValueError("noninteger adjacency row")
        if row & ~full or row & (1 << vertex):
            raise ValueError("adjacency has an invalid bit or loop")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(
                adjacency[other] & (1 << vertex)
            ):
                raise ValueError("asymmetric adjacency")


def required_clique(adjacency: Sequence[int], vertices: Sequence[int]) -> bool:
    return (
        len(vertices) == len(set(vertices))
        and all(0 <= vertex < len(adjacency) for vertex in vertices)
        and all(
            adjacency[first] & (1 << second)
            for first, second in combinations(vertices, 2)
        )
    )


def basis(root: Sequence[int]) -> dict[int, Coordinate]:
    return {
        vertex: tuple(
            Fraction(int(position == coordinate)) for coordinate in range(7)
        )
        for position, vertex in enumerate(root)
    }


def reflect(
    coordinates: dict[int, Coordinate], old: int, facet: Sequence[int]
) -> Coordinate:
    if old not in coordinates or any(vertex not in coordinates for vertex in facet):
        raise ValueError("reflection dependency is unassigned")
    result = tuple(
        sum((coordinates[vertex][entry] for vertex in facet), Fraction()) / 3
        - coordinates[old][entry]
        for entry in range(7)
    )
    if sum(result, Fraction()) != 1:
        raise ValueError("reflection result is not affine")
    return result


def distance_squared(first: Coordinate, second: Coordinate) -> Fraction:
    if sum(first, Fraction()) != 1 or sum(second, Fraction()) != 1:
        raise ValueError("coordinate is not affine")
    return sum(
        ((left - right) ** 2 for left, right in zip(first, second)), Fraction()
    ) / 2


def replay_step(
    adjacency: Sequence[int], coordinates: dict[int, Coordinate], step: object
) -> tuple[int, Coordinate]:
    if not isinstance(step, dict) or set(step) != {
        "parent_clique",
        "child_clique",
        "facet",
        "old_apex",
        "new_apex",
    }:
        raise ValueError("malformed reflection step")
    parent = step["parent_clique"]
    child = step["child_clique"]
    facet = step["facet"]
    old = step["old_apex"]
    new = step["new_apex"]
    if (
        not isinstance(parent, list)
        or not isinstance(child, list)
        or not isinstance(facet, list)
        or len(parent) != 7
        or len(child) != 7
        or len(facet) != 6
    ):
        raise ValueError("reflection step has wrong set sizes")
    if not required_clique(adjacency, parent) or not required_clique(
        adjacency, child
    ):
        raise ValueError("reflection endpoint is not a required K7")
    if set(parent) & set(child) != set(facet):
        raise ValueError("reflection facet is not the clique intersection")
    if set(parent) - set(facet) != {old} or set(child) - set(facet) != {new}:
        raise ValueError("reflection apex labels are inconsistent")
    return new, reflect(coordinates, old, facet)


def verify(adjacency: Sequence[int], certificate: object) -> None:
    """Raise ``ValueError`` unless ``certificate`` is an exact contradiction."""

    validate_graph(adjacency)
    if not isinstance(certificate, dict) or set(certificate) != {
        "schema",
        "kind",
        "root_clique",
        "trace",
        "terminal",
    }:
        raise ValueError("malformed reflection certificate")
    if certificate["schema"] != 1 or certificate["kind"] != (
        "K7_K6_reflection_contradiction"
    ):
        raise ValueError("unsupported reflection certificate")
    root = certificate["root_clique"]
    if not isinstance(root, list) or len(root) != 7 or not required_clique(
        adjacency, root
    ):
        raise ValueError("invalid root K7")
    coordinates = basis(root)

    trace = certificate["trace"]
    if not isinstance(trace, list):
        raise ValueError("certificate trace is not a list")
    for step in trace:
        new, candidate = replay_step(adjacency, coordinates, step)
        if new in coordinates or candidate in coordinates.values():
            raise ValueError("nonterminal trace assignment is inconsistent")
        for other, coordinate in coordinates.items():
            if adjacency[new] & (1 << other):
                if distance_squared(candidate, coordinate) != 1:
                    raise ValueError("nonterminal trace violates a required edge")
        coordinates[new] = candidate

    terminal = certificate["terminal"]
    if not isinstance(terminal, dict):
        raise ValueError("terminal record is not an object")
    reason = terminal.get("reason")
    new, candidate = replay_step(adjacency, coordinates, terminal.get("step"))
    if reason == "coordinate_conflict":
        if terminal.get("vertex") != new:
            raise ValueError("wrong coordinate-conflict label")
        if new not in coordinates or coordinates[new] == candidate:
            raise ValueError("terminal record is not a coordinate conflict")
        return
    if reason == "collision":
        other = terminal.get("existing_vertex")
        if terminal.get("new_vertex") != new or new in coordinates:
            raise ValueError("collision labels are inconsistent")
        if other == new or other not in coordinates or coordinates[other] != candidate:
            raise ValueError("terminal record is not a collision")
        return
    if reason == "wrong_required_distance":
        edge = terminal.get("edge")
        if (
            not isinstance(edge, list)
            or len(edge) != 2
            or edge[0] == edge[1]
            or new not in edge
        ):
            raise ValueError("wrong-distance edge is malformed")
        other = edge[0] if edge[1] == new else edge[1]
        if other not in coordinates or not (adjacency[new] & (1 << other)):
            raise ValueError("wrong-distance pair is not a required edge")
        value = distance_squared(candidate, coordinates[other])
        if value == 1:
            raise ValueError("wrong-distance pair is actually unit")
        if terminal.get("squared_distance") != {
            "numerator": value.numerator,
            "denominator": value.denominator,
        }:
            raise ValueError("recorded distance is arithmetically wrong")
        return
    raise ValueError("unknown terminal contradiction")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path, help="JSON object with adjacency")
    parser.add_argument("certificate", type=Path, help="certificate JSON")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    graph = json.loads(args.graph.read_text(encoding="utf-8"))
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    verify(graph["adjacency"], certificate)
    print("K7 reflection-overlap certificate: PASS")


if __name__ == "__main__":
    main()
