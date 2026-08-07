#!/usr/bin/env python3
"""Exact K7 full-support pinning obstruction from odd coordinate cycles.

For nonzero-factor vertices, put ``w=sqrt(7)u/t``.  At a required edge whose
two propagated masks intersect only in coordinate ``i``, the unit equation is
``w_i(x)w_i(y)=1``.  Along a nonbipartite component of this coordinate graph,
an odd cycle forces every component value to ``+1`` or ``-1``.

If all coordinates of one vertex's propagated mask are pinned this way, all
are actually nonzero and are signs.  Thus ``T=sum w_i`` is an integer and
``Q=sum w_i^2`` is the mask size, at most seven.  The normalized point
identity

    Q - T^2 + 2 sqrt(7) T - 8 = 0

separates over Q(sqrt(7)) to ``T=0`` and ``Q=8``, a contradiction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


def bits(mask: int):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate(graph: Sequence[int], masks: Sequence[int]):
    adjacency = tuple(map(int, graph))
    supports = tuple(map(int, masks))
    if len(adjacency) != len(supports):
        raise ValueError("graph and mask lengths disagree")
    full = (1 << len(adjacency)) - 1
    for vertex, row in enumerate(adjacency):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError("invalid graph")
        for neighbour in bits(row):
            if not adjacency[neighbour] & (1 << vertex):
                raise ValueError("asymmetric graph")
    if any(mask <= 0 or mask >= 128 for mask in supports):
        raise ValueError("masks must be nonempty and seven-bit")
    return adjacency, supports


def coordinate_graph(
    required: Sequence[int], masks: Sequence[int], coordinate: int
) -> tuple[int, ...]:
    bit = 1 << coordinate
    count = len(required)
    answer = [0] * count
    for first in range(count):
        for second in range(first):
            if (
                required[first] & (1 << second)
                and masks[first] & masks[second] == bit
            ):
                answer[first] |= 1 << second
                answer[second] |= 1 << first
    return tuple(answer)


def nonbipartite_components(graph: Sequence[int]):
    """Return deterministic ``(vertices, conflict_edge)`` certificates."""

    count = len(graph)
    unseen = set(range(count))
    output = []
    while unseen:
        root = min(unseen)
        colors = {root: 0}
        queue = [root]
        unseen.remove(root)
        conflict = None
        for vertex in queue:
            for neighbour in bits(graph[vertex]):
                if neighbour not in colors:
                    colors[neighbour] = colors[vertex] ^ 1
                    if neighbour in unseen:
                        unseen.remove(neighbour)
                    queue.append(neighbour)
                elif colors[neighbour] == colors[vertex] and conflict is None:
                    conflict = (min(vertex, neighbour), max(vertex, neighbour))
        if conflict is not None:
            output.append((tuple(sorted(colors)), conflict))
    return tuple(output)


@dataclass(frozen=True)
class FullPinCertificate:
    central: int
    central_mask: int
    coordinate_components: tuple[
        tuple[int, tuple[int, ...], tuple[int, int]], ...
    ]


def find_full_pin(
    graph: Sequence[int], masks: Sequence[int]
) -> FullPinCertificate | None:
    required, supports = validate(graph, masks)
    pin_data: dict[tuple[int, int], tuple[tuple[int, ...], tuple[int, int]]] = {}
    pinned = [0] * len(required)
    for coordinate in range(7):
        local = coordinate_graph(required, supports, coordinate)
        for component, conflict in nonbipartite_components(local):
            for vertex in component:
                pinned[vertex] |= 1 << coordinate
                pin_data[(vertex, coordinate)] = (component, conflict)
    for central, mask in enumerate(supports):
        if pinned[central] & mask == mask:
            entries = tuple(
                (coordinate, *pin_data[(central, coordinate)])
                for coordinate in bits(mask)
            )
            return FullPinCertificate(central, mask, entries)
    return None


def verify_certificate(
    graph: Sequence[int], masks: Sequence[int], certificate: FullPinCertificate
) -> None:
    required, supports = validate(graph, masks)
    central = int(certificate.central)
    if not 0 <= central < len(required) or supports[central] != certificate.central_mask:
        raise ValueError("bad central vertex or mask")
    expected_coordinates = tuple(bits(supports[central]))
    if tuple(entry[0] for entry in certificate.coordinate_components) != expected_coordinates:
        raise ValueError("certificate does not cover every central coordinate")
    for coordinate, component, conflict in certificate.coordinate_components:
        local = coordinate_graph(required, supports, coordinate)
        components = dict(nonbipartite_components(local))
        if central not in component or components.get(component) != conflict:
            raise ValueError("bad nonbipartite coordinate component certificate")


def certificate_json(certificate: FullPinCertificate) -> dict:
    return {
        "central": certificate.central,
        "central_mask": certificate.central_mask,
        "coordinate_components": [
            {
                "coordinate": coordinate,
                "component": list(component),
                "conflict_edge": list(conflict),
            }
            for coordinate, component, conflict in certificate.coordinate_components
        ],
    }


def certificate_from_json(value: dict) -> FullPinCertificate:
    return FullPinCertificate(
        central=int(value["central"]),
        central_mask=int(value["central_mask"]),
        coordinate_components=tuple(
            (
                int(item["coordinate"]),
                tuple(map(int, item["component"])),
                tuple(map(int, item["conflict_edge"])),
            )
            for item in value["coordinate_components"]
        ),
    )

