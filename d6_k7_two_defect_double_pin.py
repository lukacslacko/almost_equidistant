#!/usr/bin/env python3
"""Exact K7 obstruction from two singleton-intersection pinning triangles.

Work on the nonzero-factor vertices after a K7 zero-factor cover has been
fixed.  Normalize a defect vector by ``w=sqrt(7)u/t``.  Every required unit
edge ``xy`` then satisfies ``w(x).w(y)=1``.

If three required-unit vertices form a triangle and every pair of their
allowed support masks intersects in the same singleton coordinate ``i``, all
three dot products reduce to products of their ``i``-coordinates.  Hence

    ab=ac=bc=1,

so all three coordinates are the same sign in ``{+1,-1}``.

Now suppose a vertex ``x`` has a propagated allowed mask ``{i,j}`` and lies
in one such pinning triangle at ``i`` and another at ``j``.  Both coordinates
of ``x`` are nonzero, so its actual support is exactly ``{i,j}``, and each is
``+1`` or ``-1``.  But the exact two-defect point equation is

    (w_i-sqrt(7))(w_j-sqrt(7)) = 3,

which none of the four sign pairs satisfies.  This is a required-edge-only
contradiction.  Candidate nonedges are never assigned a distance.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence


def bits(mask: int):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate(graph: Sequence[int], masks: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    adjacency = tuple(map(int, graph))
    supports = tuple(map(int, masks))
    if len(adjacency) != len(supports):
        raise ValueError("graph and mask lengths disagree")
    full = (1 << len(adjacency)) - 1
    for vertex, row in enumerate(adjacency):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError("invalid graph row")
        for neighbour in bits(row):
            if not adjacency[neighbour] & (1 << vertex):
                raise ValueError("asymmetric graph")
    if any(mask <= 0 or mask >= 128 for mask in supports):
        raise ValueError("masks must be nonempty and seven-bit")
    return adjacency, supports


@dataclass(frozen=True)
class DoublePinCertificate:
    central: int
    central_mask: int
    first_coordinate: int
    first_triangle: tuple[int, int, int]
    second_coordinate: int
    second_triangle: tuple[int, int, int]


def singleton_pinning_triangles(
    graph: Sequence[int], masks: Sequence[int]
) -> dict[tuple[int, int], tuple[tuple[int, int, int], ...]]:
    """Return all pinning triangles keyed by ``(vertex, coordinate)``."""

    adjacency, supports = validate(graph, masks)
    found: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
    for triangle in combinations(range(len(adjacency)), 3):
        first, second, third = triangle
        if not (
            adjacency[first] & (1 << second)
            and adjacency[first] & (1 << third)
            and adjacency[second] & (1 << third)
        ):
            continue
        intersections = (
            supports[first] & supports[second],
            supports[first] & supports[third],
            supports[second] & supports[third],
        )
        if not (
            intersections[0]
            and intersections[0] == intersections[1] == intersections[2]
            and intersections[0].bit_count() == 1
        ):
            continue
        coordinate = intersections[0].bit_length() - 1
        for vertex in triangle:
            found.setdefault((vertex, coordinate), []).append(triangle)
    return {key: tuple(value) for key, value in sorted(found.items())}


def find_double_pin(
    graph: Sequence[int], masks: Sequence[int]
) -> DoublePinCertificate | None:
    """Return the first deterministic two-defect double-pin obstruction."""

    adjacency, supports = validate(graph, masks)
    pins = singleton_pinning_triangles(adjacency, supports)
    for central, mask in enumerate(supports):
        if mask.bit_count() != 2:
            continue
        first, second = tuple(bits(mask))
        left = pins.get((central, first), ())
        right = pins.get((central, second), ())
        if left and right:
            return DoublePinCertificate(
                central=central,
                central_mask=mask,
                first_coordinate=first,
                first_triangle=left[0],
                second_coordinate=second,
                second_triangle=right[0],
            )
    return None


def verify_certificate(
    graph: Sequence[int], masks: Sequence[int], certificate: DoublePinCertificate
) -> None:
    """Check a serialized obstruction without trusting the locator."""

    adjacency, supports = validate(graph, masks)
    central = int(certificate.central)
    if not 0 <= central < len(adjacency):
        raise ValueError("central vertex out of range")
    if supports[central] != certificate.central_mask or supports[central].bit_count() != 2:
        raise ValueError("central mask is not the recorded two-set")
    coordinates = tuple(bits(supports[central]))
    if (certificate.first_coordinate, certificate.second_coordinate) != coordinates:
        raise ValueError("recorded coordinates do not match central mask order")
    for coordinate, triangle in (
        (certificate.first_coordinate, certificate.first_triangle),
        (certificate.second_coordinate, certificate.second_triangle),
    ):
        if tuple(sorted(triangle)) != triangle or len(set(triangle)) != 3:
            raise ValueError("pinning triangle labels are not canonical")
        if central not in triangle:
            raise ValueError("pinning triangle omits the central vertex")
        for first, second in combinations(triangle, 2):
            if not adjacency[first] & (1 << second):
                raise ValueError("pinning triangle has a missing required edge")
            if supports[first] & supports[second] != 1 << coordinate:
                raise ValueError("triangle intersection is not the recorded singleton")


def certificate_json(certificate: DoublePinCertificate) -> dict:
    return {
        "central": certificate.central,
        "central_mask": certificate.central_mask,
        "first_coordinate": certificate.first_coordinate,
        "first_triangle": list(certificate.first_triangle),
        "second_coordinate": certificate.second_coordinate,
        "second_triangle": list(certificate.second_triangle),
    }


def certificate_from_json(value: dict) -> DoublePinCertificate:
    return DoublePinCertificate(
        central=int(value["central"]),
        central_mask=int(value["central_mask"]),
        first_coordinate=int(value["first_coordinate"]),
        first_triangle=tuple(map(int, value["first_triangle"])),
        second_coordinate=int(value["second_coordinate"]),
        second_triangle=tuple(map(int, value["second_triangle"])),
    )

