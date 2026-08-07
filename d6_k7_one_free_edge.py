#!/usr/bin/env python3
"""Exact K7 obstruction from a singleton edge with one free coordinate per end.

On a nonzero-factor K7 branch put ``w=sqrt(7)u/t``.  A nonbipartite
singleton-intersection coordinate component pins each component value to a
sign.  If a vertex has a propagated mask of size ``k`` and all but one of its
components are pinned, let ``P`` be the sum of the ``k-1`` pinned signs.  The
diagonal simplex equation forces the remaining component to

    F(k,P) = (9-k+P^2-2*sqrt(7)*P) / (2*(sqrt(7)-P)).

For ``1 <= k,l <= 7``, exact arithmetic in ``Q(sqrt(7))`` shows that
``F(k,P)F(l,Q)=1`` is possible only when one type is ``(4,+/-3)`` and the
other is ``(7,0)``.  But masks of sizes four and seven in a seven-coordinate
universe intersect in at least four coordinates.  Consequently, every
required edge whose propagated masks meet in exactly the two endpoints' sole
unpinned coordinate is impossible.

Candidate nonedges are never constrained, and propagated masks remain
supersets of actual supports.  The equation for ``F`` forces the last allowed
component to be nonzero, so optional zeros cause no gap.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


Quadratic = tuple[Fraction, Fraction]
PinComponent = tuple[int, tuple[int, ...], tuple[int, int]]


def bits(mask: int):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate(
    graph: Sequence[int], masks: Sequence[int]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
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


def sign_sums(mask_size: int) -> tuple[int, ...]:
    """All sums of ``mask_size-1`` independent signs."""

    if not 1 <= mask_size <= 7:
        raise ValueError("mask size must lie in [1,7]")
    return tuple(range(-(mask_size - 1), mask_size, 2))


def one_free_value(mask_size: int, sign_sum: int) -> Quadratic:
    """Return ``F(k,P)`` as rational and sqrt(7) coefficients.

    Rationalizing the displayed formula in the module docstring gives

      F(k,P) = [P(P^2-k-5) + (9-k-P^2)sqrt(7)] / [2(7-P^2)].
    """

    if sign_sum not in sign_sums(mask_size):
        raise ValueError("sign sum has the wrong range or parity")
    denominator = 2 * (7 - sign_sum * sign_sum)
    return (
        Fraction(sign_sum * (sign_sum * sign_sum - mask_size - 5), denominator),
        Fraction(9 - mask_size - sign_sum * sign_sum, denominator),
    )


def quadratic_product(first: Quadratic, second: Quadratic) -> Quadratic:
    """Multiply two ``a+b*sqrt(7)`` values exactly."""

    a, b = first
    c, d = second
    return a * c + 7 * b * d, a * d + b * c


def reciprocal_possible(first_size: int, second_size: int) -> bool:
    """Whether some two admissible one-free values have product one."""

    one = (Fraction(1), Fraction(0))
    return any(
        quadratic_product(
            one_free_value(first_size, first_sum),
            one_free_value(second_size, second_sum),
        )
        == one
        for first_sum in sign_sums(first_size)
        for second_sum in sign_sums(second_size)
    )


def coordinate_graph(
    required: Sequence[int], masks: Sequence[int], coordinate: int
) -> tuple[int, ...]:
    bit = 1 << coordinate
    answer = [0] * len(required)
    for first in range(len(required)):
        for second in range(first):
            if (
                required[first] & (1 << second)
                and masks[first] & masks[second] == bit
            ):
                answer[first] |= 1 << second
                answer[second] |= 1 << first
    return tuple(answer)


def nonbipartite_components(
    graph: Sequence[int],
) -> tuple[tuple[tuple[int, ...], tuple[int, int]], ...]:
    """Return deterministic ``(vertices, conflict_edge)`` certificates."""

    unseen = set(range(len(graph)))
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
                    unseen.discard(neighbour)
                    queue.append(neighbour)
                elif colors[neighbour] == colors[vertex] and conflict is None:
                    conflict = (min(vertex, neighbour), max(vertex, neighbour))
        if conflict is not None:
            output.append((tuple(sorted(colors)), conflict))
    return tuple(output)


def pinning_data(
    required: Sequence[int], masks: Sequence[int]
) -> tuple[tuple[int, ...], dict[tuple[int, int], PinComponent]]:
    pinned = [0] * len(required)
    data: dict[tuple[int, int], PinComponent] = {}
    for coordinate in range(7):
        local = coordinate_graph(required, masks, coordinate)
        for component, conflict in nonbipartite_components(local):
            for vertex in component:
                pinned[vertex] |= 1 << coordinate
                data[(vertex, coordinate)] = (
                    coordinate,
                    component,
                    conflict,
                )
    return tuple(pinned), data


@dataclass(frozen=True)
class OneFreeEdgeCertificate:
    first: int
    second: int
    shared_coordinate: int
    first_mask: int
    second_mask: int
    first_pin_components: tuple[PinComponent, ...]
    second_pin_components: tuple[PinComponent, ...]


def find_one_free_edge(
    graph: Sequence[int], masks: Sequence[int]
) -> OneFreeEdgeCertificate | None:
    required, supports = validate(graph, masks)
    pinned, data = pinning_data(required, supports)
    for second in range(len(required)):
        for first in bits(required[second] & ((1 << second) - 1)):
            intersection = supports[first] & supports[second]
            if intersection.bit_count() != 1:
                continue
            shared = intersection.bit_length() - 1
            shared_bit = 1 << shared
            first_needed = supports[first] ^ shared_bit
            second_needed = supports[second] ^ shared_bit
            if (
                pinned[first] & supports[first] != first_needed
                or pinned[second] & supports[second] != second_needed
                or reciprocal_possible(
                    supports[first].bit_count(), supports[second].bit_count()
                )
            ):
                continue
            return OneFreeEdgeCertificate(
                first=first,
                second=second,
                shared_coordinate=shared,
                first_mask=supports[first],
                second_mask=supports[second],
                first_pin_components=tuple(
                    data[(first, coordinate)] for coordinate in bits(first_needed)
                ),
                second_pin_components=tuple(
                    data[(second, coordinate)] for coordinate in bits(second_needed)
                ),
            )
    return None


def verify_certificate(
    graph: Sequence[int],
    masks: Sequence[int],
    certificate: OneFreeEdgeCertificate,
) -> None:
    required, supports = validate(graph, masks)
    first = int(certificate.first)
    second = int(certificate.second)
    if not 0 <= first < second < len(required):
        raise ValueError("bad edge endpoints")
    if not required[first] & (1 << second):
        raise ValueError("recorded pair is not a required edge")
    if (
        supports[first] != certificate.first_mask
        or supports[second] != certificate.second_mask
    ):
        raise ValueError("recorded masks differ")
    intersection = supports[first] & supports[second]
    if (
        intersection.bit_count() != 1
        or intersection.bit_length() - 1 != certificate.shared_coordinate
    ):
        raise ValueError("edge intersection is not the recorded singleton")
    shared_bit = intersection
    if reciprocal_possible(
        supports[first].bit_count(), supports[second].bit_count()
    ):
        raise ValueError("mask sizes admit the exceptional reciprocal pair")

    for vertex, entries in (
        (first, certificate.first_pin_components),
        (second, certificate.second_pin_components),
    ):
        expected = tuple(bits(supports[vertex] ^ shared_bit))
        if tuple(entry[0] for entry in entries) != expected:
            raise ValueError("pin certificates do not cover exactly k-1 bits")
        for coordinate, component, conflict in entries:
            local = coordinate_graph(required, supports, coordinate)
            components = dict(nonbipartite_components(local))
            if vertex not in component or components.get(component) != conflict:
                raise ValueError("bad nonbipartite component certificate")


def certificate_json(certificate: OneFreeEdgeCertificate) -> dict:
    def pins_json(entries: tuple[PinComponent, ...]) -> list[dict]:
        return [
            {
                "coordinate": coordinate,
                "component": list(component),
                "conflict_edge": list(conflict),
            }
            for coordinate, component, conflict in entries
        ]

    return {
        "first": certificate.first,
        "second": certificate.second,
        "shared_coordinate": certificate.shared_coordinate,
        "first_mask": certificate.first_mask,
        "second_mask": certificate.second_mask,
        "first_pin_components": pins_json(certificate.first_pin_components),
        "second_pin_components": pins_json(certificate.second_pin_components),
    }


def certificate_from_json(value: dict) -> OneFreeEdgeCertificate:
    def pins(entries: list[dict]) -> tuple[PinComponent, ...]:
        return tuple(
            (
                int(entry["coordinate"]),
                tuple(map(int, entry["component"])),
                tuple(map(int, entry["conflict_edge"])),
            )
            for entry in entries
        )

    return OneFreeEdgeCertificate(
        first=int(value["first"]),
        second=int(value["second"]),
        shared_coordinate=int(value["shared_coordinate"]),
        first_mask=int(value["first_mask"]),
        second_mask=int(value["second_mask"]),
        first_pin_components=pins(value["first_pin_components"]),
        second_pin_components=pins(value["second_pin_components"]),
    )
