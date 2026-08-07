#!/usr/bin/env python3
"""Exact K7 one-free-edge obstruction with correlated pinned products.

This strictly extends the singleton-overlap lemma in
``d6_k7_one_free_edge.py``.  Each endpoint of a required edge again has one
unpinned coordinate, the same coordinate at both ends, but the propagated
masks may also overlap in pinned coordinates.  Their edge equation is

    F(k,P) F(l,Q) + S = 1,

where ``P,Q`` are the sums of the pinned endpoint signs and ``S`` is the sum
of their products on pinned shared coordinates.  Signs are variables of
nonbipartite coordinate components: occurrences in the same component share
one sign, so the finite enumeration retains all locally known correlations.
If no assignment satisfies the equation, the branch is impossible.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

import d6_k7_one_free_edge as base


PinComponent = base.PinComponent
ComponentKey = tuple[int, tuple[int, ...]]


@dataclass(frozen=True)
class CompatibilityWitness:
    component_signs: tuple[int, ...]
    first_sign_sum: int
    second_sign_sum: int
    pinned_shared_product_sum: int


@dataclass(frozen=True)
class CompatibilityResult:
    feasible: bool
    sign_variables: int
    assignments_checked: int
    first_witness: CompatibilityWitness | None


def component_key(entry: PinComponent) -> ComponentKey:
    coordinate, component, _conflict = entry
    return coordinate, component


def compatibility(
    first_mask: int,
    second_mask: int,
    free_coordinate: int,
    first_pins: Sequence[PinComponent],
    second_pins: Sequence[PinComponent],
) -> CompatibilityResult:
    """Exhaust the locally correlated pinned signs exactly."""

    free_bit = 1 << free_coordinate
    if (
        not (first_mask & free_bit)
        or not (second_mask & free_bit)
        or first_mask <= 0
        or second_mask <= 0
        or first_mask >= 128
        or second_mask >= 128
    ):
        raise ValueError("bad masks or free coordinate")
    first_expected = tuple(base.bits(first_mask ^ free_bit))
    second_expected = tuple(base.bits(second_mask ^ free_bit))
    if tuple(entry[0] for entry in first_pins) != first_expected:
        raise ValueError("first pins do not cover exactly the nonfree mask")
    if tuple(entry[0] for entry in second_pins) != second_expected:
        raise ValueError("second pins do not cover exactly the nonfree mask")

    entries = tuple(first_pins) + tuple(second_pins)
    keys = tuple(sorted({component_key(entry) for entry in entries}))
    key_index = {key: position for position, key in enumerate(keys)}
    first_keys = {
        entry[0]: key_index[component_key(entry)] for entry in first_pins
    }
    second_keys = {
        entry[0]: key_index[component_key(entry)] for entry in second_pins
    }
    shared_pinned = (first_mask & second_mask) ^ free_bit
    if not shared_pinned:
        raise ValueError("correlated layer requires a pinned shared coordinate")
    assignments = 1 << len(keys)
    one = (Fraction(1), Fraction(0))
    for assignment in range(assignments):
        signs = tuple(
            1 if assignment & (1 << position) else -1
            for position in range(len(keys))
        )
        first_sum = sum(signs[position] for position in first_keys.values())
        second_sum = sum(signs[position] for position in second_keys.values())
        shared_sum = sum(
            signs[first_keys[coordinate]] * signs[second_keys[coordinate]]
            for coordinate in base.bits(shared_pinned)
        )
        product = base.quadratic_product(
            base.one_free_value(first_mask.bit_count(), first_sum),
            base.one_free_value(second_mask.bit_count(), second_sum),
        )
        if (product[0] + shared_sum, product[1]) == one:
            return CompatibilityResult(
                feasible=True,
                sign_variables=len(keys),
                assignments_checked=assignment + 1,
                first_witness=CompatibilityWitness(
                    component_signs=signs,
                    first_sign_sum=first_sum,
                    second_sign_sum=second_sum,
                    pinned_shared_product_sum=shared_sum,
                ),
            )
    return CompatibilityResult(
        feasible=False,
        sign_variables=len(keys),
        assignments_checked=assignments,
        first_witness=None,
    )


@dataclass(frozen=True)
class CorrelatedOneFreeEdgeCertificate:
    first: int
    second: int
    free_coordinate: int
    first_mask: int
    second_mask: int
    first_pin_components: tuple[PinComponent, ...]
    second_pin_components: tuple[PinComponent, ...]
    sign_variables: int
    assignments_checked: int


def find_correlated_one_free_edge(
    graph: Sequence[int], masks: Sequence[int]
) -> CorrelatedOneFreeEdgeCertificate | None:
    required, supports = base.validate(graph, masks)
    pinned, data = base.pinning_data(required, supports)
    for second in range(len(required)):
        for first in base.bits(required[second] & ((1 << second) - 1)):
            first_free = supports[first] & ~pinned[first]
            second_free = supports[second] & ~pinned[second]
            if (
                first_free != second_free
                or first_free.bit_count() != 1
                or not ((supports[first] & supports[second]) & first_free)
            ):
                continue
            shared_pinned = (supports[first] & supports[second]) ^ first_free
            if not shared_pinned:
                continue
            free_coordinate = first_free.bit_length() - 1
            first_pins = tuple(
                data[(first, coordinate)]
                for coordinate in base.bits(supports[first] ^ first_free)
            )
            second_pins = tuple(
                data[(second, coordinate)]
                for coordinate in base.bits(supports[second] ^ second_free)
            )
            result = compatibility(
                supports[first],
                supports[second],
                free_coordinate,
                first_pins,
                second_pins,
            )
            if result.feasible:
                continue
            return CorrelatedOneFreeEdgeCertificate(
                first=first,
                second=second,
                free_coordinate=free_coordinate,
                first_mask=supports[first],
                second_mask=supports[second],
                first_pin_components=first_pins,
                second_pin_components=second_pins,
                sign_variables=result.sign_variables,
                assignments_checked=result.assignments_checked,
            )
    return None


def verify_certificate(
    graph: Sequence[int],
    masks: Sequence[int],
    certificate: CorrelatedOneFreeEdgeCertificate,
) -> None:
    required, supports = base.validate(graph, masks)
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
    free_bit = 1 << certificate.free_coordinate
    if (
        not (supports[first] & free_bit)
        or not (supports[second] & free_bit)
        or not ((supports[first] & supports[second]) ^ free_bit)
    ):
        raise ValueError("bad free coordinate or no pinned overlap")

    for vertex, entries in (
        (first, certificate.first_pin_components),
        (second, certificate.second_pin_components),
    ):
        expected = tuple(base.bits(supports[vertex] ^ free_bit))
        if tuple(entry[0] for entry in entries) != expected:
            raise ValueError("pin certificates do not cover exactly k-1 bits")
        for coordinate, component, conflict in entries:
            local = base.coordinate_graph(required, supports, coordinate)
            components = dict(base.nonbipartite_components(local))
            if vertex not in component or components.get(component) != conflict:
                raise ValueError("bad nonbipartite component certificate")

    result = compatibility(
        supports[first],
        supports[second],
        certificate.free_coordinate,
        certificate.first_pin_components,
        certificate.second_pin_components,
    )
    if result.feasible:
        raise ValueError("recorded edge has a compatible sign assignment")
    if (
        result.sign_variables != certificate.sign_variables
        or result.assignments_checked != certificate.assignments_checked
    ):
        raise ValueError("recorded exhaustive-enumeration counts differ")


def certificate_json(certificate: CorrelatedOneFreeEdgeCertificate) -> dict:
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
        "free_coordinate": certificate.free_coordinate,
        "first_mask": certificate.first_mask,
        "second_mask": certificate.second_mask,
        "first_pin_components": pins_json(certificate.first_pin_components),
        "second_pin_components": pins_json(certificate.second_pin_components),
        "sign_variables": certificate.sign_variables,
        "assignments_checked": certificate.assignments_checked,
    }


def certificate_from_json(value: dict) -> CorrelatedOneFreeEdgeCertificate:
    def pins(entries: list[dict]) -> tuple[PinComponent, ...]:
        return tuple(
            (
                int(entry["coordinate"]),
                tuple(map(int, entry["component"])),
                tuple(map(int, entry["conflict_edge"])),
            )
            for entry in entries
        )

    return CorrelatedOneFreeEdgeCertificate(
        first=int(value["first"]),
        second=int(value["second"]),
        free_coordinate=int(value["free_coordinate"]),
        first_mask=int(value["first_mask"]),
        second_mask=int(value["second_mask"]),
        first_pin_components=pins(value["first_pin_components"]),
        second_pin_components=pins(value["second_pin_components"]),
        sign_variables=int(value["sign_variables"]),
        assignments_checked=int(value["assignments_checked"]),
    )
