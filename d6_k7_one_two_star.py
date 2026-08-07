#!/usr/bin/env python3
"""Exact K7 obstruction from a two-free center and one-free neighbours.

After a K7 seed and a propagated support branch are fixed, normalize the
nonzero-factor defect vectors by ``w=sqrt(7)u/t``.  Nonbipartite components
of singleton-intersection required edges pin their coordinates to correlated
signs.  If a vertex of mask size ``k`` has exactly two unpinned coordinates
``x,y`` and pinned-sign sum ``P``, its diagonal identity is

    (x - (sqrt(7)-P)) (y - (sqrt(7)-P))
        = (k+4+P^2-2 sqrt(7)P)/2.

A required edge from this center to an exactly one-free neighbour becomes an
affine line in ``x,y`` because the neighbour's final coordinate is uniquely
determined in Q(sqrt(7)).  This module conjoins every such neighbour line and
tests its real intersection with the center hyperbola, exhaustively over the
correlated pinned-component signs.  Allowed unpinned coordinates may be
zero; candidate nonedges are never constrained.

The caller supplies only the nonzero-factor vertex set ``N`` of a fixed
zero-factor cover.  Thus every normalized endpoint used here has
``t=1+sum(u_i) != 0``; optional zero defect coordinates do not change that
cover classification.  The masks are propagated supersets of actual support:
singleton-overlap orthogonality deletes only a coordinate proved zero, and no
remaining allowed coordinate is assumed nonzero unless an odd component pins
it.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

import d6_k7_one_free_edge as base


Quadratic = tuple[Fraction, Fraction]
ComponentKey = tuple[int, tuple[int, ...]]
PinComponent = base.PinComponent
ZERO: Quadratic = (Fraction(0), Fraction(0))
ONE: Quadratic = (Fraction(1), Fraction(0))


def q(value: int | Fraction) -> Quadratic:
    return Fraction(value), Fraction(0)


def qadd(first: Quadratic, second: Quadratic) -> Quadratic:
    return first[0] + second[0], first[1] + second[1]


def qneg(value: Quadratic) -> Quadratic:
    return -value[0], -value[1]


def qsub(first: Quadratic, second: Quadratic) -> Quadratic:
    return qadd(first, qneg(second))


def qmul(first: Quadratic, second: Quadratic) -> Quadratic:
    return base.quadratic_product(first, second)


def qscale(value: Quadratic, scalar: int | Fraction) -> Quadratic:
    factor = Fraction(scalar)
    return value[0] * factor, value[1] * factor


def qinv(value: Quadratic) -> Quadratic:
    a, b = value
    denominator = a * a - 7 * b * b
    if denominator == 0:
        raise ZeroDivisionError("zero Q(sqrt(7)) denominator")
    return a / denominator, -b / denominator


def qdiv(first: Quadratic, second: Quadratic) -> Quadratic:
    return qmul(first, qinv(second))


def qsign(value: Quadratic) -> int:
    """Return the exact sign under the positive embedding of sqrt(7)."""

    a, b = value
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if (a > 0) == (b > 0):
        return 1 if a > 0 else -1
    comparison = a * a - 7 * b * b
    if comparison == 0:
        raise AssertionError("sqrt(7) would be rational")
    if a > 0:
        return 1 if comparison > 0 else -1
    return -1 if comparison > 0 else 1


def component_key(entry: PinComponent) -> ComponentKey:
    return entry[0], entry[1]


def diagonal_two_free(
    mask_size: int, pinned_sum: int
) -> tuple[Quadratic, Quadratic]:
    """Return ``A,R`` for the two-free diagonal ``(x-A)(y-A)=R``."""

    if not 2 <= mask_size <= 7:
        raise ValueError("two-free mask size must lie in [2,7]")
    if pinned_sum not in range(-(mask_size - 2), mask_size - 1, 2):
        raise ValueError("pinned sign sum has wrong range or parity")
    return (
        (Fraction(-pinned_sum), Fraction(1)),
        (
            Fraction(mask_size + 4 + pinned_sum * pinned_sum, 2),
            Fraction(-pinned_sum),
        ),
    )


def assigned_signs(
    keys: Sequence[ComponentKey], assignment: int
) -> dict[ComponentKey, int]:
    if not 0 <= assignment < 1 << len(keys):
        raise ValueError("sign assignment out of range")
    return {
        key: 1 if assignment & (1 << position) else -1
        for position, key in enumerate(keys)
    }


def endpoint_pin_keys(
    entries: Sequence[PinComponent],
) -> dict[int, ComponentKey]:
    answer = {entry[0]: component_key(entry) for entry in entries}
    if len(answer) != len(entries):
        raise ValueError("duplicate pinned coordinate")
    return answer


def endpoint_values(
    mask: int,
    free_mask: int,
    pin_keys: dict[int, ComponentKey],
    signs: dict[ComponentKey, int],
) -> dict[int, Quadratic]:
    if free_mask.bit_count() != 1 or mask & free_mask != free_mask:
        raise ValueError("endpoint is not exactly one-free")
    if set(pin_keys) != set(base.bits(mask ^ free_mask)):
        raise ValueError("pins do not cover the nonfree mask")
    values = {
        coordinate: q(signs[key]) for coordinate, key in pin_keys.items()
    }
    free_coordinate = free_mask.bit_length() - 1
    pinned_sum = sum(signs[key] for key in pin_keys.values())
    values[free_coordinate] = base.one_free_value(mask.bit_count(), pinned_sum)
    return values


def one_two_edge_line(
    center: int,
    neighbour: int,
    masks: Sequence[int],
    free_masks: Sequence[int],
    pin_keys: dict[int, dict[int, ComponentKey]],
    signs: dict[ComponentKey, int],
) -> tuple[Quadratic, Quadratic, Quadratic]:
    """Return the required-edge equation ``L*x+M*y=D`` at the center."""

    if free_masks[center].bit_count() != 2:
        raise ValueError("center is not exactly two-free")
    neighbour_values = endpoint_values(
        masks[neighbour], free_masks[neighbour], pin_keys[neighbour], signs
    )
    free_coordinates = tuple(base.bits(free_masks[center]))
    coefficients = {coordinate: ZERO for coordinate in free_coordinates}
    constant = ZERO
    for coordinate in base.bits(masks[center] & masks[neighbour]):
        neighbour_value = neighbour_values[coordinate]
        if free_masks[center] & (1 << coordinate):
            coefficients[coordinate] = qadd(
                coefficients[coordinate], neighbour_value
            )
        else:
            constant = qadd(
                constant,
                qscale(neighbour_value, signs[pin_keys[center][coordinate]]),
            )
    return (
        coefficients[free_coordinates[0]],
        coefficients[free_coordinates[1]],
        qsub(ONE, constant),
    )


@dataclass(frozen=True)
class LineSystemDecision:
    feasible: bool
    reason: str
    rank: int
    discriminant: Quadratic | None = None
    solution: tuple[Quadratic, Quadratic] | None = None


def line_hyperbola_decision(
    left: Quadratic,
    right: Quadratic,
    target: Quadratic,
    center: Quadratic,
    radius_product: Quadratic,
) -> LineSystemDecision:
    """Classify one affine line against ``(x-A)(y-A)=R`` exactly."""

    if left == ZERO and right == ZERO:
        return LineSystemDecision(
            feasible=target == ZERO,
            reason=("rank0_free_hyperbola" if target == ZERO else "zero_line_inconsistent"),
            rank=0,
        )
    if left == ZERO or right == ZERO:
        coefficient = right if left == ZERO else left
        fixed_value = qdiv(target, coefficient)
        feasible = fixed_value != center or radius_product == ZERO
        return LineSystemDecision(
            feasible=feasible,
            reason=(
                "rank1_fixed_coordinate"
                if feasible
                else "rank1_axis_asymptote"
            ),
            rank=1,
        )
    shifted_target = qsub(target, qmul(qadd(left, right), center))
    discriminant = qsub(
        qmul(shifted_target, shifted_target),
        qscale(qmul(qmul(left, right), radius_product), 4),
    )
    sign = qsign(discriminant)
    return LineSystemDecision(
        feasible=sign >= 0,
        reason=(
            "rank1_negative_discriminant"
            if sign < 0
            else "rank1_tangent"
            if sign == 0
            else "rank1_secant"
        ),
        rank=1,
        discriminant=discriminant,
    )


def lines_hyperbola_decision(
    lines: Sequence[tuple[Quadratic, Quadratic, Quadratic]],
    center: Quadratic,
    radius_product: Quadratic,
) -> LineSystemDecision:
    """Classify the common real intersection of lines and one hyperbola."""

    nontrivial = []
    for left, right, target in lines:
        if left == ZERO and right == ZERO:
            if target != ZERO:
                return LineSystemDecision(False, "zero_line_inconsistent", 0)
            continue
        nontrivial.append((left, right, target))
    if not nontrivial:
        return LineSystemDecision(True, "rank0_free_hyperbola", 0)

    first_left, first_right, first_target = nontrivial[0]
    rank_two_solution = None
    for left, right, target in nontrivial[1:]:
        determinant = qsub(qmul(first_left, right), qmul(first_right, left))
        if determinant != ZERO and rank_two_solution is None:
            rank_two_solution = (
                qdiv(
                    qsub(qmul(first_target, right), qmul(first_right, target)),
                    determinant,
                ),
                qdiv(
                    qsub(qmul(first_left, target), qmul(first_target, left)),
                    determinant,
                ),
            )
            continue
        if determinant == ZERO and (
            qsub(qmul(first_left, target), qmul(first_target, left)) != ZERO
            or qsub(qmul(first_right, target), qmul(first_target, right)) != ZERO
        ):
            return LineSystemDecision(
                False,
                (
                    "rank2_inconsistent"
                    if rank_two_solution is not None
                    else "parallel_inconsistent"
                ),
                2 if rank_two_solution is not None else 1,
                solution=rank_two_solution,
            )

    if rank_two_solution is None:
        return line_hyperbola_decision(
            first_left,
            first_right,
            first_target,
            center,
            radius_product,
        )

    x_value, y_value = rank_two_solution
    for left, right, target in nontrivial:
        if qadd(qmul(left, x_value), qmul(right, y_value)) != target:
            return LineSystemDecision(
                False,
                "rank2_inconsistent",
                2,
                solution=rank_two_solution,
            )
    on_diagonal = (
        qmul(qsub(x_value, center), qsub(y_value, center))
        == radius_product
    )
    return LineSystemDecision(
        feasible=on_diagonal,
        reason=("rank2_on_diagonal" if on_diagonal else "rank2_diagonal_mismatch"),
        rank=2,
        solution=rank_two_solution,
    )


@dataclass(frozen=True)
class AssignmentFailure:
    assignment: int
    reason: str
    rank: int
    discriminant: Quadratic | None
    solution: tuple[Quadratic, Quadratic] | None


@dataclass(frozen=True)
class StarCompatibility:
    feasible: bool
    sign_variables: int
    assignments_checked: int
    first_feasible_assignment: int | None
    failures: tuple[AssignmentFailure, ...]


def compatibility(
    center: int,
    neighbours: Sequence[int],
    masks: Sequence[int],
    free_masks: Sequence[int],
    pin_components: dict[int, tuple[PinComponent, ...]],
) -> StarCompatibility:
    """Exhaust every correlated pinned-component sign assignment."""

    required, supports = base.validate([0] * len(masks), masks)
    del required
    if not 0 <= center < len(supports):
        raise ValueError("center out of range")
    neighbour_tuple = tuple(map(int, neighbours))
    if (
        len(neighbour_tuple) < 2
        or len(set(neighbour_tuple)) != len(neighbour_tuple)
        or tuple(sorted(neighbour_tuple)) != neighbour_tuple
        or center in neighbour_tuple
    ):
        raise ValueError("neighbours must be a canonical distinct tuple")
    if free_masks[center].bit_count() != 2:
        raise ValueError("center is not exactly two-free")
    if any(free_masks[vertex].bit_count() != 1 for vertex in neighbour_tuple):
        raise ValueError("a neighbour is not exactly one-free")
    relevant = (center, *neighbour_tuple)
    pin_keys = {
        vertex: endpoint_pin_keys(pin_components[vertex]) for vertex in relevant
    }
    for vertex in relevant:
        if set(pin_keys[vertex]) != set(
            base.bits(supports[vertex] ^ free_masks[vertex])
        ):
            raise ValueError("pin components do not cover the nonfree mask")
    keys = tuple(
        sorted(
            {
                key
                for vertex in relevant
                for key in pin_keys[vertex].values()
            }
        )
    )
    failures = []
    for assignment in range(1 << len(keys)):
        signs = assigned_signs(keys, assignment)
        lines = tuple(
            one_two_edge_line(
                center,
                neighbour,
                supports,
                free_masks,
                pin_keys,
                signs,
            )
            for neighbour in neighbour_tuple
        )
        pinned_sum = sum(signs[key] for key in pin_keys[center].values())
        diagonal_center, radius_product = diagonal_two_free(
            supports[center].bit_count(), pinned_sum
        )
        decision = lines_hyperbola_decision(
            lines, diagonal_center, radius_product
        )
        if decision.feasible:
            return StarCompatibility(
                feasible=True,
                sign_variables=len(keys),
                assignments_checked=assignment + 1,
                first_feasible_assignment=assignment,
                failures=tuple(failures),
            )
        failures.append(
            AssignmentFailure(
                assignment=assignment,
                reason=decision.reason,
                rank=decision.rank,
                discriminant=decision.discriminant,
                solution=decision.solution,
            )
        )
    return StarCompatibility(
        feasible=False,
        sign_variables=len(keys),
        assignments_checked=1 << len(keys),
        first_feasible_assignment=None,
        failures=tuple(failures),
    )


@dataclass(frozen=True)
class OneTwoStarCertificate:
    center: int
    center_mask: int
    center_free_mask: int
    one_free_neighbours: tuple[int, ...]
    pin_components: tuple[tuple[int, tuple[PinComponent, ...]], ...]
    sign_variables: int
    assignments_checked: int
    failures: tuple[AssignmentFailure, ...]


def find_one_two_star(
    graph: Sequence[int], masks: Sequence[int]
) -> OneTwoStarCertificate | None:
    required, supports = base.validate(graph, masks)
    pinned, data = base.pinning_data(required, supports)
    free_masks = tuple(
        mask & ~pinned_mask for mask, pinned_mask in zip(supports, pinned)
    )
    for center, free_mask in enumerate(free_masks):
        if free_mask.bit_count() != 2:
            continue
        neighbours = tuple(
            neighbour
            for neighbour in base.bits(required[center])
            if free_masks[neighbour].bit_count() == 1
        )
        if len(neighbours) < 2:
            continue
        relevant = (center, *neighbours)
        components = {
            vertex: tuple(
                data[(vertex, coordinate)]
                for coordinate in base.bits(supports[vertex] ^ free_masks[vertex])
            )
            for vertex in relevant
        }
        result = compatibility(
            center, neighbours, supports, free_masks, components
        )
        if result.feasible:
            continue
        return OneTwoStarCertificate(
            center=center,
            center_mask=supports[center],
            center_free_mask=free_mask,
            one_free_neighbours=neighbours,
            pin_components=tuple(
                (vertex, components[vertex]) for vertex in relevant
            ),
            sign_variables=result.sign_variables,
            assignments_checked=result.assignments_checked,
            failures=result.failures,
        )
    return None


def verify_certificate(
    graph: Sequence[int],
    masks: Sequence[int],
    certificate: OneTwoStarCertificate,
) -> None:
    required, supports = base.validate(graph, masks)
    pinned, data = base.pinning_data(required, supports)
    free_masks = tuple(
        mask & ~pinned_mask for mask, pinned_mask in zip(supports, pinned)
    )
    center = int(certificate.center)
    if not 0 <= center < len(required):
        raise ValueError("center out of range")
    expected_neighbours = tuple(
        neighbour
        for neighbour in base.bits(required[center])
        if free_masks[neighbour].bit_count() == 1
    )
    if (
        supports[center] != certificate.center_mask
        or free_masks[center] != certificate.center_free_mask
        or free_masks[center].bit_count() != 2
        or expected_neighbours != certificate.one_free_neighbours
        or len(expected_neighbours) < 2
    ):
        raise ValueError("recorded star structure differs")
    relevant = (center, *expected_neighbours)
    recorded_components = dict(certificate.pin_components)
    if tuple(recorded_components) != relevant:
        raise ValueError("pin-component vertex order differs")
    for vertex in relevant:
        expected = tuple(
            data[(vertex, coordinate)]
            for coordinate in base.bits(supports[vertex] ^ free_masks[vertex])
        )
        if recorded_components[vertex] != expected:
            raise ValueError("pin-component certificate differs")
        for coordinate, component, conflict in expected:
            local = base.coordinate_graph(required, supports, coordinate)
            components = dict(base.nonbipartite_components(local))
            if vertex not in component or components.get(component) != conflict:
                raise ValueError("invalid nonbipartite component certificate")
    result = compatibility(
        center,
        expected_neighbours,
        supports,
        free_masks,
        recorded_components,
    )
    if result.feasible:
        raise ValueError("recorded star has a feasible sign assignment")
    if (
        result.sign_variables != certificate.sign_variables
        or result.assignments_checked != certificate.assignments_checked
        or result.failures != certificate.failures
    ):
        raise ValueError("exhaustive sign replay differs")


def quadratic_json(value: Quadratic) -> list[list[int]]:
    return [
        [value[0].numerator, value[0].denominator],
        [value[1].numerator, value[1].denominator],
    ]


def quadratic_from_json(value: Sequence[Sequence[int]]) -> Quadratic:
    if len(value) != 2 or any(len(coefficient) != 2 for coefficient in value):
        raise ValueError("bad quadratic-field encoding")
    return (
        Fraction(int(value[0][0]), int(value[0][1])),
        Fraction(int(value[1][0]), int(value[1][1])),
    )


def failure_json(failure: AssignmentFailure) -> dict:
    return {
        "assignment": failure.assignment,
        "reason": failure.reason,
        "rank": failure.rank,
        "discriminant": (
            None
            if failure.discriminant is None
            else quadratic_json(failure.discriminant)
        ),
        "solution": (
            None
            if failure.solution is None
            else [quadratic_json(value) for value in failure.solution]
        ),
    }


def failure_from_json(value: dict) -> AssignmentFailure:
    return AssignmentFailure(
        assignment=int(value["assignment"]),
        reason=str(value["reason"]),
        rank=int(value["rank"]),
        discriminant=(
            None
            if value["discriminant"] is None
            else quadratic_from_json(value["discriminant"])
        ),
        solution=(
            None
            if value["solution"] is None
            else tuple(quadratic_from_json(item) for item in value["solution"])
        ),
    )


def certificate_json(certificate: OneTwoStarCertificate) -> dict:
    return {
        "center": certificate.center,
        "center_mask": certificate.center_mask,
        "center_free_mask": certificate.center_free_mask,
        "one_free_neighbours": list(certificate.one_free_neighbours),
        "pin_components": [
            {
                "vertex": vertex,
                "entries": [
                    {
                        "coordinate": coordinate,
                        "component": list(component),
                        "conflict_edge": list(conflict),
                    }
                    for coordinate, component, conflict in entries
                ],
            }
            for vertex, entries in certificate.pin_components
        ],
        "sign_variables": certificate.sign_variables,
        "assignments_checked": certificate.assignments_checked,
        "failures": [failure_json(failure) for failure in certificate.failures],
    }


def certificate_from_json(value: dict) -> OneTwoStarCertificate:
    return OneTwoStarCertificate(
        center=int(value["center"]),
        center_mask=int(value["center_mask"]),
        center_free_mask=int(value["center_free_mask"]),
        one_free_neighbours=tuple(map(int, value["one_free_neighbours"])),
        pin_components=tuple(
            (
                int(item["vertex"]),
                tuple(
                    (
                        int(entry["coordinate"]),
                        tuple(map(int, entry["component"])),
                        tuple(map(int, entry["conflict_edge"])),
                    )
                    for entry in item["entries"]
                ),
            )
            for item in value["pin_components"]
        ),
        sign_variables=int(value["sign_variables"]),
        assignments_checked=int(value["assignments_checked"]),
        failures=tuple(failure_from_json(item) for item in value["failures"]),
    )
