#!/usr/bin/env python3
"""Exact three-disjoint-support refinement after the frozen K7 value layer.

The frozen sparse-value checker deliberately retained a six-cycle of exact
two-defect types.  This independent post-layer adds the stronger and simpler
set-packing consequence proved in ``d6_k7_small_support_matching.md``:
among assigned actual one/two-defect supports, no three may be pairwise
disjoint.

Only vertices whose propagated mask has size at most two are assigned here.
Larger masks remain completely unconstrained.  Thus an infeasible result is
a sound partial-assignment obstruction; a feasible result is not a claimed
realization.
"""

from __future__ import annotations

from itertools import combinations
from typing import Sequence

import d6_k7_small_support_value as frozen_value


PACKING_FAILURE = "three_pairwise_disjoint_small_supports"


def disjoint_support_triple(
    supports: Sequence[int],
) -> tuple[int, int, int] | None:
    """Return positions of three pairwise-disjoint nonempty small supports.

    The input is an actual partial assignment, not a collection of allowed
    masks.  Supports of equal type remain separate points, hence positions
    rather than support values form the witness.
    """

    if any(mask <= 0 or mask >= 128 or mask.bit_count() not in (1, 2)
           for mask in supports):
        raise ValueError("packing checker expects nonempty one/two-bit masks")
    for first, second, third in combinations(range(len(supports)), 3):
        a, b, c = supports[first], supports[second], supports[third]
        if not (a & b or a & c or b & c):
            return first, second, third
    return None


def refined_value_failure(supports: Sequence[int]) -> str | None:
    """Apply the frozen exact rules and then the new set-packing rule."""

    failure = frozen_value.value_failure(supports)
    if failure is not None:
        return failure
    if disjoint_support_triple(supports) is not None:
        return PACKING_FAILURE
    return None


def check_small_support_masks(
    graph_n: Sequence[int],
    propagated_masks: Sequence[int],
) -> frozen_value.SmallSupportResult:
    """Existentially assign small supports with the exact packing refinement.

    This is a deliberately local copy of the frozen checker's small CSP,
    with ``refined_value_failure`` replacing its leaf/partial predicate.  It
    does not mutate or monkey-patch the frozen module.
    """

    if len(graph_n) != len(propagated_masks):
        raise ValueError("graph and mask orders disagree")
    count = len(graph_n)
    full = (1 << count) - 1
    for vertex, row in enumerate(graph_n):
        if row & ~full or row & (1 << vertex):
            raise ValueError("bad graph adjacency")
        for neighbour in frozen_value.bits(row):
            if not (graph_n[neighbour] & (1 << vertex)):
                raise ValueError("asymmetric graph adjacency")
    if any(mask <= 0 or mask >= 128 for mask in propagated_masks):
        raise ValueError("propagated masks must be nonempty seven-bit masks")

    selected = [
        vertex for vertex, mask in enumerate(propagated_masks)
        if mask.bit_count() <= 2
    ]
    if not selected:
        return frozen_value.SmallSupportResult(True, 1, (), ())

    domains: dict[int, tuple[int, ...]] = {}
    for vertex in selected:
        mask = int(propagated_masks[vertex])
        coordinates = tuple(1 << coordinate for coordinate in
                            frozen_value.bits(mask))
        if len(coordinates) == 1:
            domains[vertex] = coordinates
        else:
            domains[vertex] = (coordinates[0], coordinates[1], mask)

    selected_set = set(selected)
    small_degrees = {
        vertex: sum(
            1 for neighbour in frozen_value.bits(graph_n[vertex])
            if neighbour in selected_set
        )
        for vertex in selected
    }
    order = tuple(sorted(
        selected,
        key=lambda vertex: (
            len(domains[vertex]),
            -small_degrees[vertex],
            vertex,
        ),
    ))
    assignment: dict[int, int] = {}
    checked = 0
    failure_counts: dict[str, int] = {}

    def record(failure: str) -> None:
        failure_counts[failure] = failure_counts.get(failure, 0) + 1

    def visit(depth: int) -> tuple[int, ...] | None:
        nonlocal checked
        if depth == len(order):
            checked += 1
            ordered_supports = tuple(assignment[vertex] for vertex in order)
            failure = refined_value_failure(ordered_supports)
            if failure is None:
                return ordered_supports
            record(failure)
            return None

        vertex = order[depth]
        for support in domains[vertex]:
            # A required unit edge between two assigned small-support points
            # requires intersecting actual supports: the normalized dot
            # product must equal one and so cannot be zero.
            if any(
                graph_n[vertex] & (1 << prior)
                and not (support & prior_support)
                for prior, prior_support in assignment.items()
            ):
                record("disjoint_required_small_support")
                continue
            assignment[vertex] = support
            partial_supports = tuple(
                assignment[prior] for prior in order[:depth + 1]
            )
            failure = refined_value_failure(partial_supports)
            if failure is None:
                witness = visit(depth + 1)
                if witness is not None:
                    return witness
            else:
                record(failure)
            del assignment[vertex]
        return None

    witness = visit(0)
    return frozen_value.SmallSupportResult(
        feasible=witness is not None,
        assignments_checked=checked,
        first_witness=witness,
        failures=tuple(sorted(failure_counts.items())),
    )
