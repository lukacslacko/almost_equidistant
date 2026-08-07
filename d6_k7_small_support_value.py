#!/usr/bin/env python3
"""Exact one/two-defect value constraints after K7 support propagation.

This module is deliberately a small, independent refinement.  It consumes a
required-edge graph on the nonzero-factor vertices and propagated seven-bit
support masks.  Vertices whose mask has size at most two must have an actual
one- or two-defect support.  We enumerate those finitely many support choices
and apply exact value consequences proved in
``d6_k7_small_support_value.md``.

Masks of size at least three are left completely unconstrained.  Thus a
``False`` result is a sound obstruction, while a ``True`` result merely says
that this deliberately incomplete screen found no obstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Sequence


def bits(mask: int):
    """Yield the set bits of a seven-bit support mask."""

    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def simple_cycle_lengths(simple_edges: set[tuple[int, int]]) -> set[int]:
    """Return all simple-cycle lengths in a graph on seed coordinates.

    There are only seven seed coordinates.  Exhausting vertex subsets and
    cyclic orders is smaller and easier to audit than a general graph-cycle
    implementation.  Parallel two-defect types have already been collapsed
    in ``simple_edges``.
    """

    adjacency = [[False] * 7 for _ in range(7)]
    used_vertices = set()
    for first, second in simple_edges:
        if not (0 <= first < second < 7):
            raise ValueError("bad two-defect support edge")
        adjacency[first][second] = adjacency[second][first] = True
        used_vertices.update((first, second))

    found: set[int] = set()
    ordered_vertices = sorted(used_vertices)
    for size in range(3, len(ordered_vertices) + 1):
        for chosen in combinations(ordered_vertices, size):
            # The least vertex is fixed first.  Requiring the second vertex
            # to be smaller than the last removes reversal duplicates.
            first = chosen[0]
            for tail in permutations(chosen[1:]):
                if tail[0] > tail[-1]:
                    continue
                cycle = (first, *tail)
                if all(
                    adjacency[cycle[i]][cycle[(i + 1) % size]]
                    for i in range(size)
                ):
                    found.add(size)
                    break
            if size in found:
                # We need only the set of lengths, not every witness.
                continue
    return found


def _components(simple_edges: set[tuple[int, int]]) -> tuple[set[int], ...]:
    adjacency = [set() for _ in range(7)]
    for first, second in simple_edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = {vertex for vertex in range(7) if adjacency[vertex]}
    answer = []
    while unseen:
        root = min(unseen)
        stack = [root]
        component = {root}
        unseen.remove(root)
        while stack:
            vertex = stack.pop()
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        answer.append(component)
    return tuple(answer)


def value_failure(supports: Sequence[int]) -> str | None:
    """Return an exact sparse-value obstruction for fixed actual supports.

    Every input support must have cardinality one or two.  Only consequences
    involving this selected subset of points are used, so the function is
    also sound on a partial assignment.
    """

    if any(mask <= 0 or mask >= 128 or mask.bit_count() not in (1, 2)
           for mask in supports):
        raise ValueError("value checker expects nonempty one/two-bit masks")

    singleton_coordinates = [next(bits(mask)) for mask in supports
                             if mask.bit_count() == 1]
    if len(singleton_coordinates) > 2:
        return "three_one_defects"
    if len(singleton_coordinates) != len(set(singleton_coordinates)):
        return "duplicate_one_defect"

    two_types = [tuple(bits(mask)) for mask in supports
                 if mask.bit_count() == 2]
    multiplicities: dict[tuple[int, int], int] = {}
    for edge in two_types:
        multiplicities[edge] = multiplicities.get(edge, 0) + 1
        if multiplicities[edge] > 2:
            return "three_same_two_defects"

    singleton_set = set(singleton_coordinates)
    incident_to_singleton = {coordinate: 0 for coordinate in singleton_set}
    for first, second in two_types:
        if first in incident_to_singleton:
            incident_to_singleton[first] += 1
        if second in incident_to_singleton:
            incident_to_singleton[second] += 1
    if any(count > 1 for count in incident_to_singleton.values()):
        return "two_defects_repeat_at_one_defect"

    simple_edges = set(multiplicities)
    cycle_lengths = simple_cycle_lengths(simple_edges)
    if any(length % 6 for length in cycle_lengths):
        return "forbidden_two_defect_cycle"

    degrees = [0] * 7
    for first, second in simple_edges:
        degrees[first] += 1
        degrees[second] += 1

    # If two points have the same two-defect type and a third point has a
    # distinct type incident at either endpoint, unit products with the third
    # point force the two copies to have an equal endpoint value.  The phi
    # relation then makes both normalized vectors equal, impossible for their
    # own required unit edge.
    if any(
        multiplicity == 2 and (degrees[first] >= 2 or degrees[second] >= 2)
        for (first, second), multiplicity in multiplicities.items()
    ):
        return "parallel_two_defect_type_touches_another"

    for component in _components(simple_edges):
        component_singletons = component & singleton_set
        if len(component_singletons) > 1:
            return "one_defects_joined_by_two_defects"
        if component_singletons and any(degrees[v] >= 3 for v in component):
            return "one_defect_component_branches"

    # At a coordinate incident with at least three distinct types, all their
    # endpoint values are a common sign sigma in {+1,-1}.  A maximal path of
    # length ell between two such branch coordinates transports the sign by
    # T^ell.  The exact order-six action maps a sign to a sign only for
    # ell=0 mod 3, preserving it for 0 mod 6 and reversing it for 3 mod 6.
    branch_vertices = {v for v in range(7) if degrees[v] >= 3}
    sign_constraints: list[tuple[int, int, int]] = []
    seen_paths: set[tuple[int, int, int]] = set()
    for start in sorted(branch_vertices):
        neighbours = {
            second if first == start else first
            for first, second in simple_edges
            if first == start or second == start
        }
        for neighbour in sorted(neighbours):
            previous = start
            current = neighbour
            length = 1
            while current not in branch_vertices and degrees[current] == 2:
                following = {
                    second if first == current else first
                    for first, second in simple_edges
                    if first == current or second == current
                } - {previous}
                if len(following) != 1:
                    raise AssertionError("bad degree-two path traversal")
                previous, current = current, next(iter(following))
                length += 1
            if current not in branch_vertices or current == start:
                continue
            key = (min(start, current), max(start, current), length)
            if key in seen_paths:
                continue
            seen_paths.add(key)
            if length % 3:
                return "incompatible_two_defect_branch_distance"
            parity = (length // 3) & 1
            sign_constraints.append((start, current, parity))

    signs: dict[int, int] = {}
    constraint_graph: dict[int, list[tuple[int, int]]] = {
        vertex: [] for vertex in branch_vertices
    }
    for first, second, parity in sign_constraints:
        constraint_graph[first].append((second, parity))
        constraint_graph[second].append((first, parity))
    for root in sorted(branch_vertices):
        if root in signs:
            continue
        signs[root] = 0
        stack = [root]
        while stack:
            vertex = stack.pop()
            for neighbour, parity in constraint_graph[vertex]:
                required = signs[vertex] ^ parity
                if neighbour in signs:
                    if signs[neighbour] != required:
                        return "inconsistent_two_defect_branch_signs"
                else:
                    signs[neighbour] = required
                    stack.append(neighbour)
    return None


@dataclass(frozen=True)
class SmallSupportResult:
    feasible: bool
    assignments_checked: int
    first_witness: tuple[int, ...] | None
    failures: tuple[tuple[str, int], ...]


def check_small_support_masks(
    graph_n: Sequence[int],
    propagated_masks: Sequence[int],
    *,
    apply_value_constraints: bool = True,
) -> SmallSupportResult:
    """Existentially assign every support forced to have size at most two.

    A required edge between two such vertices also requires their actual
    supports to intersect: their normalized dot product is one and therefore
    cannot have disjoint support.  Larger masks are intentionally omitted.
    """

    if len(graph_n) != len(propagated_masks):
        raise ValueError("graph and mask orders disagree")
    count = len(graph_n)
    full = (1 << count) - 1
    for vertex, row in enumerate(graph_n):
        if row & ~full or row & (1 << vertex):
            raise ValueError("bad graph adjacency")
        for neighbour in bits(row):
            if not (graph_n[neighbour] & (1 << vertex)):
                raise ValueError("asymmetric graph adjacency")
    if any(mask <= 0 or mask >= 128 for mask in propagated_masks):
        raise ValueError("propagated masks must be nonempty seven-bit masks")

    selected = [vertex for vertex, mask in enumerate(propagated_masks)
                if mask.bit_count() <= 2]
    if not selected:
        return SmallSupportResult(True, 1, (), ())

    domains: dict[int, tuple[int, ...]] = {}
    for vertex in selected:
        mask = int(propagated_masks[vertex])
        coordinates = tuple(1 << coordinate for coordinate in bits(mask))
        if len(coordinates) == 1:
            domains[vertex] = coordinates
        else:
            domains[vertex] = (coordinates[0], coordinates[1], mask)

    selected_set = set(selected)
    small_degrees = {
        vertex: sum(
            1 for neighbour in bits(graph_n[vertex])
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

    def visit(depth: int) -> tuple[int, ...] | None:
        nonlocal checked
        if depth == len(order):
            checked += 1
            ordered_supports = tuple(assignment[vertex] for vertex in order)
            failure = (
                value_failure(ordered_supports)
                if apply_value_constraints
                else None
            )
            if failure is None:
                return ordered_supports
            failure_counts[failure] = failure_counts.get(failure, 0) + 1
            return None

        vertex = order[depth]
        for support in domains[vertex]:
            # Every already assigned required neighbour must have an actual
            # support intersection.  This is necessary before any value
            # identity is considered.
            if any(
                graph_n[vertex] & (1 << prior)
                and not (support & prior_support)
                for prior, prior_support in assignment.items()
            ):
                name = "disjoint_required_small_support"
                failure_counts[name] = failure_counts.get(name, 0) + 1
                continue
            assignment[vertex] = support
            partial_supports = tuple(assignment[v] for v in order[:depth + 1])
            failure = (
                value_failure(partial_supports)
                if apply_value_constraints
                else None
            )
            if failure is None:
                witness = visit(depth + 1)
                if witness is not None:
                    return witness
            else:
                failure_counts[failure] = failure_counts.get(failure, 0) + 1
            del assignment[vertex]
        return None

    witness = visit(0)
    return SmallSupportResult(
        feasible=witness is not None,
        assignments_checked=checked,
        first_witness=witness,
        failures=tuple(sorted(failure_counts.items())),
    )
