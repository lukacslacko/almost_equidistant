#!/usr/bin/env python3
"""Exact defect-support multiplicity relaxation after a required K7.

For every nonzero-factor outside vertex ``x`` the frozen singleton
propagation layer supplies a nonempty seven-bit mask ``E_x`` containing its
unknown *actual* defect support ``S_x``.  This module searches the finite
relaxation

    empty != S_x subseteq E_x,

subject to intersection on required edges and the exact multiplicity bound

    number of outside points of exact type S <= |S|.

Already fixed zero-factor supports consume the same capacities.  Candidate
nonedges impose no condition.  A proved infeasibility is therefore a sound
obstruction; feasibility is only a witness for this deliberately incomplete
support relaxation.

The search uses capacity-aware backtracking.  At every node an exact
bipartite b-matching test (variables to support types) applies the residual
capacity Hall relaxation while deliberately forgetting required-edge
intersection.  Its failure is a sound prune.  A node limit produces the
explicit status ``UNRESOLVED`` and is never converted into a rejection.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


STATUS_FEASIBLE = "FEASIBLE"
STATUS_INFEASIBLE = "INFEASIBLE"
STATUS_UNRESOLVED = "UNRESOLVED"


def bits(mask: int):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def nonempty_submasks(mask: int) -> tuple[int, ...]:
    """Return every nonempty submask, larger supports first."""

    if mask <= 0 or mask >= 128:
        raise ValueError("support masks must be nonempty seven-bit masks")
    answer = []
    submask = mask
    while submask:
        answer.append(submask)
        submask = (submask - 1) & mask
    answer.sort(key=lambda support: (-support.bit_count(), support))
    return tuple(answer)


def validate_graph(graph: Sequence[int]) -> tuple[int, ...]:
    adjacency = tuple(int(row) for row in graph)
    full = (1 << len(adjacency)) - 1
    for vertex, row in enumerate(adjacency):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError("invalid required-edge graph row")
        for neighbour in bits(row):
            if not adjacency[neighbour] & (1 << vertex):
                raise ValueError("asymmetric required-edge graph")
    return adjacency


@dataclass(frozen=True)
class CapacitySearchResult:
    status: str
    witness: tuple[int, ...] | None
    nodes: int
    flow_checks: int
    flow_prunes: int
    empty_domain_prunes: int
    edge_prunes: int
    capacity_prunes: int
    maximum_depth: int
    reason: str | None

    @property
    def feasible(self) -> bool:
        return self.status == STATUS_FEASIBLE


def _capacity_matching(
    domains: Sequence[tuple[int, ...]], capacities: Sequence[int]
) -> bool:
    """Exact bipartite b-matching feasibility for a list-domain family.

    This is a standard augmenting-path matcher with support vertices having
    integral capacities.  It is used only as a relaxation: the caller has
    already filtered domains against assigned required neighbours, while
    edges among two unassigned variables are intentionally forgotten.
    """

    if not domains:
        return True
    holders: list[list[int]] = [[] for _ in range(128)]
    order = sorted(range(len(domains)), key=lambda v: (len(domains[v]), v))

    def augment(variable: int, seen_types: set[int], seen_variables: set[int]) -> bool:
        for support in domains[variable]:
            if support in seen_types or capacities[support] <= 0:
                continue
            seen_types.add(support)
            if len(holders[support]) < capacities[support]:
                holders[support].append(variable)
                return True
            # Re-route one current occupant through an alternating path.
            for position, occupant in enumerate(tuple(holders[support])):
                if occupant in seen_variables:
                    continue
                seen_variables.add(occupant)
                if augment(occupant, seen_types, seen_variables):
                    holders[support][position] = variable
                    return True
        return False

    for variable in order:
        if not augment(variable, set(), {variable}):
            return False
    return True


def verify_witness(
    graph: Sequence[int],
    allowed_masks: Sequence[int],
    fixed_supports: Sequence[int],
    witness: Sequence[int],
) -> None:
    """Independently validate the defining finite constraints of a witness."""

    adjacency = validate_graph(graph)
    allowed = tuple(int(mask) for mask in allowed_masks)
    fixed = tuple(int(mask) for mask in fixed_supports)
    assigned = tuple(int(mask) for mask in witness)
    if len(adjacency) != len(allowed) or len(assigned) != len(allowed):
        raise ValueError("graph, mask, and witness orders disagree")
    if any(mask <= 0 or mask >= 128 for mask in (*allowed, *fixed, *assigned)):
        raise ValueError("all support masks must be nonempty and seven-bit")
    for support, mask in zip(assigned, allowed):
        if support & ~mask:
            raise ValueError("witness support is not contained in its mask")
    for first in range(len(adjacency)):
        for second in bits(adjacency[first] & ((1 << first) - 1)):
            if not assigned[first] & assigned[second]:
                raise ValueError("required-edge supports are disjoint")
    counts = [0] * 128
    for support in (*fixed, *assigned):
        counts[support] += 1
    for support in range(1, 128):
        if counts[support] > support.bit_count():
            raise ValueError("exact support type exceeds its multiplicity cap")


def solve_support_capacity(
    graph: Sequence[int],
    allowed_masks: Sequence[int],
    fixed_supports: Sequence[int] = (),
    *,
    node_limit: int | None = 200_000,
) -> CapacitySearchResult:
    """Decide the finite support-capacity relaxation, or report unresolved.

    ``INFEASIBLE`` is returned only after an exact exhaustive search or an
    exact Hall/capacity contradiction.  ``node_limit=None`` requests an
    unbounded exhaustive search.  A nonnegative finite limit is a resource
    guard and yields ``UNRESOLVED`` when exhausted.
    """

    adjacency = validate_graph(graph)
    allowed = tuple(int(mask) for mask in allowed_masks)
    fixed = tuple(int(mask) for mask in fixed_supports)
    if len(adjacency) != len(allowed):
        raise ValueError("graph and allowed-mask orders disagree")
    if any(mask <= 0 or mask >= 128 for mask in (*allowed, *fixed)):
        raise ValueError("all support masks must be nonempty and seven-bit")
    if node_limit is not None and node_limit < 0:
        raise ValueError("node limit must be nonnegative or None")

    capacities = [0] + [support.bit_count() for support in range(1, 128)]
    fixed_counts = [0] * 128
    for support in fixed:
        fixed_counts[support] += 1
        capacities[support] -= 1
        if capacities[support] < 0:
            return CapacitySearchResult(
                STATUS_INFEASIBLE, None, 0, 0, 0, 0, 0, 1, 0,
                "fixed_support_capacity",
            )

    base_domains = tuple(nonempty_submasks(mask) for mask in allowed)
    assignment = [0] * len(allowed)
    nodes = 0
    flow_checks = 0
    flow_prunes = 0
    empty_domain_prunes = 0
    edge_prunes = 0
    capacity_prunes = 0
    maximum_depth = 0
    exhausted = False

    def filtered_domain(vertex: int) -> tuple[int, ...]:
        nonlocal edge_prunes, capacity_prunes
        answer = []
        for support in base_domains[vertex]:
            if capacities[support] <= 0:
                capacity_prunes += 1
                continue
            compatible = True
            assigned_neighbours = adjacency[vertex]
            while assigned_neighbours:
                bit = assigned_neighbours & -assigned_neighbours
                assigned_neighbours ^= bit
                neighbour = bit.bit_length() - 1
                other = assignment[neighbour]
                if other and not support & other:
                    compatible = False
                    edge_prunes += 1
                    break
            if compatible:
                answer.append(support)
        return tuple(answer)

    def visit(depth: int) -> tuple[int, ...] | None:
        nonlocal nodes, flow_checks, flow_prunes, empty_domain_prunes
        nonlocal maximum_depth, exhausted
        maximum_depth = max(maximum_depth, depth)
        if depth == len(allowed):
            return tuple(assignment)
        if node_limit is not None and nodes >= node_limit:
            exhausted = True
            return None
        nodes += 1

        remaining = [v for v, support in enumerate(assignment) if not support]
        domains: dict[int, tuple[int, ...]] = {}
        for vertex in remaining:
            domain = filtered_domain(vertex)
            if not domain:
                empty_domain_prunes += 1
                return None
            domains[vertex] = domain

        flow_checks += 1
        flow_domains = tuple(domains[vertex] for vertex in remaining)
        if not _capacity_matching(flow_domains, capacities):
            flow_prunes += 1
            return None

        # Minimum remaining values, then most unassigned required neighbours.
        vertex = min(
            remaining,
            key=lambda v: (
                len(domains[v]),
                -sum(1 for u in bits(adjacency[v]) if not assignment[u]),
                v,
            ),
        )
        for support in domains[vertex]:
            assignment[vertex] = support
            capacities[support] -= 1
            witness = visit(depth + 1)
            capacities[support] += 1
            assignment[vertex] = 0
            if witness is not None:
                return witness
            if exhausted:
                return None
        return None

    witness = visit(0)
    if witness is not None:
        verify_witness(adjacency, allowed, fixed, witness)
        status = STATUS_FEASIBLE
        reason = None
    elif exhausted:
        status = STATUS_UNRESOLVED
        reason = "node_limit"
    else:
        status = STATUS_INFEASIBLE
        reason = "exhaustive_capacity_search"
    return CapacitySearchResult(
        status=status,
        witness=witness,
        nodes=nodes,
        flow_checks=flow_checks,
        flow_prunes=flow_prunes,
        empty_domain_prunes=empty_domain_prunes,
        edge_prunes=edge_prunes,
        capacity_prunes=capacity_prunes,
        maximum_depth=maximum_depth,
        reason=reason,
    )
