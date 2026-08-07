#!/usr/bin/env python3
"""Focused exact controls for the K7 support-capacity relaxation."""

from __future__ import annotations

import itertools
import random
import unittest

import d6_k7_support_capacity as capacity


def graph_from_edges(order: int, edges: tuple[tuple[int, int], ...]):
    adjacency = [0] * order
    for first, second in edges:
        adjacency[first] |= 1 << second
        adjacency[second] |= 1 << first
    return tuple(adjacency)


def brute_force(graph, allowed, fixed=()):
    """Tiny independent product enumeration used only by controls."""

    domains = []
    for mask in allowed:
        domains.append(tuple(
            support for support in range(1, 128)
            if not support & ~mask
        ))
    for assignment in itertools.product(*domains):
        try:
            capacity.verify_witness(graph, allowed, fixed, assignment)
        except ValueError:
            continue
        return assignment
    return None


class CapacityKernelControls(unittest.TestCase):
    def test_feasible_witness_obeys_edges_and_caps(self):
        graph = graph_from_edges(3, ((0, 1), (1, 2)))
        allowed = (0b0011, 0b0110, 0b1100)
        result = capacity.solve_support_capacity(graph, allowed, node_limit=None)
        self.assertEqual(result.status, capacity.STATUS_FEASIBLE)
        self.assertIsNotNone(result.witness)
        capacity.verify_witness(graph, allowed, (), result.witness)

    def test_fixed_supports_consume_capacity(self):
        graph = graph_from_edges(2, ())
        # The two N vertices have only type {0}; a fixed point already uses
        # its cardinality-one capacity.
        result = capacity.solve_support_capacity(
            graph, (0b1, 0b1), (0b1,), node_limit=None
        )
        self.assertEqual(result.status, capacity.STATUS_INFEASIBLE)
        self.assertGreaterEqual(
            result.empty_domain_prunes + result.flow_prunes, 1
        )

    def test_fixed_family_can_already_exceed_a_type_cap(self):
        result = capacity.solve_support_capacity((), (), (0b1, 0b1))
        self.assertEqual(result.status, capacity.STATUS_INFEASIBLE)
        self.assertEqual(result.reason, "fixed_support_capacity")

    def test_required_edge_intersection_is_not_imposed_on_nonedge(self):
        disjoint = (0b1, 0b10)
        edge = capacity.solve_support_capacity(
            graph_from_edges(2, ((0, 1),)), disjoint, node_limit=None
        )
        nonedge = capacity.solve_support_capacity(
            graph_from_edges(2, ()), disjoint, node_limit=None
        )
        self.assertEqual(edge.status, capacity.STATUS_INFEASIBLE)
        self.assertEqual(nonedge.status, capacity.STATUS_FEASIBLE)

    def test_exact_type_cap_is_cardinality(self):
        # Type {0,1} has capacity two, so two copies pass and three fail.
        two = capacity.solve_support_capacity(
            graph_from_edges(2, ()), (0b11, 0b11), node_limit=None
        )
        three = capacity.solve_support_capacity(
            graph_from_edges(3, ()), (0b11, 0b11, 0b11),
            # Consume all proper subtypes to force exact type {0,1}.
            (0b1, 0b10), node_limit=None,
        )
        self.assertEqual(two.status, capacity.STATUS_FEASIBLE)
        self.assertEqual(three.status, capacity.STATUS_INFEASIBLE)

    def test_node_limit_is_an_explicit_nonclaim(self):
        result = capacity.solve_support_capacity(
            graph_from_edges(1, ()), (0b111,), node_limit=0
        )
        self.assertEqual(result.status, capacity.STATUS_UNRESOLVED)
        self.assertEqual(result.reason, "node_limit")

    def test_random_tiny_instances_match_full_product_search(self):
        rng = random.Random(0xD6C7)
        for _ in range(120):
            order = rng.randrange(0, 6)
            edges = tuple(
                (first, second)
                for first in range(order)
                for second in range(first)
                if rng.randrange(2)
            )
            graph = graph_from_edges(order, edges)
            # Keep masks tiny enough for independent Cartesian enumeration.
            allowed = tuple(rng.randrange(1, 8) for _ in range(order))
            fixed = tuple(rng.randrange(1, 8) for _ in range(rng.randrange(3)))
            expected = brute_force(graph, allowed, fixed)
            observed = capacity.solve_support_capacity(
                graph, allowed, fixed, node_limit=None
            )
            self.assertEqual(
                observed.status == capacity.STATUS_FEASIBLE,
                expected is not None,
                (graph, allowed, fixed, observed),
            )


if __name__ == "__main__":
    unittest.main()
