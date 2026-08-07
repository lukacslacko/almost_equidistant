#!/usr/bin/env python3
"""Exact controls for the K7 three-disjoint-small-support refinement."""

from __future__ import annotations

import unittest
from functools import lru_cache
from itertools import combinations

import d6_k7_small_support_matching as matching
import d6_k7_small_support_value as frozen_value


def edge_mask(first: int, second: int) -> int:
    return (1 << first) | (1 << second)


def intersection_graph(supports: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for vertex, support in enumerate(supports):
        row = 0
        for neighbour, other in enumerate(supports):
            if vertex != neighbour and support & other:
                row |= 1 << neighbour
        rows.append(row)
    return tuple(rows)


def add_edges(
    graph: tuple[int, ...], edges: tuple[tuple[int, int], ...]
) -> tuple[int, ...]:
    rows = list(graph)
    for first, second in edges:
        rows[first] |= 1 << second
        rows[second] |= 1 << first
    return tuple(rows)


def independence_number(graph: tuple[int, ...]) -> int:
    answer = 0
    for mask in range(1 << len(graph)):
        if mask.bit_count() <= answer:
            continue
        if all(not (graph[vertex] & mask)
               for vertex in range(len(graph)) if mask & (1 << vertex)):
            answer = mask.bit_count()
    return answer


def independent_matching_number(
    edges: tuple[tuple[int, int], ...], vertex_count: int
) -> int:
    adjacency = [0] * vertex_count
    for first, second in edges:
        adjacency[first] |= 1 << second
        adjacency[second] |= 1 << first

    @lru_cache(maxsize=None)
    def solve(vertices: int) -> int:
        if not vertices:
            return 0
        first_bit = vertices & -vertices
        first = first_bit.bit_length() - 1
        remainder = vertices ^ first_bit
        answer = solve(remainder)
        neighbours = adjacency[first] & remainder
        while neighbours:
            second_bit = neighbours & -neighbours
            answer = max(answer, 1 + solve(remainder ^ second_bit))
            neighbours ^= second_bit
        return answer

    return solve((1 << vertex_count) - 1)


def is_forest(edges: tuple[tuple[int, int], ...], vertex_count: int) -> bool:
    parent = list(range(vertex_count))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for first, second in edges:
        first_root, second_root = root(first), root(second)
        if first_root == second_root:
            return False
        parent[first_root] = second_root
    return True


class DirectPackingControls(unittest.TestCase):
    def test_named_two_defect_patterns(self):
        p5 = tuple(edge_mask(i, i + 1) for i in range(4))
        p6 = tuple(edge_mask(i, i + 1) for i in range(5))
        three_k2 = tuple(edge_mask(2 * i, 2 * i + 1) for i in range(3))
        c6 = p6 + (edge_mask(0, 5),)

        self.assertIsNone(matching.refined_value_failure(p5))
        for supports in (p6, three_k2, c6):
            self.assertEqual(
                matching.refined_value_failure(supports),
                matching.PACKING_FAILURE,
            )

    def test_mixed_one_and_two_defect_packings(self):
        singleton_plus_two_k2 = (1, edge_mask(1, 2), edge_mask(3, 4))
        two_singletons_plus_k2 = (1, 2, edge_mask(2, 3))
        for supports in (singleton_plus_two_k2, two_singletons_plus_k2):
            self.assertIsNone(frozen_value.value_failure(supports))
            self.assertEqual(
                matching.refined_value_failure(supports),
                matching.PACKING_FAILURE,
            )

    def test_positions_form_a_checkable_witness(self):
        supports = (3, 6, 12, 24, 48)
        witness = matching.disjoint_support_triple(supports)
        self.assertIsNotNone(witness)
        selected = tuple(supports[position] for position in witness or ())
        self.assertEqual(len(selected), 3)
        self.assertTrue(all(a & b == 0 for a, b in combinations(selected, 2)))

    def test_frozen_failure_precedes_new_refinement(self):
        triangle = (edge_mask(0, 1), edge_mask(1, 2), edge_mask(0, 2))
        self.assertEqual(
            matching.refined_value_failure(triangle),
            "forbidden_two_defect_cycle",
        )

    def test_all_six_vertex_forests_match_independent_dp(self):
        possible = tuple(combinations(range(6), 2))
        for graph_mask in range(1 << len(possible)):
            edges = tuple(
                edge for position, edge in enumerate(possible)
                if graph_mask & (1 << position)
            )
            if not is_forest(edges, 6):
                continue
            supports = tuple(edge_mask(*edge) for edge in edges)
            expected_failure = independent_matching_number(edges, 6) >= 3
            observed_failure = (
                matching.disjoint_support_triple(supports) is not None
            )
            self.assertEqual(observed_failure, expected_failure, edges)


class RefinedCspControls(unittest.TestCase):
    def test_six_cycle_is_now_infeasible(self):
        masks = (3, 6, 12, 24, 48, 33)
        graph = intersection_graph(masks)
        self.assertEqual(independence_number(graph), 3)
        self.assertTrue(
            frozen_value.check_small_support_masks(graph, masks).feasible
        )
        result = matching.check_small_support_masks(graph, masks)
        self.assertFalse(result.feasible)
        self.assertGreater(dict(result.failures)[matching.PACKING_FAILURE], 0)

    def test_alpha_two_completion_is_already_rejected_by_frozen_csp(self):
        masks = (3, 6, 12, 24, 48, 33)
        # The only independent triples in the synthetic point-C6 are its two
        # parity classes.  One chord in each class makes alpha=2.  Both new
        # required edges join disjoint *allowed* masks, so even the frozen
        # support-intersection-only layer rejects before value constraints.
        graph = add_edges(intersection_graph(masks), ((0, 2), (1, 3)))
        self.assertEqual(independence_number(graph), 2)
        support_only = frozen_value.check_small_support_masks(
            graph, masks, apply_value_constraints=False
        )
        self.assertFalse(support_only.feasible)
        self.assertGreater(
            dict(support_only.failures)["disjoint_required_small_support"],
            0,
        )

    def test_matching_at_most_two_forest_is_retained(self):
        masks = (3, 6, 12, 24)  # P5, with matching number two.
        result = matching.check_small_support_masks(
            intersection_graph(masks), masks
        )
        self.assertTrue(result.feasible)

    def test_large_masks_remain_unassigned(self):
        result = matching.check_small_support_masks((0, 0), (7, 112))
        self.assertTrue(result.feasible)
        self.assertEqual(result.first_witness, ())


if __name__ == "__main__":
    unittest.main()
