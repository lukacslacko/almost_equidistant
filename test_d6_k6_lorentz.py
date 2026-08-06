#!/usr/bin/env python3
"""Exact controls for the independent K6 odd-component Lorentz CSP."""

from __future__ import annotations

import unittest

from d6_k6_lorentz_reference import (
    add_clique,
    add_edge,
    build_instance,
    evaluate_graph,
    feasible_light_ray_assignments,
    find_clique_mask,
    has_independent_triple,
    nonbipartite_components,
    solve_k6_seed,
    support_matching,
    validate_graph,
    vertices,
)


def lower_bound_18_graph() -> list[int]:
    """The realizable 16-point half-cube plus its two unit apices."""

    base = [word for word in range(32) if word.bit_count() % 2 == 1]
    adj = [0] * 18
    for i, word in enumerate(base):
        for j, other in enumerate(base[:i]):
            if (word ^ other).bit_count() == 2:
                add_edge(adj, i, j)
    for apex in (16, 17):
        for vertex in range(16):
            add_edge(adj, apex, vertex)
    validate_graph(adj)
    return adj


def oversized_odd_component() -> list[int]:
    """A K6-only alpha-two graph with an undeletable 7-vertex odd L component."""

    adj = [0] * 13
    add_clique(adj, range(6))
    # Local outside masks and nonedges are a compact isolating control supplied
    # by an independent audit.  Every defect mask has size at most two, so Z0
    # is empty.  L is connected/non-bipartite on all seven outside vertices,
    # hence too large for one six-dimensional lightlike bin.  The full graph
    # nevertheless has omega=6 and alpha<=2, so rejection is not a disguised
    # K7/K8 or independent-triple check.
    defect_masks = (
        {0},
        {0, 4},
        {1, 3},
        {4, 5},
        {3, 4},
        {2, 3},
        {2, 5},
    )
    outside_nonedges = {(0, 2), (0, 6), (3, 5), (4, 6)}
    for local, x in enumerate(range(6, 13)):
        for q in range(6):
            if q not in defect_masks[local]:
                add_edge(adj, x, q)
    for i in range(7):
        for j in range(i):
            if (j, i) not in outside_nonedges:
                add_edge(adj, 6 + i, 6 + j)
    validate_graph(adj)
    return adj


def opposite_color_triangles() -> list[int]:
    """Two odd L components whose masks fit only in opposite bins."""

    adj = [0] * 12
    add_clique(adj, range(6))
    # Both triangles have singleton masks {0},{1},{2}.  Each triangle is
    # matchable, whereas their union is not.  Cross edges join equal types;
    # those masks overlap, so the edges do not enter L.  The outside graph is
    # the triangular prism, whose complement is a six-cycle, hence alpha<=2.
    for triangle_start in (6, 9):
        for local, x in enumerate(range(triangle_start, triangle_start + 3)):
            for q in range(6):
                if q != local:
                    add_edge(adj, x, q)
        add_clique(adj, range(triangle_start, triangle_start + 3))
    for local in range(3):
        add_edge(adj, 6 + local, 9 + local)
    validate_graph(adj)
    return adj


def z0_matching_failure_graph() -> list[int]:
    """Four eligible Z0 candidates confined to only three coordinates."""

    adj = [0] * 10
    add_clique(adj, range(6))
    for x in range(6, 10):
        for q in range(3, 6):
            add_edge(adj, x, q)  # D_x={0,1,2}, so every x is Z0-eligible.
    add_clique(adj, range(6, 10))
    validate_graph(adj)
    return adj


class K6LorentzControls(unittest.TestCase):
    def test_realizable_18_points_pass_over_every_k6_seed(self) -> None:
        adj = lower_bound_18_graph()
        decision = evaluate_graph(adj, require_k6_only=True, scan_all=True)
        self.assertTrue(decision.applicable)
        self.assertFalse(decision.rejected)
        self.assertEqual(decision.seeds_checked, 32)
        self.assertEqual(decision.impossible_seeds, 0)

    def test_oversized_robust_odd_component_is_rejected(self) -> None:
        adj = oversized_odd_component()
        self.assertFalse(has_independent_triple(adj))
        self.assertFalse(find_clique_mask(adj, 7))
        self.assertTrue(find_clique_mask(adj, 6))
        instance = build_instance(adj, tuple(range(6)))
        self.assertEqual(instance.eligible_z0_mask, 0)
        odd = nonbipartite_components(instance, 0)
        self.assertEqual(len(odd), 1)
        self.assertEqual(odd[0].bit_count(), 7)
        decision = solve_k6_seed(adj, tuple(range(6)))
        self.assertFalse(decision.feasible)
        self.assertEqual(decision.z0_subsets_considered, 1)
        self.assertEqual(decision.colorings_considered, 2)
        graph_decision = evaluate_graph(adj, require_k6_only=True, scan_all=True)
        self.assertTrue(graph_decision.applicable)
        self.assertTrue(graph_decision.rejected)

    def test_two_odd_components_require_opposite_colors(self) -> None:
        adj = opposite_color_triangles()
        self.assertFalse(has_independent_triple(adj))
        instance = build_instance(adj, tuple(range(6)))
        self.assertEqual(instance.eligible_z0_mask, 0)
        odd = nonbipartite_components(instance, 0)
        self.assertEqual(
            sorted(component.bit_count() for component in odd), [3, 3]
        )
        assignments = feasible_light_ray_assignments(instance, 0, odd)
        self.assertEqual(
            {colors for colors, _ in assignments}, {(0, 1), (1, 0)}
        )
        decision = solve_k6_seed(adj, tuple(range(6)))
        self.assertTrue(decision.feasible)
        assert decision.component_colors is not None
        self.assertNotEqual(*decision.component_colors)

    def test_candidate_z0_fails_support_matching(self) -> None:
        adj = z0_matching_failure_graph()
        self.assertFalse(has_independent_triple(adj))
        instance = build_instance(adj, tuple(range(6)))
        all_outside = (1 << len(instance.outside)) - 1
        self.assertEqual(instance.eligible_z0_mask, all_outside)
        self.assertTrue(all(mask.bit_count() == 3 for mask in instance.defects))
        self.assertIsNone(support_matching(all_outside, instance.defects))
        # Smaller candidate zero sets remain possible, so this is an isolating
        # control of Z0 support matching rather than a graph rejection claim.
        self.assertTrue(solve_k6_seed(adj, tuple(range(6))).feasible)

    def test_alpha_two_precondition_is_enforced(self) -> None:
        adj = [0] * 9
        add_clique(adj, range(6))
        self.assertTrue(has_independent_triple(adj))
        with self.assertRaisesRegex(ValueError, "alpha"):
            solve_k6_seed(adj, tuple(range(6)))

    def test_candidate_nonedges_are_not_forced_nonunit(self) -> None:
        # The actual unit-edge graph may contain an omitted candidate edge.
        # Adding it only strengthens the required graph and is accepted by the
        # reference; nowhere does the CSP demand a nonzero defect/non-unit pair.
        adj = opposite_color_triangles()
        add_edge(adj, 6, 10)
        validate_graph(adj)
        decision = solve_k6_seed(adj, tuple(range(6)))
        self.assertIsInstance(decision.feasible, bool)


if __name__ == "__main__":
    unittest.main(verbosity=2)
