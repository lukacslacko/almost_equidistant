#!/usr/bin/env python3
"""Exact controls for the K6 bipartite Lorentz-component rank filter."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    _z0_subsets,
    check_bipartite_components,
    evaluate_rank_graph,
    lorentz_components,
    solve_rank_seed,
)
from d6_k6_lorentz_reference import K6LorentzInstance, add_edge
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "d6_k6_bipartite_rank_input.json"
REPORT = ROOT / "d6_k6_bipartite_rank_report.json"


def path_graph(n: int) -> tuple[int, ...]:
    adj = [0] * n
    for vertex in range(n - 1):
        add_edge(adj, vertex, vertex + 1)
    return tuple(adj)


def cycle_graph(n: int) -> tuple[int, ...]:
    adj = list(path_graph(n))
    add_edge(adj, 0, n - 1)
    return tuple(adj)


def complete_graph(n: int) -> tuple[int, ...]:
    return tuple(((1 << n) - 1) & ~(1 << vertex) for vertex in range(n))


def artificial_bipartite_instance() -> tuple[list[int], K6LorentzInstance]:
    """A kernel control whose initial bipartite component fails at Z0=empty.

    Its side graphs are P4 and P5, while L is complete bipartite.  Making the
    internal P4 vertex 7 a zero factor makes the dimension bound pass.  This
    synthetic object tests the Z0 quantifier and is not a graph-level alpha-2
    obstruction claim.
    """

    n = 15
    adj = [0] * n
    side_a = tuple(range(6, 10))
    side_b = tuple(range(10, 15))
    for left, right in zip(side_a, side_a[1:]):
        add_edge(adj, left, right)
    for left, right in zip(side_b, side_b[1:]):
        add_edge(adj, left, right)
    l_adj = [0] * 9
    for i in range(4):
        for j in range(4, 9):
            l_adj[i] |= 1 << j
            l_adj[j] |= 1 << i
            add_edge(adj, side_a[i], side_b[j - 4])
    # Only local vertex 1 (absolute vertex 7) is eligible for Z0.
    defects = [0] * 9
    defects[1] = 0b000111
    instance = K6LorentzInstance(
        tuple(range(6)),
        side_a + side_b,
        tuple(defects),
        tuple(l_adj),
        1 << 1,
    )
    return adj, instance


class K6BipartiteRankControls(unittest.TestCase):
    def test_zero_forcing_known_graphs(self) -> None:
        solver = ZeroForcingSolver()
        self.assertEqual(solver.solve(()).number, 0)
        self.assertEqual(solver.solve((0,)).number, 1)
        self.assertEqual(solver.solve(path_graph(6)).number, 1)
        self.assertEqual(solver.solve(cycle_graph(5)).number, 2)
        self.assertEqual(solver.solve(complete_graph(6)).number, 5)

    def test_component_bipartition_and_odd_cycle(self) -> None:
        # C4, C3, and one isolated vertex.
        l_adj = [0] * 8
        for u, v in ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 4)):
            add_edge(l_adj, u, v)
        instance = K6LorentzInstance(
            tuple(range(6)),
            tuple(range(6, 14)),
            tuple(0 for _ in range(8)),
            tuple(l_adj),
            0,
        )
        components = lorentz_components(instance, 0)
        self.assertEqual(len(components), 3)
        self.assertTrue(components[0].bipartite)
        self.assertEqual(components[0].side_a.bit_count(), 2)
        self.assertEqual(components[0].side_b.bit_count(), 2)
        self.assertFalse(components[1].bipartite)
        self.assertTrue(components[2].bipartite)
        self.assertEqual(components[2].side_a.bit_count(), 1)
        self.assertEqual(components[2].side_b, 0)

    def test_rank_inequality_detects_seven_dimensions(self) -> None:
        adj, instance = artificial_bipartite_instance()
        components = lorentz_components(instance, 0)
        passed, checks = check_bipartite_components(
            adj, instance, 0, components, ZeroForcingSolver()
        )
        self.assertFalse(passed)
        self.assertEqual(len(checks), 1)
        self.assertEqual(checks[0].zero_forcing_a, 1)
        self.assertEqual(checks[0].zero_forcing_b, 1)
        self.assertEqual(checks[0].required_dimension, 7)

    def test_zero_in_initially_bipartite_component_is_enumerated(self) -> None:
        adj, instance = artificial_bipartite_instance()
        self.assertEqual(list(_z0_subsets(instance.eligible_z0_mask)), [0, 0b10])
        self.assertTrue(lorentz_components(instance, 0)[0].bipartite)
        decision = solve_rank_seed(adj, instance)
        self.assertTrue(decision.feasible)
        self.assertEqual(decision.z0_subsets_considered, 2)
        self.assertEqual(decision.bipartite_rank_failures, 1)
        self.assertEqual(decision.chosen_z0, [7])

    def test_realizable_18_points_pass_over_every_k6_seed(self) -> None:
        decision = evaluate_rank_graph(lower_bound_18_graph(), scan_all=True)
        self.assertTrue(decision.applicable)
        self.assertFalse(decision.rejected)
        self.assertEqual(decision.seeds_checked, 32)
        self.assertEqual(decision.impossible_seeds, 0)

    def test_fixed_full_profile_witnesses(self) -> None:
        if not REPORT.exists():
            self.skipTest("run the complete bipartite-rank profile first")
        with REPORT.open(encoding="utf-8") as stream:
            report = json.load(stream)
        expected = {461_363}
        self.assertEqual(report["input_graphs"], 1_098)
        self.assertEqual(report["graphs_rejected"], 1)
        self.assertEqual(report["graphs_surviving"], 1_097)
        self.assertEqual(set(report["rejected_indices"]), expected)
        self.assertEqual(report["impossible_seeds"], 1)
        with INPUT.open(encoding="utf-8") as stream:
            sample = json.load(stream)
        selected = {
            record["index"]: record["adjacency"]
            for record in sample["graphs"]
            if record["index"] in expected
        }
        self.assertEqual(set(selected), expected)
        for index, adj in selected.items():
            with self.subTest(index=index):
                decision = evaluate_rank_graph(adj, scan_all=True)
                self.assertTrue(decision.rejected)
                self.assertGreater(decision.impossible_seeds, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
