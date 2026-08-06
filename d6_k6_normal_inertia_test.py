#!/usr/bin/env python3
"""Exact controls for the K6 normal-coordinate inertia prototype."""

from __future__ import annotations

import json
import random
import unittest
from itertools import combinations
from pathlib import Path

from d6_k6_normal_inertia import (
    evaluate_normal_graph,
    exact_inertia_rows,
    normal_dimension_lower_bound,
)
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "d6_k6_bipartite_rank_input.json"
REPORT = ROOT / "d6_k6_normal_inertia_sample.json"


def i_plus_adjacency(size: int, edges: tuple[tuple[int, int], ...]) -> list[list[int]]:
    matrix = [[int(row == column) for column in range(size)] for row in range(size)]
    for left, right in edges:
        matrix[left][right] = 1
        matrix[right][left] = 1
    return matrix


class K6NormalInertiaControls(unittest.TestCase):
    def test_exact_inertia_known_unit_patterns(self) -> None:
        self.assertEqual(exact_inertia_rows([]), (0, 0, 0))
        self.assertEqual(exact_inertia_rows([[1, 0], [0, 1]]), (2, 0, 0))
        self.assertEqual(exact_inertia_rows([[1, 1], [1, 1]]), (1, 0, 1))
        path_three = i_plus_adjacency(3, ((0, 1), (1, 2)))
        cycle_four = i_plus_adjacency(
            4, ((0, 1), (1, 2), (2, 3), (3, 0))
        )
        self.assertEqual(exact_inertia_rows(path_three), (2, 1, 0))
        self.assertEqual(exact_inertia_rows(cycle_four), (3, 1, 0))

    def test_inertia_is_preserved_by_exact_congruence(self) -> None:
        randomizer = random.Random(0xD6_06)
        for size in range(1, 9):
            for _ in range(20):
                signs = [randomizer.choice((-1, 0, 1)) for _ in range(size)]
                transform = [
                    [
                        int(row == column)
                        if column >= row
                        else randomizer.randint(-3, 3)
                        for column in range(size)
                    ]
                    for row in range(size)
                ]
                matrix = [
                    [
                        sum(
                            transform[pivot][row]
                            * signs[pivot]
                            * transform[pivot][column]
                            for pivot in range(size)
                        )
                        for column in range(size)
                    ]
                    for row in range(size)
                ]
                expected = (
                    signs.count(1),
                    signs.count(-1),
                    signs.count(0),
                )
                self.assertEqual(exact_inertia_rows(matrix), expected)

    def test_all_labeled_graphs_through_five_against_sympy(self) -> None:
        try:
            import sympy as sp
        except ImportError:
            self.skipTest("SymPy is unavailable for the independent cross-check")
        variable = sp.symbols("x")
        checked = 0
        for size in range(1, 6):
            edges = tuple(combinations(range(size), 2))
            for encoded in range(1 << len(edges)):
                matrix = i_plus_adjacency(
                    size,
                    tuple(
                        edge
                        for bit, edge in enumerate(edges)
                        if encoded & (1 << bit)
                    ),
                )
                polynomial = sp.Matrix(matrix).charpoly(variable).as_poly()
                positive = 0
                negative = 0
                zero = 0
                for (left, right), multiplicity in polynomial.intervals(
                    eps=sp.Rational(1, 10**30)
                ):
                    if right < 0:
                        negative += multiplicity
                    elif left > 0:
                        positive += multiplicity
                    elif left == right == 0:
                        zero += multiplicity
                    else:
                        self.fail(f"root interval straddles zero: {left},{right}")
                self.assertEqual(
                    exact_inertia_rows(matrix), (positive, negative, zero)
                )
                checked += 1
        self.assertEqual(checked, 1_099)

    def test_structured_bound_can_require_seven_dimensions(self) -> None:
        # These are two exact I+Adj inertia triples occurring in a corpus
        # witness.  The older arbitrary-pattern zero-forcing lower bounds of
        # the two sides are only two apiece, while either Lorentz orientation
        # of the structured Gram matrices needs seven dimensions in total.
        self.assertEqual(
            normal_dimension_lower_bound(4, (2, 1, 1), 6, (2, 1, 3)),
            (7, 7, 7),
        )

    def test_realizable_18_points_pass_every_k6_seed(self) -> None:
        decision = evaluate_normal_graph(lower_bound_18_graph(), scan_all=True)
        self.assertTrue(decision.applicable)
        self.assertFalse(decision.rejected)
        self.assertEqual(decision.seeds_checked, 32)
        self.assertEqual(decision.impossible_seeds, 0)

    def test_fixed_residue_witness_exhausts_every_z0(self) -> None:
        with INPUT.open(encoding="utf-8") as stream:
            payload = json.load(stream)
        record = next(item for item in payload["graphs"] if item["index"] == 3_279)
        decision = evaluate_normal_graph(record["adjacency"], scan_all=True)
        self.assertTrue(decision.rejected)
        self.assertEqual(decision.impossible_seeds, 1)
        self.assertIsNotNone(decision.first_witness)
        assert decision.first_witness is not None
        self.assertEqual(decision.first_witness["seed"], [4, 11, 13, 14, 15, 18])
        self.assertEqual(
            decision.first_witness["counts"]["Z0_subsets_considered"], 127
        )
        self.assertEqual(
            decision.first_witness["counts"]["Z0_inertia_passed"], 0
        )

    def test_frozen_sample_report_if_present(self) -> None:
        if not REPORT.exists():
            self.skipTest("run the 64-graph normal-inertia sample first")
        with REPORT.open(encoding="utf-8") as stream:
            report = json.load(stream)
        self.assertEqual(report["schema"], "d6-k6-normal-inertia-sample-v1")
        self.assertEqual(report["status"], "COMPLETE")
        self.assertEqual(report["sample_size"], 64)
        self.assertEqual(report["graphs_rejected"], 8)
        self.assertEqual(
            set(report["rejected_indices"]),
            {
                95_658,
                420_602,
                502_025,
                1_196_109,
                1_334_813,
                1_860_514,
                2_097_618,
                2_343_525,
            },
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
