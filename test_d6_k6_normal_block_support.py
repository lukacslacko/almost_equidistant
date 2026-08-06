#!/usr/bin/env python3
"""Exact controls for the K6 normal block-support refinement."""

from __future__ import annotations

import json
import unittest

from d6_k6_bipartite_rank_reference import lorentz_components
from d6_k6_lorentz_reference import build_instance, vertices
from d6_k6_normal_block_support import (
    ROOT,
    InertiaCache,
    OrthogonalBlock,
    check_block_subsets,
    check_component,
    evaluate_graph,
    sha256,
)
from test_d6_k6_lorentz import lower_bound_18_graph


INPUT = ROOT / "d6_k6_bipartite_rank_input.json"
REPORT = ROOT / "d6_k6_normal_block_support_report.json"
SOURCE = ROOT / "d6_k6_normal_block_support.py"


class BlockSubsetKernelControls(unittest.TestCase):
    def test_proper_subset_can_strengthen_the_ambient_bound(self) -> None:
        # The full set has total rank and coordinate-union size six, so the
        # old ambient test passes.  The A block alone needs five dimensions
        # but is allowed in only four coordinates.
        blocks = (
            OrthogonalBlock("A", 5, 0b001111),
            OrthogonalBlock("B", 1, 0b110000),
        )
        self.assertEqual(sum(block.rank_lower for block in blocks), 6)
        self.assertEqual(
            (blocks[0].allowed_coordinates | blocks[1].allowed_coordinates)
            .bit_count(),
            6,
        )
        result = check_block_subsets(blocks)
        self.assertFalse(result.passed)
        self.assertEqual(result.first_failure["blocks"], ["A"])
        self.assertEqual(result.first_failure["rank_lower"], 5)
        self.assertEqual(result.first_failure["coordinate_capacity"], 4)

    def test_allowed_mask_enlargement_cannot_create_a_failure(self) -> None:
        blocks = (
            OrthogonalBlock("A", 2, 0b000011),
            OrthogonalBlock("B", 2, 0b001100),
            OrthogonalBlock("Z0:7", 1, 0b010000),
        )
        self.assertTrue(check_block_subsets(blocks).passed)
        enlarged = tuple(
            OrthogonalBlock(
                block.label,
                block.rank_lower,
                block.allowed_coordinates | 0b100000,
            )
            for block in blocks
        )
        self.assertTrue(check_block_subsets(enlarged).passed)

        # An upper mask may contain unused coordinates.  Rank zero in an
        # empty allowed mask is legal; nonedges are never forced non-unit.
        self.assertTrue(
            check_block_subsets((OrthogonalBlock("optional-zero", 0, 0),)).passed
        )

    def test_at_most_eight_blocks_are_admitted(self) -> None:
        blocks = tuple(OrthogonalBlock(str(i), 0, 0) for i in range(9))
        with self.assertRaises(ValueError):
            check_block_subsets(blocks)


class CorpusAndSemanticControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        payload = json.loads(INPUT.read_text(encoding="utf-8"))
        cls.by_index = {
            record["index"]: record["adjacency"] for record in payload["graphs"]
        }

    def test_fixed_orientation_level_strengthening(self) -> None:
        adjacency = self.by_index[2_552_186]
        instance = build_instance(adjacency, (2, 5, 7, 10, 13, 17))
        wanted = {0, 4, 8, 9, 12, 14, 15, 16, 18}
        component = next(
            item
            for item in lorentz_components(instance, 0)
            if {
                instance.outside[local] for local in vertices(item.component)
            }
            == wanted
        )
        self.assertEqual(
            [instance.outside[local] for local in vertices(component.side_a)],
            [0, 9, 12, 14, 16, 18],
        )
        self.assertEqual(
            [instance.outside[local] for local in vertices(component.side_b)],
            [4, 8, 15],
        )

        decision = check_component(
            adjacency, instance, 0, component, InertiaCache()
        )
        positive, negative = decision.generic_cases
        self.assertTrue(positive["ambient_dimension_passed"])
        self.assertEqual(positive["ambient_dimension_lower"], 6)
        self.assertFalse(positive["passed"])
        self.assertEqual(positive["first_failure"]["blocks"], ["A"])
        self.assertEqual(positive["first_failure"]["rank_lower"], 5)
        self.assertEqual(positive["first_failure"]["coordinate_capacity"], 4)
        self.assertEqual(
            positive["first_failure"]["allowed_coordinates"], [2, 3, 4, 5]
        )
        self.assertTrue(negative["passed"])
        self.assertFalse(decision.lightlike_case["passed"])
        self.assertTrue(decision.passed)

    def test_known_realizable_eighteen_point_graph_passes(self) -> None:
        decision = evaluate_graph(lower_bound_18_graph())
        self.assertFalse(decision["rejected"])
        self.assertEqual(decision["seeds_checked"], 32)
        self.assertEqual(decision["impossible_seeds"], 0)
        self.assertEqual(decision["generic_cases_support_new_failures"], 0)

    def test_alpha_two_hypothesis_is_checked(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_graph([0, 0, 0])

    def test_frozen_full_report(self) -> None:
        report = json.loads(REPORT.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "COMPLETE")
        self.assertEqual(report["input_graphs"], 990)
        self.assertEqual(report["graphs_rejected"], 0)
        self.assertEqual(report["graphs_surviving"], 990)
        self.assertEqual(report["seeds_checked"], 31_654)
        self.assertEqual(report["z0_considered"], 31_673)
        self.assertEqual(report["z0_matchable"], 31_673)
        self.assertEqual(report["z0_failed"], 19)
        self.assertEqual(report["bipartite_components_checked"], 218_377)
        self.assertEqual(report["generic_cases_checked"], 436_754)
        self.assertEqual(report["generic_cases_support_new_failures"], 2)
        self.assertEqual(report["block_subsets_checked"], 1_310_900)
        self.assertEqual(report["components_old_passed_new_failed"], 0)
        self.assertEqual(report["z0_old_passed_new_failed"], 0)
        strengthened = [
            (item["index"], item["decision"]["generic_cases_support_new_failures"])
            for item in report["graph_results"]
            if item["decision"]["generic_cases_support_new_failures"]
        ]
        self.assertEqual(strengthened, [(2_552_186, 1), (2_870_953, 1)])
        self.assertEqual(
            report["sources"][SOURCE.name], sha256(SOURCE)
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
