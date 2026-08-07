#!/usr/bin/env python3
"""Source-bound controls for exact K6 arbitrary-subset rank Hall."""

from __future__ import annotations

import ast
import json
import tempfile
import unittest
from pathlib import Path

import d6_k6_arbitrary_subset_hall as production
import verify_d6_k6_arbitrary_subset_hall as independent
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent


class ArbitrarySubsetHallTests(unittest.TestCase):
    def test_production_and_independent_synthetic_controls(self) -> None:
        produced = production.synthetic_controls()
        checked = independent.kernel_controls()
        self.assertEqual(produced["C5_selected_P3_rank"], 2)
        self.assertEqual(checked["C5_selected_P3_rank"], 2)
        self.assertTrue(produced["arbitrary_subset_failed"])
        self.assertTrue(checked["arbitrary_C5_Hall_failed"])
        self.assertTrue(produced["one_vertex_extension_omitted"])
        self.assertTrue(checked["one_vertex_extension_omitted"])

    def test_independent_container_hall_equals_raw_cartesian(self) -> None:
        spans = (
            independent.prior.Span(
                "A",
                (
                    independent.prior.Option(0, 0, (), "empty"),
                    independent.prior.Option(1, 1, (0,), "proper"),
                    independent.prior.Option(2, 3, (0, 1), "proper"),
                ),
                False,
            ),
            independent.prior.Span(
                "B",
                (
                    independent.prior.Option(0, 0, (), "empty"),
                    independent.prior.Option(1, 2, (2,), "proper"),
                    independent.prior.Option(2, 6, (2, 3), "proper"),
                ),
                False,
            ),
        )
        raw_passed, _ = independent.prior.exact_hall(spans)
        container_passed, _ = independent.independent_hall(spans)
        self.assertEqual(container_passed, raw_passed)
        self.assertFalse(container_passed)

    def test_one_subset_per_original_span(self) -> None:
        overlap = production.parent.component_group(
            "one-span", (0, 1, 2), (1, 2, 4), 2, True
        )
        self.assertTrue(production.parent.check_groups((overlap,)).passed)
        other = production.parent.component_group(
            "other-span", (3, 4, 5), (1, 2, 4), 2, True
        )
        self.assertFalse(production.parent.check_groups((overlap, other)).passed)

    def test_pinned_831_boundary(self) -> None:
        _, records, indices = production.verify_inputs()
        self.assertEqual(len(records), 831)
        self.assertEqual(len(indices), 831)
        self.assertEqual(production.stable_hash(indices), production.EXPECTED_INPUT_SHA256)

    def test_fixed_pilot_rejection_agrees_independently(self) -> None:
        _, records, _ = production.verify_inputs()
        record = next(row for row in records if row["index"] == 650_158)
        produced = production.evaluate_graph(record["adjacency"])
        checked = independent.evaluate_graph(record["adjacency"])
        seed = [0, 2, 5, 6, 14, 16]
        self.assertTrue(produced["rejected"])
        self.assertTrue(checked["rejected"])
        self.assertEqual(produced["first_impossible_seed"], seed)
        self.assertEqual(checked["first_impossible_seed"], seed)
        self.assertEqual(len(produced["certificate"]["choices"]), 32)

    def test_known_realizable_positive_control(self) -> None:
        result = production.evaluate_graph(lower_bound_18_graph())
        self.assertFalse(result["rejected"])
        self.assertEqual(result["seeds_checked"], 32)

    def test_import_boundary_checkpoint_and_no_extension_kernel(self) -> None:
        tree = ast.parse(
            (ROOT / "verify_d6_k6_arbitrary_subset_hall.py").read_text()
        )
        imports = {
            node.module for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        imports.update(
            alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names
        )
        self.assertNotIn("d6_k6_arbitrary_subset_hall", imports)
        self.assertNotIn("probe_d6_k6_arbitrary_subset_hall", imports)
        production_tree = ast.parse(
            (ROOT / "d6_k6_arbitrary_subset_hall.py").read_text()
        )
        function_names = {
            node.name for node in ast.walk(production_tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        self.assertNotIn("has_one_vertex_bipartite_extension", function_names)
        payload = production.checkpoint_payload("source", (3, 5), [{"index": 3}])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            production.atomic_json(path, payload)
            self.assertEqual(json.loads(path.read_text()), payload)


if __name__ == "__main__":
    unittest.main()
