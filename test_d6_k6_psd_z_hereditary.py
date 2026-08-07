#!/usr/bin/env python3
"""Source-bound controls for hereditary PSD--Z support Hall."""

from __future__ import annotations

import ast
import json
import tempfile
import unittest
from pathlib import Path

import sympy as sp

import d6_k6_psd_z_hereditary as production
import verify_d6_k6_psd_z_hereditary as independent
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent


class HereditaryPsdZTests(unittest.TestCase):
    def test_irreducible_psd_z_proper_principal_controls(self) -> None:
        singular = sp.Matrix(((1, -1, 0), (-1, 2, -1), (0, -1, 1)))
        positive_definite = sp.Matrix(((2, -1, 0), (-1, 2, -1), (0, -1, 2)))
        self.assertEqual(singular.det(), 0)
        self.assertGreater(positive_definite.det(), 0)
        for matrix in (singular, positive_definite):
            for selected in ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2)):
                self.assertGreater(matrix.extract(selected, selected).det(), 0)

        # Bipartite signature switching changes the signs of all P3 edges and
        # preserves every principal determinant.
        signature = sp.diag(1, -1, 1)
        switched = signature * singular * signature
        self.assertEqual(switched.det(), 0)
        self.assertEqual(switched[0, 1], 1)
        self.assertEqual(switched[1, 2], 1)
        for selected in ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2)):
            self.assertEqual(
                switched.extract(selected, selected).det(),
                singular.extract(selected, selected).det(),
            )

        # Connectedness is necessary: a reducible PSD Z-matrix may have a
        # singular proper principal submatrix.
        reducible = sp.diag(0, 1)
        self.assertEqual(reducible.extract((0,), (0,)).det(), 0)

    def test_one_subset_per_span_prevents_false_overlap_sum(self) -> None:
        controls = production.synthetic_controls()
        self.assertTrue(controls["one_span_overlapping_subsets_passed"])
        self.assertEqual(controls["forbidden_naive_overlap_rank"], 4)
        self.assertEqual(controls["forbidden_naive_overlap_capacity"], 3)
        self.assertTrue(controls["two_distinct_spans_failed"])
        self.assertEqual(independent.kernel_controls()["one_span_overlap_passed"], True)

    def test_independent_checker_import_boundary(self) -> None:
        tree = ast.parse(
            (ROOT / "verify_d6_k6_psd_z_hereditary.py").read_text(encoding="utf-8")
        )
        imported = {
            node.module for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        imported.update(
            alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names
        )
        self.assertNotIn("d6_k6_psd_z_hereditary", imported)
        self.assertNotIn("probe_d6_k6_psd_z_hereditary", imported)

    def test_pinned_861_input_boundary(self) -> None:
        _, records, indices = production.verify_inputs()
        self.assertEqual(len(records), 861)
        self.assertEqual(len(indices), 861)
        self.assertEqual(production.stable_hash(indices), production.EXPECTED_INPUT_SHA256)

    def test_fixed_pilot_rejection_agrees_independently(self) -> None:
        _, records, _ = production.verify_inputs()
        record = next(row for row in records if row["index"] == 95_496)
        result = production.evaluate_graph(record["adjacency"])
        fresh = independent.evaluate_graph(record["adjacency"])
        self.assertTrue(result["rejected"])
        self.assertTrue(fresh["rejected"])
        self.assertEqual(result["first_impossible_seed"], [0, 1, 3, 6, 14, 18])
        self.assertEqual(fresh["first_impossible_seed"], result["first_impossible_seed"])
        self.assertEqual(len(result["certificate"]["choices"]), 127)

    def test_known_realizable_positive_control(self) -> None:
        graph = lower_bound_18_graph()
        result = production.evaluate_graph(graph)
        fresh = independent.evaluate_graph(graph)
        self.assertFalse(result["rejected"])
        self.assertFalse(fresh["rejected"])
        self.assertEqual(result["seeds_checked"], 32)
        self.assertEqual(fresh["seeds_checked"], 32)

    def test_checkpoint_prefix_contract(self) -> None:
        payload = production.checkpoint_payload("source", (3, 5), [{"index": 3}])
        self.assertEqual(payload["schema"], production.CHECKPOINT_SCHEMA)
        self.assertEqual(payload["completed"], [{"index": 3}])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            production.atomic_json(path, payload)
            self.assertEqual(json.loads(path.read_text()), payload)


if __name__ == "__main__":
    unittest.main()
