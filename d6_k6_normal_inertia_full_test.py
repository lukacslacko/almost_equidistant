#!/usr/bin/env python3
"""Controls for the full K6 normal-inertia runner and independent verifier."""

from __future__ import annotations

import ast
import json
import unittest
from pathlib import Path

from d6_k6_normal_inertia import evaluate_normal_graph
from d6_k6_normal_inertia_full_runner import (
    ROOT,
    select_residue,
    verify_dependencies,
)
from d6_k6_normal_inertia_full_verifier import (
    SCALAR_FIELDS,
    SympyInertiaCache,
    evaluate_graph,
)
from test_d6_k6_lorentz import lower_bound_18_graph


INPUT = ROOT / "d6_k6_bipartite_rank_input.json"
PRIOR = ROOT / "d6_k6_bipartite_rank_report.json"
VERIFIER = ROOT / "d6_k6_normal_inertia_full_verifier.py"


class K6NormalFullControls(unittest.TestCase):
    def test_runner_hash_pins_and_selects_exact_residue(self) -> None:
        hashes = verify_dependencies(INPUT, PRIOR)
        self.assertEqual(len(hashes), 9)
        selected = select_residue(INPUT, PRIOR)
        self.assertEqual(len(selected), 1_097)
        self.assertNotIn(461_363, {record["index"] for record in selected})

    def test_verifier_does_not_import_production_modules(self) -> None:
        syntax = ast.parse(VERIFIER.read_text(encoding="utf-8"))
        imported = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
        self.assertNotIn("d6_k6_normal_inertia", imported)
        self.assertNotIn("d6_k6_normal_inertia_full_runner", imported)

    def test_sympy_inertia_has_independent_exact_controls(self) -> None:
        controls = {
            (): (0, 0, 0),
            (1,): (1, 0, 0),
            (0b11, 0b11): (1, 0, 1),
            # I+Adj(P3)
            (0b011, 0b111, 0b110): (2, 1, 0),
        }
        for rows, expected in controls.items():
            with self.subTest(rows=rows):
                self.assertEqual(SympyInertiaCache().solve(rows), expected)

    def test_independent_decisions_match_fixed_corpus_controls(self) -> None:
        payload = json.loads(INPUT.read_text(encoding="utf-8"))
        by_index = {record["index"]: record for record in payload["graphs"]}
        indices = (
            3_278,
            3_279,
            95_658,
            420_602,
            502_025,
            1_196_109,
            1_334_813,
            1_860_514,
            2_097_618,
            2_343_525,
        )
        for index in indices:
            with self.subTest(index=index):
                adjacency = by_index[index]["adjacency"]
                production = evaluate_normal_graph(adjacency, scan_all=True)
                independent = evaluate_graph(adjacency)
                for name in SCALAR_FIELDS:
                    self.assertEqual(getattr(independent, name), getattr(production, name))
                witness = production.first_witness
                production_seed = tuple(witness["seed"]) if witness else None
                self.assertEqual(independent.first_impossible_seed, production_seed)

    def test_positive_control_matches_both_engines(self) -> None:
        adjacency = lower_bound_18_graph()
        production = evaluate_normal_graph(adjacency, scan_all=True)
        independent = evaluate_graph(adjacency)
        self.assertFalse(production.rejected)
        self.assertFalse(independent.rejected)
        self.assertEqual(production.seeds_checked, 32)
        self.assertEqual(independent.seeds_checked, 32)


if __name__ == "__main__":
    unittest.main(verbosity=2)
