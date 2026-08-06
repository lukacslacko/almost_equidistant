#!/usr/bin/env python3
"""Controls for the required-edge K7 Mobius-cycle obstructions."""

from __future__ import annotations

import ast
import unittest

import d6_k7_mobius_cycle_containment as production
import verify_d6_k7_mobius_cycle_containment as independent


class K7MobiusCycleControls(unittest.TestCase):
    def test_exact_mobius_order_and_discriminants(self) -> None:
        result = independent.verify_mobius_algebra()
        self.assertEqual(
            result["fixed_point_discriminants_T1_through_T5"],
            [-3, -27, -108, -243, -243],
        )
        self.assertEqual(result["retained_negative_control_cycle_length"], 6)

    def test_canonical_required_edge_patterns(self) -> None:
        triangle = production.PATTERNS[3]
        quadrilateral = production.PATTERNS[4]
        self.assertEqual((triangle["order"], triangle["edges"]), (10, 39))
        self.assertEqual(
            (quadrilateral["order"], quadrilateral["edges"]), (11, 45)
        )
        self.assertFalse(triangle["nonedges_used_as_distance_constraints"])
        self.assertFalse(
            quadrilateral["nonedges_used_as_distance_constraints"]
        )
        self.assertEqual(triangle, independent.expected_pattern(3))
        self.assertEqual(quadrilateral, independent.expected_pattern(4))

    def test_positive_and_negative_containment_controls(self) -> None:
        for length in (3, 4):
            pattern = production.PATTERNS[length]["adjacency"]
            embedded = production.embedded_target(pattern)
            mapping = production.find_embedding(embedded, length)
            self.assertIsNotNone(mapping)
            production.verify_mapping(pattern, embedded, mapping)
            other, _, _ = independent.independent_embedding(embedded, length)
            self.assertIsNotNone(other)
            independent.verify_mapping(pattern, embedded, other)
            self.assertIsNone(
                production.find_embedding((0,) * production.TARGET_ORDER, length)
            )
            self.assertIsNone(
                independent.independent_embedding(
                    (0,) * production.TARGET_ORDER, length
                )[0]
            )

    def test_extra_unit_edges_do_not_hurt_embedding(self) -> None:
        complete = production.complete_target()
        for length in (3, 4):
            mapping = production.find_embedding(complete, length)
            self.assertIsNotNone(mapping)
            production.verify_mapping(
                production.PATTERNS[length]["adjacency"], complete, mapping
            )

    def test_tampered_mapping_is_rejected(self) -> None:
        pattern = production.PATTERNS[3]["adjacency"]
        embedded = list(production.embedded_target(pattern))
        embedded[0] &= ~(1 << 1)
        embedded[1] &= ~(1 << 0)
        with self.assertRaisesRegex(ValueError, "loses required edge"):
            production.verify_mapping(pattern, embedded, list(range(10)))

    def test_independent_checker_does_not_import_production(self) -> None:
        syntax = ast.parse(
            independent.Path(independent.__file__).read_text(encoding="utf-8")
        )
        imports = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("d6_k7_mobius_cycle_containment", imports)


if __name__ == "__main__":
    unittest.main(verbosity=2)
