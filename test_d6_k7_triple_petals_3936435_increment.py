#!/usr/bin/env python3
"""Focused controls for the graph-3936435 triple-petal package."""

from __future__ import annotations

import ast
import copy
import hashlib
import json
import tempfile
import unittest
from fractions import Fraction
from itertools import combinations

import build_d6_k7_triple_petals_3936435_increment as builder
import verify_d6_k7_triple_petals_3936435_increment as checker


class TriplePetalIncrementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest, cls.parent, cls.parent_check = builder.load_upstream()
        cls.builder_certificate = builder.build_certificate(cls.manifest, cls.parent)
        independent_manifest, independent_parent, _ = checker.load_upstream()
        cls.checker_rebuild = checker.independently_rebuild(
            independent_manifest, independent_parent
        )

    def test_exact_q_sqrt7_arithmetic(self) -> None:
        radical = builder.SQRT7
        self.assertEqual(radical * radical, builder.Q7.of(7))
        for epsilon in (-1, 1):
            leaf = (radical - epsilon) / 2
            self.assertEqual(
                (builder.Q7.of(epsilon) - radical) * (leaf - radical),
                builder.Q7.of(3),
            )
            self.assertNotEqual(leaf * leaf, builder.Q7.of(3))

    def test_full_hub_implication_in_producer(self) -> None:
        algebra = self.builder_certificate["algebra"]
        one = builder.Q7.of(1).json()
        for row in algebra["cases"]:
            leaf = builder.Q7(
                Fraction(*row["leaf_value"][0]),
                Fraction(*row["leaf_value"][1]),
            )
            inverse = builder.Q7(
                Fraction(*row["leaf_inverse"][0]),
                Fraction(*row["leaf_inverse"][1]),
            )
            self.assertNotEqual(leaf, builder.Q7())
            self.assertEqual(leaf * inverse, builder.Q7.of(1))
            self.assertEqual(row["six_hub_petal_products"], [one] * 6)
            self.assertEqual(
                row["leaf_square_times_hub_dot_minus_one"],
                row["three_minus_leaf_square"],
            )

    def test_two_cover_quantifier(self) -> None:
        quantifier = self.builder_certificate["quantifier"]
        self.assertEqual(quantifier["raw_eligible_covers"], 255)
        self.assertEqual(quantifier["prior_layer_eliminated_raw_covers"], 253)
        self.assertEqual(quantifier["upstream_current_covers"], [0, 2048])
        self.assertEqual(quantifier["new_families_checked"], 2)
        self.assertEqual(quantifier["new_families_rejected"], 2)
        self.assertEqual(
            [row["zmask"] for row in self.builder_certificate["covers"]],
            [0, 2048],
        )

    def test_independent_cover_replay_agrees(self) -> None:
        self.assertEqual(len(self.checker_rebuild["covers"]), 255)
        self.assertEqual(
            checker.stable_hash(list(self.checker_rebuild["covers"])),
            checker.EXPECTED_RAW_COVERS_SHA256,
        )
        self.assertEqual(
            [row["zmask"] for row in self.checker_rebuild["current_cover_rows"]],
            [0, 2048],
        )
        self.assertTrue(
            all(
                row["five_masks"] == [60, 10, 6, 45, 34]
                for row in self.checker_rebuild["current_cover_rows"]
            )
        )

    def test_independent_algebra(self) -> None:
        algebra = checker.independent_algebra()
        self.assertEqual(
            algebra["status"], "EXACT_CONTRADICTION_FOR_BOTH_SIGNS"
        )
        self.assertEqual(
            [row["conjugate_product_mod_r_squared_minus_7"] for row in algebra["petal_diagonal_cases"]],
            ["3", "3"],
        )
        self.assertEqual(algebra["hub_ideal_identity"], "0")
        self.assertEqual(algebra["hub_equations_checked"], 7)
        self.assertEqual(
            [row["lambda_conjugate_product"] for row in algebra["petal_diagonal_cases"]],
            ["-3/2", "-3/2"],
        )
        self.assertTrue(
            all(
                row["six_hub_petal_substitutions"] == ["0"] * 6
                for row in algebra["petal_diagonal_cases"]
            )
        )

    def test_required_edge_not_optional_nonedge(self) -> None:
        adjacency = [0] * 10
        for first, second in combinations(builder.FIVE_LOCALS, 2):
            adjacency[first] |= 1 << second
            adjacency[second] |= 1 << first
        masks = {
            0: 60,
            1: 10,
            3: 6,
            6: 45,
            9: 34,
        }
        self.assertTrue(
            builder.triple_petal_pattern_holds(adjacency, tuple(range(10)), masks)
        )
        first, second = builder.HUB_LOCALS
        adjacency[first] &= ~(1 << second)
        adjacency[second] &= ~(1 << first)
        self.assertFalse(
            builder.triple_petal_pattern_holds(adjacency, tuple(range(10)), masks)
        )

    def test_known_realizable_18_control(self) -> None:
        self.assertEqual(
            builder.positive_and_synthetic_controls()["known_realizable_18"],
            {"passed": True, "K7_seeds": 0},
        )
        self.assertEqual(
            checker.independent_positive_controls()["known_realizable_18"],
            {"passed": True, "K7_seeds": 0},
        )

    def test_checker_import_independence(self) -> None:
        checker.assert_import_independence()
        with open(checker.__file__, encoding="utf-8") as stream:
            tree = ast.parse(stream.read())
        imports = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        imports.update(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        self.assertNotIn(builder.__name__, imports)

    def test_source_only_report_boundary(self) -> None:
        sources = {
            name: builder.sha256(builder.ROOT / name)
            for name in builder.PACKAGE_SOURCES
        }
        fake_boundary = {
            "commit": "0" * 40,
            "branch": "codex/dimension6",
            "source_sha256": sources,
            "porcelain_lines": [],
            "porcelain_sha256": hashlib.sha256(b"").hexdigest(),
            "proof_and_checker_sources_equal_committed_blobs": True,
            "tracked_clean": True,
        }
        report = builder.build_report(source_boundary=fake_boundary)
        checker.validate_report_schema(report)
        checker.validate_certificate(
            report, self.checker_rebuild, checker.independent_algebra()
        )
        self.assertEqual(report["summary"]["rejected_indices"], [3936435])
        self.assertEqual(report["summary"]["graphs_surviving"], 11)

    def test_checker_rejects_tampered_archived_algebra(self) -> None:
        sources = {
            name: builder.sha256(builder.ROOT / name)
            for name in builder.PACKAGE_SOURCES
        }
        fake_boundary = {
            "commit": "0" * 40,
            "branch": "codex/dimension6",
            "source_sha256": sources,
            "porcelain_lines": [],
            "porcelain_sha256": hashlib.sha256(b"").hexdigest(),
            "proof_and_checker_sources_equal_committed_blobs": True,
            "tracked_clean": True,
        }
        report = builder.build_report(source_boundary=fake_boundary)
        tampered = copy.deepcopy(report)
        tampered["certificate"]["algebra"]["cases"][0]["leaf_inverse"][0][0] += 1
        tampered["certificate_sha256"] = checker.stable_hash(tampered["certificate"])
        checker.validate_report_schema(tampered)
        with self.assertRaisesRegex(ValueError, "archived Q"):
            checker.validate_certificate(
                tampered, self.checker_rebuild, checker.independent_algebra()
            )

    def test_fake_provenance_end_to_end_verification(self) -> None:
        sources = {
            name: builder.sha256(builder.ROOT / name)
            for name in builder.PACKAGE_SOURCES
        }
        fake_boundary = {
            "commit": "0" * 40,
            "branch": "codex/dimension6",
            "source_sha256": sources,
            "porcelain_lines": [],
            "porcelain_sha256": hashlib.sha256(b"").hexdigest(),
            "proof_and_checker_sources_equal_committed_blobs": True,
            "tracked_clean": True,
        }
        report = builder.build_report(source_boundary=fake_boundary)
        with tempfile.TemporaryDirectory() as directory:
            path = builder.Path(directory) / "fake-report.json"
            path.write_text(
                json.dumps(report, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            verification = checker.verify_report(
                path,
                builder.sha256(path),
                enforce_source_boundary=False,
            )
        self.assertEqual(verification["status"], "PASS")
        self.assertEqual(
            verification["independent_algebra"]["hub_ideal_identity"], "0"
        )

    def test_fraction_serialization_is_canonical(self) -> None:
        value = builder.Q7(Fraction(-3, 4), Fraction(5, 6))
        self.assertEqual(value.json(), [[-3, 4], [5, 6]])


if __name__ == "__main__":
    unittest.main()
