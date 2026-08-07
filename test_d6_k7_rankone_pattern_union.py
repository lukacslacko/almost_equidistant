#!/usr/bin/env python3
"""Controls for the exact K7 tetrad/pattern union accounting."""

from __future__ import annotations

import ast
import json
import unittest

import build_d6_k7_rankone_pattern_union as union
import verify_d6_k7_rankone_pattern_union as checker


MANIFEST = union.OUTPUT
MANIFEST_SHA256 = (
    "1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599"
)


def row(**updates: int) -> dict:
    value = {
        "index": 7,
        "prior_dual_passing_covers": 0,
        "prior_passing_no_saturating_clique_covers": 0,
        "covers_with_near_clique": 0,
        "covers_without_near_clique": 0,
        "tetrad_failed_covers": 0,
        "tetrad_passing_covers": 0,
    }
    value.update(updates)
    return value


class RankonePatternUnionControls(unittest.TestCase):
    def test_fresh_rebuild_is_identical(self) -> None:
        frozen = json.loads(MANIFEST.read_text(encoding="utf-8"))
        gates = frozen["gates"]["tetrad"]
        self.assertEqual(
            union.build_manifest(
                gates["report_sha256"], gates["verification_sha256"]
            ),
            frozen,
        )
        self.assertEqual(union.sha256(MANIFEST), MANIFEST_SHA256)

    def test_independent_full_checker_passes(self) -> None:
        result = checker.verify(MANIFEST, MANIFEST_SHA256)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["graphs"], 12_839)
        self.assertEqual(result["tetrad_rejected"], 11_902)
        self.assertEqual(result["pattern_954_rejected"], 3_403)
        self.assertEqual(result["intersection"], 2_724)
        self.assertEqual(result["exact_union"], 12_581)
        self.assertEqual(result["exact_residue"], 258)
        self.assertEqual(
            result["residue_indices_sha256"],
            "55ab329dbe3ae4dbfa10378ff023a168fc6ab306ffc7af14bbd9790a68108a09",
        )

    def test_exact_set_partition_and_residue_profiles(self) -> None:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        sets = manifest["sets"]
        tetrad = set(sets["tetrad_rejected"]["indices"])
        pattern = set(sets["pattern_954_rejected"]["indices"])
        universe = set(
            json.loads(union.SELECTION.read_text(encoding="utf-8"))[
                "selected_indices"
            ]
        )
        self.assertEqual(tetrad & pattern, set(sets["both"]["indices"]))
        self.assertEqual(tetrad | pattern, set(sets["exact_union"]["indices"]))
        self.assertEqual(universe - tetrad - pattern,
                         set(sets["exact_residue"]["indices"]))
        profiles = manifest["cover_structure"]["residue_profiles"]
        self.assertEqual(len(profiles), 258)
        self.assertEqual(
            [row["index"] for row in profiles],
            sets["exact_residue"]["indices"],
        )

    def test_cover_profile_three_structures(self) -> None:
        saturating = union.cover_profile(row(
            prior_dual_passing_covers=2,
            tetrad_passing_covers=2,
        ))
        self.assertEqual(saturating["structure"], "saturating_only")
        self.assertTrue(saturating["has_saturating_cover"])

        mixed = union.cover_profile(row(
            prior_dual_passing_covers=5,
            prior_passing_no_saturating_clique_covers=3,
            covers_with_near_clique=2,
            covers_without_near_clique=1,
            tetrad_failed_covers=1,
            tetrad_passing_covers=4,
        ))
        self.assertEqual(mixed["structure"], "mixed")
        self.assertEqual(mixed["saturating_covers"], 2)
        self.assertEqual(mixed["tetrad_resistant_near_clique_covers"], 1)

        no_saturating = union.cover_profile(row(
            prior_dual_passing_covers=2,
            prior_passing_no_saturating_clique_covers=2,
            covers_with_near_clique=2,
            tetrad_failed_covers=2,
        ))
        self.assertEqual(no_saturating["structure"], "no_saturating_only")
        self.assertEqual(no_saturating["tetrad_passing_covers"], 0)

    def test_cover_profile_rejects_bad_partition(self) -> None:
        with self.assertRaisesRegex(ValueError, "partition"):
            union.cover_profile(row(
                prior_dual_passing_covers=2,
                prior_passing_no_saturating_clique_covers=2,
                covers_with_near_clique=1,
            ))

    def test_frozen_pattern_gate_passes(self) -> None:
        selection = json.loads(
            union.SELECTION.read_text(encoding="utf-8")
        )["selected_indices"]
        hits = union.validate_pattern_gate(
            json.loads(union.PATTERN_REPORT.read_text(encoding="utf-8")),
            json.loads(
                union.PATTERN_VERIFICATION.read_text(encoding="utf-8")
            ),
            selection,
        )
        self.assertEqual(len(selection), 12_839)
        self.assertEqual(len(hits), 3_403)

    def test_independent_checker_does_not_import_builder(self) -> None:
        syntax = ast.parse(union.VERIFIER.read_text(encoding="utf-8"))
        imports = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_k7_rankone_pattern_union", imports)


if __name__ == "__main__":
    unittest.main(verbosity=2)
