#!/usr/bin/env python3
"""Controls for the cross-method dimension-six v6 residue boundary."""

from __future__ import annotations

import ast
import copy
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_current_residue_manifest_v6 as builder
import verify_d6_current_residue_manifest_v6 as checker


class CurrentResidueV6Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = builder.build_manifest()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.path = Path(cls.temporary.name) / "v6.json"
        cls.write(cls.manifest, cls.path)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    @staticmethod
    def write(payload: object, path: Path) -> None:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def test_all_primary_and_transitive_hashes_are_bound(self) -> None:
        for name, expected in builder.PRIMARY_ARTIFACTS.items():
            self.assertEqual(builder.sha256(builder.ROOT / name), expected, name)
        self.assertEqual(
            checker.sha256(checker.ROOT / checker.BUILDER),
            checker.EXPECTED_BUILDER_SHA256,
        )
        for name, expected in self.manifest["source_boundary"][
            "transitive_files"
        ].items():
            self.assertEqual(builder.sha256(builder.ROOT / name), expected, name)

    def test_exact_counts_hashes_and_trust_tiers(self) -> None:
        k7 = self.manifest["classes"]["K7"]
        k6 = self.manifest["classes"]["K6_only"]
        self.assertEqual(k7["count"], 19)
        self.assertEqual(k7["indices"], builder.EXPECTED_K7_UNION_SURVIVORS)
        self.assertEqual(
            k7["indices_sha256"], builder.EXPECTED_K7_UNION_SURVIVORS_SHA256
        )
        self.assertEqual(k7["exact_algebra_rejection_count"], 133)
        self.assertEqual(k7["interval_new_rejected_indices"], [423661, 424226, 3936176])
        self.assertEqual(k6["count"], 625)
        self.assertEqual(k6["indices_sha256"], builder.EXPECTED_K6_SURVIVORS_SHA256)
        self.assertEqual(k6["exact_rejection_count"], 131)
        self.assertEqual(self.manifest["combined"]["count"], 644)
        accounting = self.manifest["certificate_accounting"]
        self.assertFalse(
            accounting["exact_algebra_and_graph_logic"][
                "floating_point_enters_rejection"
            ]
        )
        self.assertEqual(accounting["interval"]["incremental_rejections"], 3)
        self.assertFalse(
            accounting["interval"]["non_KILLED_statuses_used_for_rejection"]
        )
        self.assertIn(
            "8 ulps", accounting["interval"]["trust_assumptions"]["transcendentals"]
        )

    def test_independent_checker_reconstructs_every_record(self) -> None:
        result = checker.verify(self.path, None)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(
            result["counts"],
            {
                "K7": 19,
                "K6_only": 625,
                "combined": 644,
                "exact_rejections_from_v5": 264,
                "interval_incremental_rejections": 3,
            },
        )
        self.assertTrue(all(result["checks"].values()))

    def test_graph_or_trust_tier_tamper_is_rejected(self) -> None:
        graph_tamper = copy.deepcopy(self.manifest)
        graph_tamper["classes"]["K6_only"]["graphs"][0]["adjacency"][0] ^= 2
        graph_path = Path(self.temporary.name) / "graph-tamper.json"
        self.write(graph_tamper, graph_path)
        with self.assertRaises(ValueError):
            checker.verify(graph_path, None)

        trust_tamper = copy.deepcopy(self.manifest)
        trust_tamper["certificate_accounting"]["interval"][
            "non_KILLED_statuses_used_for_rejection"
        ] = True
        trust_path = Path(self.temporary.name) / "trust-tamper.json"
        self.write(trust_tamper, trust_path)
        with self.assertRaises(ValueError):
            checker.verify(trust_path, None)

    def test_checker_imports_no_builder_or_production_kernel(self) -> None:
        tree = ast.parse(
            (checker.ROOT / "verify_d6_current_residue_manifest_v6.py").read_text(
                encoding="utf-8"
            )
        )
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest_v6", imports)
        self.assertFalse(
            any(
                name.startswith(
                    ("d6_k6_", "d6_k7_", "verify_d6_k6_", "verify_d6_k7_")
                )
                for name in imports
            )
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
