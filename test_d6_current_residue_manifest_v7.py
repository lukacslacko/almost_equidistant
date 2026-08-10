#!/usr/bin/env python3
"""Controls for the source-bound dimension-six v7 residue manifest."""

from __future__ import annotations

import ast
import copy
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_current_residue_manifest_v7 as builder
import verify_d6_current_residue_manifest_v7 as checker


class CurrentResidueV7Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        package = {
            "commit": "0" * 40,
            "branch": "TEST",
            "package_sources": {
                name: builder.sha256(builder.ROOT / name)
                for name in builder.PACKAGE_FILES
            },
            "porcelain_lines": [],
            "porcelain_sha256": builder.stable_hash("TEST"),
            "package_sources_equal_committed_blobs": False,
        }
        cls.manifest = builder.build_manifest(package)
        cls.temporary = tempfile.TemporaryDirectory()
        cls.path = Path(cls.temporary.name) / "v7.json"
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

    def test_final_counts_hashes_and_disjoint_layers(self) -> None:
        k7 = self.manifest["classes"]["K7"]
        k6 = self.manifest["classes"]["K6_only"]
        self.assertEqual(k7["count"], 12)
        self.assertEqual(k7["indices"], builder.EXPECTED_K7)
        self.assertEqual(k7["indices_sha256"], builder.EXPECTED_K7_SHA256)
        self.assertEqual(k7["post_v6_exact_rejection_count"], 5)
        self.assertEqual(k7["post_v6_interval_rejection_count"], 2)
        self.assertEqual(k7["post_v6_union_rejection_count"], 7)
        self.assertFalse(
            set(k7["post_v6_exact_rejected_indices"])
            & set(k7["post_v6_interval_rejected_indices"])
        )
        self.assertEqual(k6["count"], 251)
        self.assertEqual(k6["post_v6_exact_rejection_count"], 374)
        self.assertEqual(k6["indices_sha256"], builder.EXPECTED_K6_SHA256)
        self.assertEqual(self.manifest["combined"]["count"], 263)

    def test_all_primary_and_layer_source_roots_are_present(self) -> None:
        boundary = self.manifest["source_boundary"]
        self.assertEqual(
            boundary["primary_artifacts"], dict(sorted(builder.PRIMARY_ARTIFACTS.items()))
        )
        self.assertEqual(
            set(boundary["layer_source_roots"]),
            {
                "inherited_v6",
                "interval_cap500000",
                "exact_star",
                "interval_focused_3936177",
                "exact_schur_3949382",
                "exact_2593240",
                "exact_K6_chain",
            },
        )
        self.assertTrue(boundary["layer_source_roots"]["exact_2593240"])

    def test_trust_tiers_remain_separate(self) -> None:
        accounting = self.manifest["certificate_accounting"]
        self.assertEqual(accounting["exact"]["total_rejections_from_v5"], 643)
        self.assertFalse(accounting["exact"]["floating_point_enters_rejection"])
        self.assertEqual(accounting["interval"]["total_incremental_rejections_from_v5"], 5)
        self.assertTrue(accounting["interval"]["only_KILLED_used_for_rejection"])
        self.assertFalse(
            accounting["interval"]["ABORT_UNRESOLVED_INFRA_ERROR_used_for_rejection"]
        )
        self.assertIn(
            "8 ulps", accounting["interval"]["trust_assumptions"]["transcendentals"]
        )
        self.assertEqual(accounting["union_rejections_from_v5"], 648)

    def test_import_independent_checker_reconstructs_manifest(self) -> None:
        result = checker.verify(
            self.path,
            checker.sha256(self.path),
            enforce_source_boundary=False,
        )
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(
            result["counts"],
            {
                "K7": 12,
                "K6_only": 251,
                "combined": 263,
                "exact_rejections_from_v5": 643,
                "interval_incremental_rejections": 5,
            },
        )
        self.assertTrue(all(result["checks"].values()))

    def test_graph_trust_and_source_root_tampering_are_rejected(self) -> None:
        mutations = []
        graph = copy.deepcopy(self.manifest)
        graph["classes"]["K6_only"]["graphs"][0]["adjacency"][0] ^= 2
        mutations.append(graph)
        trust = copy.deepcopy(self.manifest)
        trust["certificate_accounting"]["interval"][
            "ABORT_UNRESOLVED_INFRA_ERROR_used_for_rejection"
        ] = True
        mutations.append(trust)
        source = copy.deepcopy(self.manifest)
        source["source_boundary"]["layer_source_roots"]["exact_star"][
            "launch_commit"
        ] = "0" * 40
        mutations.append(source)
        for number, payload in enumerate(mutations):
            path = Path(self.temporary.name) / f"tamper-{number}.json"
            self.write(payload, path)
            with self.assertRaises(ValueError):
                checker.verify(path, None, enforce_source_boundary=False)

    def test_checker_imports_no_builder_or_production_module(self) -> None:
        source = (checker.ROOT / "verify_d6_current_residue_manifest_v7.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest_v7", imports)
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
