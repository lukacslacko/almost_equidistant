#!/usr/bin/env python3
"""Controls for the exact dimension-six current-residue manifest."""

from __future__ import annotations

import ast
import json
import unittest
from pathlib import Path

from build_d6_current_residue_manifest import ROOT, build_manifest, sha256, stable_hash
from verify_d6_current_residue_manifest import verify


MANIFEST = ROOT / "d6_current_residue_manifest.json"
VERIFIER = ROOT / "verify_d6_current_residue_manifest.py"


class CurrentResidueManifestControls(unittest.TestCase):
    def test_fresh_rebuild_is_identical(self) -> None:
        frozen = json.loads(MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(build_manifest(), frozen)
        self.assertEqual(
            sha256(MANIFEST),
            "d06cdb23b03b1f8523fbe3c3ede5bc297ace21a844282830c5cf51a6cee3512d",
        )

    def test_independent_checker_passes(self) -> None:
        result = verify(MANIFEST)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["K7"], 12_839)
        self.assertEqual(result["K6_only"], 990)
        self.assertEqual(result["combined"], 13_829)
        self.assertEqual(result["layers_checked"], 9)
        self.assertEqual(result["pairwise_overlaps_checked"], 36)

    def test_final_indices_are_explicit_disjoint_and_hash_bound(self) -> None:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        k7 = manifest["classes"]["K7"]
        k6 = manifest["classes"]["K6_only"]
        self.assertEqual(len(k7["final_indices"]), k7["final_count"])
        self.assertEqual(len(k6["final_indices"]), k6["final_count"])
        self.assertFalse(set(k7["final_indices"]) & set(k6["final_indices"]))
        self.assertEqual(stable_hash(k7["final_indices"]), k7["final_indices_sha256"])
        self.assertEqual(stable_hash(k6["final_indices"]), k6["final_indices_sha256"])
        combined = sorted(set(k7["final_indices"]) | set(k6["final_indices"]))
        self.assertEqual(combined, manifest["combined"]["indices"])
        self.assertEqual(stable_hash(combined), manifest["combined"]["indices_sha256"])

    def test_layer_counts_and_overlaps(self) -> None:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        layers = {(row["population"], row["layer"]): row for row in manifest["layers"]}
        self.assertEqual(layers[("K7", "positive_dual_full")]["incremental_rejections"], 102)
        self.assertEqual(layers[("K7", "positive_dual_full")]["overlap_with_prior_layers"], 0)
        self.assertEqual(layers[("K7", "reflection_overlap_full")]["incremental_rejections"], 0)
        self.assertEqual(layers[("K6-only", "normal_inertia_full")]["incremental_rejections"], 107)
        self.assertEqual(layers[("K6-only", "normal_inertia_full")]["overlap_with_prior_layers"], 0)
        self.assertEqual(len(manifest["pairwise_layer_overlaps"]), 36)

    def test_checker_does_not_import_manifest_builder(self) -> None:
        syntax = ast.parse(VERIFIER.read_text(encoding="utf-8"))
        imports = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest", imports)


if __name__ == "__main__":
    unittest.main(verbosity=2)
