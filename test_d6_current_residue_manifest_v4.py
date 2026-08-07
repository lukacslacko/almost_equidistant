#!/usr/bin/env python3
"""Controls for the exact dimension-six v4 residue boundary."""

from __future__ import annotations

import ast
import copy
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_current_residue_manifest_v4 as builder
import verify_d6_current_residue_manifest_v4 as checker


class CurrentResidueV4Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = builder.build_manifest()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.path = Path(cls.temporary.name) / "v4.json"
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

    def test_all_source_and_result_hashes_are_pinned(self) -> None:
        for name, expected in builder.EXPECTED_FILES.items():
            self.assertEqual(builder.sha256(builder.ROOT / name), expected, name)
        self.assertEqual(
            checker.sha256(checker.ROOT / checker.BUILDER),
            checker.EXPECTED_BUILDER_SHA256,
        )

    def test_exact_counts_and_rejection_partition(self) -> None:
        self.assertEqual(self.manifest["classes"]["K7"]["count"], 155)
        k6 = self.manifest["classes"]["K6_only"]
        self.assertEqual(k6["v3_count"], 822)
        self.assertEqual(k6["empty_support_rejection_count"], 17)
        self.assertEqual(k6["count"], 805)
        self.assertEqual(self.manifest["combined"]["count"], 960)
        self.assertEqual(
            set(k6["empty_support_rejected_indices"])
            & {record["index"] for record in k6["graphs"]},
            set(),
        )

    def test_independent_checker_rebuilds_v4(self) -> None:
        result = checker.verify(self.path, None)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["counts"], {"K7": 155, "K6_only": 805, "combined": 960, "rejected": 17})
        self.assertTrue(all(result["checks"].values()))

    def test_graph_or_rejection_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["classes"]["K6_only"]["graphs"][0]["adjacency"][0] ^= 2
        path = Path(self.temporary.name) / "tampered.json"
        self.write(payload, path)
        with self.assertRaises(ValueError):
            checker.verify(path, None)

    def test_checker_imports_no_builder_or_k6_kernel(self) -> None:
        tree = ast.parse(
            (checker.ROOT / "verify_d6_current_residue_manifest_v4.py").read_text(
                encoding="utf-8"
            )
        )
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest_v4", imports)
        self.assertFalse(any(name.startswith("d6_k6_") for name in imports))


if __name__ == "__main__":
    unittest.main(verbosity=2)
