#!/usr/bin/env python3
"""Controls for the exact v5 rooted 18-deletion corpus."""

from __future__ import annotations

import ast
import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import build_d6_residue_18_deletions_v2 as builder
import verify_d6_residue_18_deletions_v2 as checker


class Residue18DeletionV2Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        executable = shutil.which("labelg")
        if executable is None:
            raise unittest.SkipTest("labelg is required")
        cls.labelg = Path(executable).resolve()
        cls.manifest = builder.build_manifest()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.manifest_path = Path(cls.temporary.name) / "deletions-v2.json"
        cls.write(cls.manifest, cls.manifest_path)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    @staticmethod
    def write(payload: object, path: Path) -> None:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def test_source_and_checker_hash_pins(self) -> None:
        for name, expected in builder.EXPECTED_FILES.items():
            self.assertEqual(builder.sha256(builder.ROOT / name), expected, name)
        self.assertEqual(
            checker.file_sha256(checker.ROOT / checker.BUILDER),
            checker.EXPECTED_BUILDER_SHA256,
        )

    def test_exact_counts_and_rooted_partition(self) -> None:
        summary = self.manifest["summary"]
        for key, expected in builder.EXPECTED_COUNTS.items():
            self.assertEqual(summary[key], expected, key)
        self.assertEqual(
            summary["parent_class_counts"], {"K7": 155, "K6_only": 756}
        )
        self.assertEqual(summary["unique_clique_classes"]["contains_K7"], 1616)
        self.assertEqual(
            summary["unique_clique_classes"]["K6_without_K7"], 10359
        )
        self.assertEqual(summary["unique_clique_classes"]["no_K6"], 0)
        occurrence_keys = {
            (
                occurrence["parent_class"],
                occurrence["parent_index"],
                occurrence["deleted_vertex"],
            )
            for record in self.manifest["unique_deletions"]
            for occurrence in record["occurrences"]
        }
        self.assertEqual(len(occurrence_keys), 17_309)
        self.assertTrue(
            all(
                0 <= occurrence["attachment_mask"] < (1 << 18)
                for record in self.manifest["unique_deletions"]
                for occurrence in record["occurrences"]
            )
        )

    def test_v1_subset_losses_and_standard_matches(self) -> None:
        v1 = json.loads(
            (builder.ROOT / builder.V1_NAME).read_text(encoding="utf-8")
        )
        self.assertEqual(
            v1["summary"]["deletion_occurrences"]
            - self.manifest["summary"]["deletion_occurrences"],
            49 * 19,
        )
        self.assertEqual(
            v1["summary"]["unique_deletions"]
            - self.manifest["summary"]["unique_deletions"],
            737,
        )
        for key in (
            "standard_support_deletion_occurrences",
            "standard_support_parent_graphs",
            "standard_compatible_unique_deletions",
            "standard_compatible_deletion_occurrences",
            "standard_compatible_parent_graphs",
        ):
            self.assertEqual(self.manifest["summary"][key], v1["summary"][key])

    def test_independent_checker_reconstructs_every_occurrence(self) -> None:
        result = checker.verify(self.manifest_path, self.labelg, None)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["counts"], builder.EXPECTED_COUNTS)
        self.assertTrue(all(result["checks"].values()))

    def test_explicit_hash_boundary_rejects_tamper_immediately(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["unique_deletions"][0]["occurrences"][0]["attachment_mask"] ^= 1
        path = Path(self.temporary.name) / "tampered.json"
        self.write(payload, path)
        expected = checker.file_sha256(self.manifest_path)
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            checker.verify(path, self.labelg, expected)

    def test_checker_imports_neither_v2_builder_nor_geometry_kernel(self) -> None:
        tree = ast.parse(
            (checker.ROOT / "verify_d6_residue_18_deletions_v2.py").read_text(
                encoding="utf-8"
            )
        )
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_residue_18_deletions_v2", imports)
        self.assertFalse(any(name.startswith("d6_k6_") for name in imports))
        self.assertFalse(any(name.startswith("d6_k7_") for name in imports))


if __name__ == "__main__":
    unittest.main(verbosity=2)
