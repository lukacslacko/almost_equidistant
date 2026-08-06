#!/usr/bin/env python3
"""Controls for the self-contained exact dimension-six residue manifest v3."""

from __future__ import annotations

import ast
import copy
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_current_residue_manifest_v3 as builder
import verify_d6_current_residue_manifest_v3 as checker


class CurrentResidueManifestV3Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = builder.build_manifest()
        cls.directory = tempfile.TemporaryDirectory()
        cls.path = Path(cls.directory.name) / "v3.json"
        cls.write(cls.manifest, cls.path)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.directory.cleanup()

    @staticmethod
    def write(payload: dict, path: Path) -> None:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def verify_tamper(self, payload: dict) -> None:
        path = Path(self.directory.name) / "tampered.json"
        self.write(payload, path)
        with self.assertRaises(ValueError):
            checker.verify(path)

    def test_upstream_commits_and_all_hash_pins(self) -> None:
        self.assertEqual(builder.EXPECTED_FILES, checker.EXPECTED_FILES)
        self.assertEqual(
            builder.K7_RESULT_COMMIT,
            "01c5b18413b6bab5225558fa950971efeda58e0d",
        )
        self.assertEqual(
            builder.K6_RESULT_COMMIT,
            "a21db74bacdf4c9c5c841ac137af326b331d8137",
        )
        for name, expected in builder.EXPECTED_FILES.items():
            self.assertEqual(builder.sha256(builder.ROOT / name), expected, name)
        self.assertEqual(
            builder.sha256(builder.ROOT / checker.BUILDER),
            checker.EXPECTED_BUILDER_SHA256,
        )

    def test_build_is_deterministic_and_has_exact_counts(self) -> None:
        self.assertEqual(builder.build_manifest(), self.manifest)
        self.assertEqual(self.manifest["class_order"], ["K7", "K6_only"])
        self.assertEqual(self.manifest["classes"]["K7"]["residue_count"], 155)
        self.assertEqual(self.manifest["classes"]["K6_only"]["residue_count"], 822)
        self.assertEqual(self.manifest["combined"]["count"], 977)
        self.assertEqual(
            self.manifest["combined"]["ordered_indices_sha256"],
            "71af8031e783b8a87d710156cae9cc4a19c575190a7c164a65935061f8ac1d7d",
        )

    def test_independent_checker_rederives_every_boundary(self) -> None:
        result = checker.verify(self.path)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["classes"]["K7"]["residue"], 155)
        self.assertEqual(result["classes"]["K6_only"]["residue"], 822)
        self.assertEqual(result["classes"]["K6_only"]["stage_rejections"], [116, 30, 9])
        self.assertEqual(result["combined"]["count"], 977)
        self.assertTrue(all(result["checks"].values()))

    def test_embedded_graphs_are_ordered_and_structurally_classified(self) -> None:
        for class_name in self.manifest["class_order"]:
            payload = self.manifest["classes"][class_name]
            graphs = payload["graphs"]
            self.assertEqual(
                [record["index"] for record in graphs], payload["residue_indices"]
            )
            self.assertEqual(builder.stable_hash(graphs), payload["graphs_sha256"])
            for record in graphs:
                rows = builder.validate_adjacency(
                    record["adjacency"], 19, f"test {class_name}"
                )
                self.assertTrue(builder.alpha_at_most_two(rows))
                if class_name == "K7":
                    self.assertTrue(builder.contains_clique(rows, 7))
                else:
                    self.assertTrue(builder.contains_clique(rows, 6))
                    self.assertFalse(builder.contains_clique(rows, 7))

    def test_exact_stage_partitions(self) -> None:
        k7 = self.manifest["classes"]["K7"]
        self.assertEqual(k7["v2_boundary_count"], 258)
        self.assertEqual(k7["exact_rejections_since_v2"], 103)
        self.assertEqual(
            len(k7["exact_rejected_indices"]) + len(k7["residue_indices"]), 258
        )
        k6 = self.manifest["classes"]["K6_only"]
        self.assertEqual(k6["v2_boundary_count"], 977)
        self.assertEqual(
            [stage["rejection_count"] for stage in k6["stage_rejections"]],
            [116, 30, 9],
        )
        self.assertEqual(k6["exact_rejections_since_v2"], 155)
        self.assertEqual(
            len(k6["exact_rejected_indices"]) + len(k6["residue_indices"]), 977
        )

    def test_positive18_is_separate_and_exact(self) -> None:
        control = self.manifest["positive_18_control"]
        self.assertIs(control["in_level_19_corpus"], False)
        self.assertEqual(control["vertices"], 18)
        self.assertEqual(control["K6_seeds"], 32)
        self.assertEqual(control["K7_seeds"], 0)
        self.assertEqual(
            control["adjacency_sha256"],
            "7697e049810093251c017328fa1043647cc3bf12df7c61ecd5433225154497b5",
        )
        corpus_rows = sum(
            len(self.manifest["classes"][name]["graphs"])
            for name in self.manifest["class_order"]
        )
        self.assertEqual(corpus_rows, 977)

    def test_partition_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["classes"]["K7"]["residue_indices"].pop()
        self.verify_tamper(payload)

    def test_adjacency_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["classes"]["K6_only"]["graphs"][0]["adjacency"][0] ^= 1 << 1
        self.verify_tamper(payload)

    def test_positive_control_cannot_enter_corpus(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["positive_18_control"]["in_level_19_corpus"] = True
        self.verify_tamper(payload)

    def test_checker_imports_no_builder_or_production_kernel(self) -> None:
        syntax = ast.parse(
            (checker.ROOT / "verify_d6_current_residue_manifest_v3.py")
            .read_text(encoding="utf-8")
        )
        imports = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest_v3", imports)
        self.assertFalse(any(name.startswith("d6_k7_") for name in imports))
        self.assertFalse(any(name.startswith("d6_k6_") for name in imports))


if __name__ == "__main__":
    unittest.main(verbosity=2)
