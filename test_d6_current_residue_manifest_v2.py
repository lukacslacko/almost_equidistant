#!/usr/bin/env python3
"""Controls for the self-contained dimension-six residue manifest v2."""

from __future__ import annotations

import ast
import copy
import json
import tempfile
import unittest
from pathlib import Path

from build_d6_current_residue_manifest_v2 import (
    ROOT,
    atomic_json,
    build_manifest,
    sha256,
    stable_hash,
)
from verify_d6_current_residue_manifest_v2 import verify


MANIFEST = ROOT / "d6_current_residue_manifest_v2.json"
VERIFICATION = ROOT / "d6_current_residue_manifest_v2_verification.json"
VERIFIER = ROOT / "verify_d6_current_residue_manifest_v2.py"


class CurrentResidueManifestV2Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))

    def tampered_path(self, payload: dict, directory: str) -> Path:
        output = Path(directory) / "tampered.json"
        output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return output

    def test_fresh_rebuild_is_byte_identical(self) -> None:
        self.assertEqual(build_manifest(), self.manifest)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / MANIFEST.name
            atomic_json(output, build_manifest())
            self.assertEqual(output.read_bytes(), MANIFEST.read_bytes())
        self.assertEqual(
            sha256(MANIFEST),
            "961dac1f9b44bb541e2c5bd1626027826eeea7bc628bf5ff0a85ab26afa949f4",
        )

    def test_independent_checker_passes(self) -> None:
        result = verify(MANIFEST)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["classes"]["K7"]["residue"], 258)
        self.assertEqual(result["classes"]["K6_only"]["residue"], 977)
        self.assertEqual(result["combined"]["count"], 1_235)
        self.assertTrue(all(result["checks"].values()))
        frozen = json.loads(VERIFICATION.read_text(encoding="utf-8"))
        self.assertEqual(result, frozen)

    def test_exact_partition_order_and_hashes(self) -> None:
        classes = self.manifest["classes"]
        k7 = classes["K7"]
        k6 = classes["K6_only"]
        self.assertEqual(self.manifest["class_order"], ["K7", "K6_only"])
        self.assertEqual(len(k7["residue_indices"]), 258)
        self.assertEqual(len(k6["residue_indices"]), 977)
        self.assertFalse(set(k7["residue_indices"]) & set(k6["residue_indices"]))
        self.assertEqual(stable_hash(k7["residue_indices"]), k7["residue_indices_sha256"])
        self.assertEqual(stable_hash(k6["residue_indices"]), k6["residue_indices_sha256"])
        self.assertEqual(
            stable_hash(k7["residue_indices"] + k6["residue_indices"]),
            self.manifest["combined"]["ordered_indices_sha256"],
        )

    def test_embedded_graphs_match_indices_and_are_well_formed(self) -> None:
        for class_name in self.manifest["class_order"]:
            payload = self.manifest["classes"][class_name]
            graphs = payload["graphs"]
            self.assertEqual(
                [record["index"] for record in graphs], payload["residue_indices"]
            )
            self.assertEqual(stable_hash(graphs), payload["graphs_sha256"])
            for record in graphs:
                adjacency = record["adjacency"]
                self.assertEqual(len(adjacency), 19)
                for left, mask in enumerate(adjacency):
                    self.assertFalse(mask & (1 << left))
                    for right in range(19):
                        self.assertEqual(
                            (mask >> right) & 1,
                            (adjacency[right] >> left) & 1,
                        )

    def test_partition_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["classes"]["K7"]["residue_indices"].pop()
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                verify(self.tampered_path(payload, directory))

    def test_adjacency_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["classes"]["K7"]["graphs"][0]["adjacency"][0] ^= 1 << 1
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                verify(self.tampered_path(payload, directory))

    def test_positive_control_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["verification_gates"]["K6_only"]["positive_control_passed"] = False
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                verify(self.tampered_path(payload, directory))

    def test_class_order_tamper_is_rejected(self) -> None:
        payload = copy.deepcopy(self.manifest)
        payload["class_order"] = ["K6_only", "K7"]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                verify(self.tampered_path(payload, directory))

    def test_checker_imports_no_builder_or_production_engine(self) -> None:
        syntax = ast.parse(VERIFIER.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest_v2", imports)
        self.assertFalse(any(name.startswith("d6_k7_") for name in imports))
        self.assertFalse(any(name.startswith("d6_k6_") for name in imports))


if __name__ == "__main__":
    unittest.main(verbosity=2)
