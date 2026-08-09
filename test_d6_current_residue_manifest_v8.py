#!/usr/bin/env python3
"""Focused controls for the source-bound dimension-six v8 residue manifest."""

from __future__ import annotations

import ast
import copy
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_current_residue_manifest_v8 as builder
import verify_d6_current_residue_manifest_v8 as checker


class CurrentResidueV8Controls(unittest.TestCase):
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
            "porcelain_sha256": "0" * 64,
            "package_sources_equal_committed_blobs": False,
        }
        cls.manifest = builder.build_manifest(package)
        cls.v7 = json.loads(
            (builder.ROOT / "d6_current_residue_manifest_v7.json").read_text(
                encoding="utf-8"
            )
        )
        cls.temporary = tempfile.TemporaryDirectory()
        cls.path = Path(cls.temporary.name) / "v8.json"
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

    def test_exact_delta_and_unchanged_K7_class(self) -> None:
        k7 = self.manifest["classes"]["K7"]
        k6 = self.manifest["classes"]["K6_only"]
        self.assertEqual(k7, self.v7["classes"]["K7"])
        self.assertEqual(k7["count"], 12)
        self.assertEqual(k7["indices_sha256"], builder.EXPECTED_V7_K7_SHA256)

        v7_indices = self.v7["classes"]["K6_only"]["indices"]
        expected = [
            index for index in v7_indices
            if index not in set(builder.EXPECTED_REJECTED)
        ]
        self.assertEqual(k6["v7_count"], 251)
        self.assertEqual(k6["post_v7_exact_rejection_count"], 2)
        self.assertEqual(
            k6["post_v7_exact_rejected_indices"], builder.EXPECTED_REJECTED
        )
        self.assertEqual(k6["indices"], expected)
        self.assertEqual(k6["count"], 249)
        self.assertEqual(k6["indices_sha256"], builder.EXPECTED_K6_SHA256)
        self.assertEqual(
            k6["graphs_sha256"], builder.EXPECTED_K6_GRAPHS_SHA256
        )
        self.assertNotIn(3_138_618, k6["indices"])
        self.assertNotIn(3_673_988, k6["indices"])

    def test_sequential_layer_and_cumulative_accounting(self) -> None:
        k6 = self.manifest["classes"]["K6_only"]
        layer = k6["exact_layers"][-1]
        self.assertEqual(
            layer,
            {
                "name": "saturated_singleton_basis",
                "input_count": 251,
                "input_indices_sha256": builder.EXPECTED_V7_K6_SHA256,
                "rejection_count": 2,
                "rejected_indices": builder.EXPECTED_REJECTED,
                "rejected_indices_sha256": builder.EXPECTED_REJECTED_SHA256,
                "residue_count": 249,
                "residue_indices_sha256": builder.EXPECTED_K6_SHA256,
            },
        )
        self.assertEqual(k6["post_v6_exact_rejection_count"], 376)
        self.assertEqual(
            k6["post_v6_exact_rejected_indices"][-2:],
            builder.EXPECTED_REJECTED,
        )
        self.assertEqual(
            k6["post_v6_exact_rejected_indices_sha256"],
            builder.EXPECTED_POST_V6_K6_REJECTED_SHA256,
        )

        combined = self.manifest["combined"]
        self.assertEqual(combined["v7_count"], 263)
        self.assertEqual(combined["post_v7_exact_rejection_count"], 2)
        self.assertEqual(combined["count"], 261)
        self.assertEqual(combined["class_counts"], {"K7": 12, "K6_only": 249})
        self.assertEqual(
            combined["ordered_indices_sha256"], builder.EXPECTED_COMBINED_SHA256
        )
        self.assertEqual(
            combined["sorted_indices_sha256"], builder.EXPECTED_SORTED_SHA256
        )

    def test_trust_tiers_controls_and_source_roots(self) -> None:
        accounting = self.manifest["certificate_accounting"]
        self.assertEqual(accounting["exact"]["total_rejections_from_v5"], 645)
        self.assertEqual(accounting["exact"]["post_v6_rejections"], 381)
        self.assertEqual(
            accounting["exact"]["layers"]["K6_saturated_singleton_basis"], 2
        )
        self.assertFalse(accounting["exact"]["floating_point_enters_rejection"])
        self.assertEqual(
            accounting["interval"]["total_incremental_rejections_from_v5"], 5
        )
        self.assertEqual(accounting["union_rejections_from_v5"], 650)
        self.assertEqual(accounting["current_residue"], 261)

        controls = self.manifest["controls"]
        self.assertEqual(
            controls["known_realizable_18"][
                "exact_K6_saturated_singleton_basis"
            ],
            {"passed": True, "K6_seeds": 32},
        )
        self.assertTrue(
            controls["exact_K6_saturated_singleton_synthetic"][
                "odd_opposite_sign_cycle_fails"
            ]
        )
        boundary = self.manifest["source_boundary"]
        self.assertEqual(
            boundary["primary_artifacts"],
            dict(sorted(builder.PRIMARY_ARTIFACTS.items())),
        )
        self.assertEqual(
            set(boundary["layer_source_roots"]),
            {"inherited_v7", "exact_K6_saturated_singleton_basis"},
        )
        self.assertEqual(
            boundary["layer_source_roots"][
                "exact_K6_saturated_singleton_basis"
            ]["launch_commit"],
            "87c58d612eb17959e43bf920057fcff95212776c",
        )

    def test_import_independent_checker_reconstructs_everything(self) -> None:
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
                "K6_only": 249,
                "combined": 261,
                "exact_rejections_from_v5": 645,
                "interval_incremental_rejections": 5,
            },
        )
        self.assertEqual(
            result["new_exact_layer"]["rejected_indices"],
            builder.EXPECTED_REJECTED,
        )
        self.assertTrue(all(result["checks"].values()))

    def test_graph_delta_accounting_and_evidence_tampering_rejected(self) -> None:
        mutations = []

        k7 = copy.deepcopy(self.manifest)
        k7["classes"]["K7"]["indices"][0] = 3_138_618
        mutations.append(k7)

        graph = copy.deepcopy(self.manifest)
        removed_record = next(
            record
            for record in self.v7["classes"]["K6_only"]["graphs"]
            if record["index"] == 3_138_618
        )
        graph["classes"]["K6_only"]["graphs"].append(removed_record)
        mutations.append(graph)

        layer = copy.deepcopy(self.manifest)
        layer["classes"]["K6_only"]["exact_layers"][-1][
            "rejected_indices"
        ] = [3_138_618]
        mutations.append(layer)

        trust = copy.deepcopy(self.manifest)
        trust["certificate_accounting"]["interval"]["trust_assumptions"][
            "transcendentals"
        ] = "none"
        mutations.append(trust)

        source = copy.deepcopy(self.manifest)
        source["source_boundary"]["layer_source_roots"][
            "exact_K6_saturated_singleton_basis"
        ]["launch_commit"] = "0" * 40
        mutations.append(source)

        for number, payload in enumerate(mutations):
            path = Path(self.temporary.name) / f"tamper-{number}.json"
            self.write(payload, path)
            with self.subTest(number=number):
                with self.assertRaises(ValueError):
                    checker.verify(path, None, enforce_source_boundary=False)

    def test_checker_imports_no_builder_or_production_module(self) -> None:
        source = (
            checker.ROOT / "verify_d6_current_residue_manifest_v8.py"
        ).read_text(encoding="utf-8")
        tree = ast.parse(source)
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_current_residue_manifest_v8", imports)
        self.assertFalse(any(
            name.startswith(
                ("d6_k6_", "d6_k7_", "verify_d6_k6_", "verify_d6_k7_")
            )
            for name in imports
        ))


if __name__ == "__main__":
    unittest.main(verbosity=2)
