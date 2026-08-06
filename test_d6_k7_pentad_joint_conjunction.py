#!/usr/bin/env python3
"""Focused controls for the exact K7 support-or-pentad conjunction."""

from __future__ import annotations

import ast
import csv
import gzip
import json
import unittest

import build_d6_k7_pentad_joint_conjunction as builder
import verify_d6_k7_pentad_joint_conjunction as checker
from verify_profile_d6 import lower_bound_18_graph


class PentadJointConjunctionTests(unittest.TestCase):
    def test_all_root_artifact_hashes_are_pinned_and_present(self) -> None:
        self.assertEqual(builder.EXPECTED_HASHES, checker.EXPECTED_HASHES)
        for name, expected in builder.EXPECTED_HASHES.items():
            self.assertEqual(builder.sha256(builder.ROOT / name), expected, name)

    def test_v2_manifest_embeds_exact_ordered_258_graph_corpus(self) -> None:
        payload = json.loads(builder.CURRENT_RESIDUE.read_text(encoding="utf-8"))
        k7 = payload["classes"]["K7"]
        indices = [int(value) for value in k7["residue_indices"]]
        graphs = list(k7["graphs"])
        self.assertEqual(len(indices), 258)
        self.assertEqual([int(graph["index"]) for graph in graphs], indices)
        self.assertEqual(builder.stable_hash(indices), k7["residue_indices_sha256"])
        self.assertEqual(builder.stable_hash(graphs), k7["graphs_sha256"])

    def test_independent_enumerator_includes_cap_four_to_seven(self) -> None:
        # With no L-edges and four eligible vertices, all 16 subsets are
        # covers.  The older support-only verifier intentionally specialized
        # to |Z|<=3; this theorem checker must also count the direct cap-four
        # failure at Z=1111.
        covers = checker.independently_all_eligible_covers((0, 0, 0, 0), 0b1111)
        self.assertEqual(covers, tuple(range(16)))
        self.assertIn(0b1111, covers)

    def test_pentad_rejection_keys_are_36_covers_on_34_graphs(self) -> None:
        report = json.loads(builder.PENTAD_REPORT.read_text(encoding="utf-8"))
        keys = [
            builder.key_tuple(record)
            for record in report["covers"] if record["status"] == "REJECTED"
        ]
        self.assertEqual(len(keys), 36)
        self.assertEqual(len(set(keys)), 36)
        self.assertEqual(len({key[0] for key in keys}), 34)
        self.assertTrue(all(
            record.get("certificate") is not None
            for record in report["covers"] if record["status"] == "REJECTED"
        ))

    def test_first_graph_independent_cover_partition_matches_archive(self) -> None:
        payload = json.loads(builder.CURRENT_RESIDUE.read_text(encoding="utf-8"))
        graph = payload["classes"]["K7"]["graphs"][0]
        index = int(graph["index"])
        selected = {index}
        prior = checker.read_witness_keys(
            checker.PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
        )
        tetrad = checker.read_witness_keys(
            checker.TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
        )
        result = checker.verify_graph((
            graph,
            [[list(seed), zmask] for seed, zmask in prior.get(index, set())],
            [[list(seed), zmask] for seed, zmask in tetrad.get(index, set())],
        ))
        with gzip.open(
            checker.JOINT_DECISIONS, "rt", encoding="utf-8", newline=""
        ) as stream:
            row = next(csv.DictReader(stream, delimiter="\t"))
        observed = result["counts"]
        self.assertEqual(int(row["index"]), index)
        self.assertEqual(observed["seeds"], int(row["seeds"]))
        self.assertEqual(observed["covers"], int(row["covers"]))
        self.assertEqual(
            observed["current_passing_covers"],
            int(row["current_passing_covers"]),
        )
        self.assertEqual(
            observed.get("joint_pre_capacity_failing_covers", 0),
            int(row["joint_pre_capacity_failing_covers"]),
        )
        self.assertEqual(
            observed["pre_capacity_passing_covers"],
            int(row["pre_capacity_passing_covers"]),
        )

    def test_checker_has_no_builder_import(self) -> None:
        tree = ast.parse(
            (checker.ROOT / "verify_d6_k7_pentad_joint_conjunction.py")
            .read_text(encoding="utf-8")
        )
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module)
        self.assertNotIn("build_d6_k7_pentad_joint_conjunction", imported)

    def test_positive_18_control_has_no_k7(self) -> None:
        graph = tuple(int(row) for row in lower_bound_18_graph())
        self.assertEqual(sum(1 for _ in builder.reference.clique_masks(graph, 7)), 0)


if __name__ == "__main__":
    unittest.main()
