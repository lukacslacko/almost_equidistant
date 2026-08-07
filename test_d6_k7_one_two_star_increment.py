#!/usr/bin/env python3
"""Corpus and independence controls for the K7 one/two-free star layer."""

from __future__ import annotations

import ast
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import build_d6_k7_one_two_star_increment as builder
import verify_d6_k7_one_two_star_increment as verifier


def inputs() -> tuple[list[dict], dict[int, dict], list[int]]:
    manifest, current = builder.validate_upstream()
    graphs = manifest["classes"]["K7"]["graphs"]
    current_by_index = {
        int(record["index"]): record for record in current["records"]
    }
    indices = [int(graph["index"]) for graph in graphs]
    return graphs, current_by_index, indices


class StarIncrementCorpusControls(unittest.TestCase):
    def test_complete_19_graph_builder_replay_has_exact_three_graph_marginal(self) -> None:
        graphs, current_by_index, indices = inputs()
        records = [
            builder.evaluate_graph((graph, current_by_index[int(graph["index"])]))
            for graph in graphs
        ]
        rejected = [
            record["index"]
            for record in records
            if record["decision"] == "REJECTED"
        ]
        totals: Counter[str] = Counter()
        for record in records:
            totals.update(record["counts"])
        self.assertEqual(indices, builder.json.loads(
            builder.MANIFEST.read_text(encoding="utf-8")
        )["classes"]["K7"]["indices"])
        self.assertEqual(builder.stable_hash(indices), builder.EXPECTED_INPUT_SHA256)
        self.assertEqual(rejected, builder.EXPECTED_REJECTIONS)
        self.assertEqual(totals["current_passing_families"], 88)
        self.assertEqual(totals["star_failed_families"], 4)
        self.assertEqual(sum(len(record["certificates"]) for record in records), 4)

    def test_independent_full_quantifier_rebuild_matches_rejecting_graph(self) -> None:
        graphs, current_by_index, _indices = inputs()
        graph = next(item for item in graphs if int(item["index"]) == 2592657)
        index = int(graph["index"])
        selected = {index}
        prior = verifier.current_independent.base_independent.read_failure_keys(
            verifier.current_independent.base_independent.PRIOR_CERTIFICATES,
            selected,
            "dual_failure_witnesses",
        )
        tetrad = verifier.current_independent.base_independent.read_failure_keys(
            verifier.current_independent.base_independent.TETRAD_CERTIFICATES,
            selected,
            "tetrad_failure_witnesses",
        )
        rows = verifier.current_independent.base_independent.read_tetrad_rows(
            selected
        )
        current_payload = (
            graph,
            [
                [list(seed), zmask]
                for seed, zmask in sorted(prior.get(index, set()))
            ],
            [
                [list(seed), zmask]
                for seed, zmask in sorted(tetrad.get(index, set()))
            ],
            int(rows[index]["tetrad_passing_covers"]),
        )
        rebuilt = verifier.evaluate_graph_from_independent_current(
            (graph, current_payload, current_by_index[index])
        )
        production = builder.evaluate_graph((graph, current_by_index[index]))
        self.assertEqual(rebuilt["star_record"], production)
        self.assertEqual(production["decision"], "REJECTED")
        self.assertGreater(rebuilt["current_layer_counts"]["eligible_covers"], 0)

    def test_known_realizable_18_control_is_not_applicable(self) -> None:
        positive = tuple(map(int, builder.lower_bound_18_graph()))
        builder.reference.validate_graph(positive)
        self.assertEqual(
            sum(1 for _ in builder.reference.clique_masks(positive, 7)), 0
        )


class IndependentVerifierControls(unittest.TestCase):
    def test_checker_imports_neither_producer_nor_star_kernel(self) -> None:
        source = Path(verifier.__file__).read_text(encoding="utf-8")
        tree = ast.parse(source)
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        self.assertNotIn("build_d6_k7_one_two_star_increment", imports)
        self.assertNotIn("d6_k7_one_two_star", imports)

    def test_report_hash_tamper_is_rejected_before_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            path.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "report hash differs"):
                verifier.verify_report(
                    path,
                    "0" * 64,
                    workers=1,
                    enforce_source_boundary=False,
                )

    def test_frozen_builder_and_kernel_hashes_match_working_sources(self) -> None:
        self.assertEqual(
            verifier.EXPECTED_BUILDER_SHA256,
            verifier.sha256(verifier.ROOT / "build_d6_k7_one_two_star_increment.py"),
        )
        self.assertEqual(
            verifier.EXPECTED_KERNEL_SHA256,
            verifier.sha256(verifier.ROOT / "d6_k7_one_two_star.py"),
        )


if __name__ == "__main__":
    unittest.main()
