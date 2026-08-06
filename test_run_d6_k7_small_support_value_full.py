#!/usr/bin/env python3
"""Production-wrapper controls for the full K7 sparse-value campaign."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import run_d6_k7_small_support_value_full as runner


ROOT = Path(__file__).resolve().parent


class ProvenanceControls(unittest.TestCase):
    def test_dependency_hash_pins(self):
        self.assertEqual(
            runner.verify_dependencies(), runner.EXPECTED_DEPENDENCY_SHA256
        )

    def test_full_selection_partition(self):
        selected, provenance = runner.load_selection(
            ROOT / ".runs/d6_k7_rank_survivors.json",
            ROOT / "d6_k7_support_rank_survivors_decisions.tsv.gz",
            ROOT / "d6_k7_support_rank_survivors_report.json",
            ROOT / "d6_k7_support_rank_survivors_checkpoint.json",
        )
        self.assertEqual(len(selected), runner.EXPECTED_SELECTED)
        self.assertEqual(
            provenance["rank_survivors"], runner.EXPECTED_RANK_SURVIVORS
        )
        self.assertEqual(
            provenance["support_rejected"], runner.EXPECTED_SUPPORT_REJECTED
        )
        self.assertEqual(len({int(graph["index"]) for graph in selected}),
                         len(selected))


class DecisionControls(unittest.TestCase):
    def test_flatten_and_validate(self):
        result = {
            "index": 17,
            "decision": "REJECTED",
            "first_failing_seed": 127,
            "counts": {
                "seeds": 1,
                "sparse_value_failures": 3,
                "sparse_value_branch:forbidden_two_defect_cycle": 3,
            },
        }
        row = runner.flatten_result(result)
        self.assertEqual(runner.validate_row(row, 17), tuple(map(str, row)))
        self.assertEqual(row[1], "REJECTED")
        position = runner.DECISION_FIELDS.index(
            "sparse_value_branch_forbidden_two_defect_cycle"
        )
        self.assertEqual(row[position], 3)

    def test_committed_prefix_discards_torn_tail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "partial.tsv"
            first = [0] * len(runner.DECISION_FIELDS)
            first[0:3] = [10, "SURVIVOR", 0]
            second = [0] * len(runner.DECISION_FIELDS)
            second[0:3] = [20, "SURVIVOR", 0]
            runner.atomic_decisions(path, (first, second))
            completed = runner.load_completed(path, (10, 20), 1)
            self.assertEqual(len(completed), 1)
            with path.open(encoding="ascii") as stream:
                self.assertEqual(sum(1 for _ in stream), 2)  # header + row

    def test_deterministic_compact_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source = base / "decisions.tsv"
            source.write_bytes(b"index\tdecision\n1\tSURVIVOR\n")
            first = base / "first.tsv.gz"
            second = base / "second.tsv.gz"
            runner.deterministic_gzip(source, first)
            runner.deterministic_gzip(source, second)
            self.assertEqual(runner.sha256(first), runner.sha256(second))
            with gzip.open(first, "rb") as stream:
                self.assertEqual(stream.read(), source.read_bytes())


if __name__ == "__main__":
    unittest.main()

