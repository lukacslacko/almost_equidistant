#!/usr/bin/env python3
"""Focused tests for the source-only class-2100 interval increment."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import build_d6_interval_18_class2100_increment as producer
import verify_d6_interval_18_class2100_increment as checker


ROOT = Path(__file__).resolve().parent
LIVE_CAMPAIGN = ROOT / ".runs/d6_interval_18_cover_v7_cap100000_w3/campaign.json"
LIVE_RESULT = (
    ROOT
    / ".runs/d6_interval_18_cover_v7_cap100000_w3/results/"
    "result_000032_0002100.json"
)


class StructuralTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.deletions = json.loads(
            (ROOT / "d6_residue_18_deletions_v2.json").read_text(encoding="utf-8")
        )
        cls.record = cls.deletions["unique_deletions"][checker.TARGET_CLASS_INDEX]
        cls.rows = checker.validate_rows(cls.record["adjacency"], 18, "target")

    def test_graph6_and_class_id(self) -> None:
        self.assertEqual(checker.decode_graph6(checker.TARGET_GRAPH6), self.rows)
        self.assertEqual(producer.decode_graph6(producer.TARGET_GRAPH6), self.rows)
        self.assertEqual(self.record["class_id"], checker.TARGET_CLASS_ID)

    def test_independent_order_generator(self) -> None:
        orders = checker.generate_orders(self.rows, 4)
        self.assertEqual(len(orders), 2)
        self.assertEqual(
            [list(seed) for seed, _order in orders],
            [[1, 7, 8, 13, 14, 15, 16]] * 2,
        )
        self.assertEqual(
            [list(order) for _seed, order in orders],
            [
                [5, 12, 17, 6, 11, 10, 9, 4, 3, 2, 0],
                [12, 6, 17, 11, 10, 9, 5, 4, 3, 2, 0],
            ],
        )
        self.assertEqual(
            [checker.circle_positions(self.rows, *item) for item in orders],
            [[0], [0]],
        )

    def test_induced_parent_occurrences(self) -> None:
        v7 = checker.load_residue(
            "d6_current_residue_manifest_v7.json",
            checker.EXPECTED_V7_COUNTS,
            checker.EXPECTED_V7_CLASS_HASHES,
            checker.EXPECTED_V7_COMBINED_HASH,
        )
        audit = checker.reconstruct_occurrences(v7, self.rows, self.record)
        self.assertEqual(audit["edge_and_degree_candidates"], 10)
        self.assertEqual(audit["unique_parent_indices"], [3_945_564])
        self.assertEqual(
            [item["deleted_vertex"] for item in audit["exact_isomorphism_occurrences"]],
            [4, 15],
        )

    def test_independent_v8_containment_scan(self) -> None:
        v8 = checker.load_residue(
            "d6_current_residue_manifest_v8.json",
            checker.EXPECTED_V8_COUNTS,
            checker.EXPECTED_V8_CLASS_HASHES,
            checker.EXPECTED_V8_COMBINED_HASH,
        )
        audit = checker.containment_scan(v8, self.rows)
        self.assertEqual(audit["parents_scanned"], 261)
        self.assertEqual(audit["deletions_scanned"], 4_959)
        self.assertEqual(audit["eligible_deletions"], 1_085)
        self.assertEqual(audit["deleted_complement_edge_variants"], 17_344)
        self.assertEqual(audit["exact_isomorphism_comparisons"], 6)
        self.assertEqual(audit["rejected_parent_indices"], [3_945_564])


@unittest.skipUnless(
    LIVE_CAMPAIGN.exists() and LIVE_RESULT.exists(),
    "live immutable class-2100 checkpoint is not present",
)
class LiveEvidenceIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temporary.name)
        cls.evidence = cls.directory / "evidence.json.gz"
        cls.report = cls.directory / "report.json"
        producer.build_increment(
            LIVE_CAMPAIGN,
            LIVE_RESULT,
            cls.evidence,
            cls.report,
            source_boundary={"status": "TEST_ONLY"},
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_raw_hashes_and_deterministic_evidence(self) -> None:
        self.assertEqual(producer.sha256(LIVE_CAMPAIGN), producer.TARGET_RAW_CAMPAIGN_SHA256)
        self.assertEqual(producer.sha256(LIVE_RESULT), producer.TARGET_RAW_RESULT_SHA256)
        first = self.evidence.read_bytes()
        rebuilt = producer.make_evidence(LIVE_CAMPAIGN.read_bytes(), LIVE_RESULT.read_bytes())
        self.assertEqual(first, rebuilt)

    def test_import_independent_structural_verification(self) -> None:
        result = checker.verify_increment(
            self.report,
            self.evidence,
            replay=False,
            enforce_source_boundary=False,
        )
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["conclusion"]["rejected_parent_indices"], [3_945_564])
        self.assertEqual(result["interval_replay"]["status"], "SKIPPED")

    def test_report_tamper_fails_closed(self) -> None:
        directory = self.directory / "tampered_report"
        directory.mkdir()
        evidence = directory / self.evidence.name
        shutil.copyfile(self.evidence, evidence)
        value = json.loads(self.report.read_text(encoding="utf-8"))
        value["conclusion"]["rejected_parent_indices"] = [3_945_565]
        report = directory / "report.json"
        report.write_text(json.dumps(value), encoding="utf-8")
        with self.assertRaises(AssertionError):
            checker.verify_increment(
                report,
                evidence,
                replay=False,
                enforce_source_boundary=False,
            )

    def test_evidence_tamper_fails_closed(self) -> None:
        directory = self.directory / "tampered_evidence"
        directory.mkdir()
        report = directory / "report.json"
        shutil.copyfile(self.report, report)
        evidence = directory / self.evidence.name
        payload = bytearray(self.evidence.read_bytes())
        payload[-1] ^= 1
        evidence.write_bytes(payload)
        with self.assertRaises(AssertionError):
            checker.verify_increment(
                report,
                evidence,
                replay=False,
                enforce_source_boundary=False,
            )

    @unittest.skipUnless(
        os.environ.get("D6_INTERVAL_REPLAY") == "1",
        "set D6_INTERVAL_REPLAY=1 for the approximately 30-second kernel replay",
    )
    def test_full_kernel_and_control_replay(self) -> None:
        result = checker.verify_increment(
            self.report,
            self.evidence,
            replay=True,
            enforce_source_boundary=False,
        )
        self.assertEqual(result["interval_replay"]["status"], "PASS")
        self.assertEqual(result["interval_replay"]["kernel_calls"], 31)
        self.assertEqual(result["controls"]["slice_control_records_replayed"], 48)


if __name__ == "__main__":
    unittest.main()
