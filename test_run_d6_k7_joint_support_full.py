#!/usr/bin/env python3
"""Focused controls for the restartable full K7 joint-support wrapper."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import run_d6_k7_joint_support_full as runner


class FullWrapperControls(unittest.TestCase):
    def test_decision_row_keeps_joint_and_capacity_counts_separate(self):
        record = {
            "ordinal": 3,
            "index": 17,
            "status": "JOINT_PRE_CAPACITY_REJECTED",
            "error_type": "",
            "result": {"counts": {
                "seeds": 2,
                "covers": 11,
                "cover_passing": 4,
                "joint_pre_capacity_failing_covers": 3,
                "pre_capacity_passing_covers": 1,
                "capacity_passing_covers": 1,
                "capacity_incremental_failing_covers": 0,
            }},
        }
        row = dict(zip(runner.DECISION_FIELDS, runner.decision_row(record)))
        self.assertEqual(row["joint_pre_capacity_failing_covers"], 3)
        self.assertEqual(row["capacity_incremental_failing_covers"], 0)
        self.assertEqual(row["status"], "JOINT_PRE_CAPACITY_REJECTED")

    def test_certificate_archives_first_failure_for_each_current_cover(self):
        record = {
            "ordinal": 0,
            "index": 23,
            "status": "JOINT_PRE_CAPACITY_REJECTED",
            "result": {
                "pre_capacity_rejected_seeds": [7],
                "seed_records": [{
                    "seed_mask": 7,
                    "seed": [0, 1, 2, 3, 4, 5, 6],
                    "current_cover_records": [{
                        "zmask": 1,
                        "counts": {"labeled_z_families": 2},
                        "first_pre_capacity_failure": {
                            "stage": "propagation",
                            "reason": "empty_propagated_mask",
                        },
                    }],
                }],
            },
        }
        certificate = runner.certificate_entry(record)
        self.assertIsNotNone(certificate)
        failures = certificate["failing_seeds"][0]["current_cover_failures"]
        self.assertEqual(len(failures), 1)
        self.assertEqual(
            failures[0]["first_pre_capacity_failure"]["stage"],
            "propagation",
        )

    def test_checkpoint_configuration_is_strictly_bound(self):
        selected = [{"ordinal": 0, "index": 101}]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "graph_000_101.json"
            path.write_text(json.dumps({
                "schema": 1,
                "config_sha256": "right",
                "record": {
                    "ordinal": 0, "index": 101, "status": "SURVIVOR",
                },
            }))
            loaded = runner.load_checkpoints(Path(directory), "right", selected)
            self.assertEqual(set(loaded), {0})
            with self.assertRaisesRegex(ValueError, "incompatible"):
                runner.load_checkpoints(Path(directory), "wrong", selected)


if __name__ == "__main__":
    unittest.main()
