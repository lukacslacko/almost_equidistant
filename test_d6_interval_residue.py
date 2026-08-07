#!/usr/bin/env python3
"""Focused controls for the restartable dimension-six interval runner."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_d6_interval_residue as runner


def graph(ordinal: int = 0, index: int = 17) -> dict:
    return {
        "ordinal": ordinal,
        "index": index,
        "population": "K7",
        "adjacency": tuple(0 for _ in range(runner.N)),
    }


class QuantifierSemantics(unittest.TestCase):
    def test_killed_requires_every_slice_of_one_order(self) -> None:
        calls = iter((
            ("SURVIVORS", 2, 1),
            ("KILLED", 3, 0),
            ("KILLED", 5, 0),
        ))
        orders = [((0,), [1]), ((0,), [2])]
        with (
            patch.object(runner.cdriver6, "gen_orders", return_value=orders),
            patch.object(runner.cdriver6, "decide6", side_effect=lambda *a, **k: next(calls)),
            patch.object(runner, "ncircle", return_value=0),
        ):
            result = runner.analyze_graph(graph(), 12, 100, True, True, 2)
        self.assertEqual(result.status, "KILLED")
        self.assertEqual(result.winning_order, 1)
        self.assertEqual(result.kernel_killed, 2)
        self.assertEqual(result.kernel_survivors, 1)
        self.assertEqual(result.kernel_aborts, 0)
        self.assertEqual(len(result.winning_slice_records or ()), 2)

    def test_abort_and_unresolved_are_distinct(self) -> None:
        orders = [((0,), [1]), ((0,), [2])]
        calls = iter((("SURVIVORS", 2, 1), ("ABORT", 7, 3)))
        with (
            patch.object(runner.cdriver6, "gen_orders", return_value=orders),
            patch.object(runner.cdriver6, "decide6", side_effect=lambda *a, **k: next(calls)),
            patch.object(runner, "ncircle", return_value=0),
        ):
            result = runner.analyze_graph(graph(), 12, 100, True, True, 1)
        self.assertEqual(result.status, "ABORT")
        self.assertEqual(result.kernel_survivors, 1)
        self.assertEqual(result.kernel_aborts, 1)

        with patch.object(runner.cdriver6, "gen_orders", return_value=[]):
            unresolved = runner.analyze_graph(graph(), 12, 100, True, True, 1)
        self.assertEqual(unresolved.status, "UNRESOLVED")
        self.assertEqual(unresolved.kernel_calls, 0)

    def test_exception_becomes_infrastructure_error(self) -> None:
        with patch.object(runner.cdriver6, "gen_orders", side_effect=RuntimeError("boom")):
            result = runner.analyze_graph_safe(
                graph(),
                {
                    "orders": 12,
                    "cap": 100,
                    "zero_circle_only": True,
                    "include_bulk_order": True,
                    "slices": 1,
                },
            )
        self.assertEqual(result.status, "INFRA_ERROR")
        self.assertEqual(result.error_type, "RuntimeError")
        self.assertIn("boom", result.error_message or "")


class CheckpointControls(unittest.TestCase):
    def test_complete_results_resume_without_recomputation(self) -> None:
        graphs = [graph(0, 17), graph(1, 23)]
        configuration = {
            "schema": 2,
            "search": {
                "orders": 1,
                "cap": 10,
                "zero_circle_only": True,
                "include_bulk_order": True,
                "slices": 1,
            },
        }

        def fake(item: dict, _parameters: dict) -> runner.Result:
            status = "KILLED" if item["ordinal"] == 0 else "ABORT"
            result = runner.Result(
                item["ordinal"], item["index"], item["population"], status,
                0 if status == "KILLED" else -1, 0 if status == "KILLED" else -1,
                1, 1, 1, int(status == "KILLED"), 0, int(status == "ABORT"),
                4, 0, 0.01,
            )
            if status != "KILLED":
                return result
            return runner.Result(
                **{
                    **result.__dict__,
                    "winning_seed": tuple(range(7)),
                    "winning_placement_order": tuple(range(7, runner.N)),
                    "winning_slice_records": ({
                        "part": 0,
                        "lo": 0.0,
                        "hi": 2.0 * math.pi,
                        "status": "KILLED",
                        "nodes": 4,
                        "unresolved_cells": 0,
                    },),
                }
            )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(runner, "analyze_graph_safe", side_effect=fake),
                patch.object(runner, "initialize_native_thread", return_value=None),
            ):
                first, first_infra = runner.run_campaign(
                    graphs,
                    configuration=configuration,
                    workers=2,
                    checkpoint_every=1,
                    progress_every=99,
                    retry_infra_errors=False,
                    checkpoint_dir=root / "checkpoints",
                    report_path=root / "first.json",
                    decisions_path=root / "first.tsv",
                )
            self.assertFalse(first_infra)
            self.assertEqual(first["status_counts"]["KILLED"], 1)
            self.assertEqual(first["status_counts"]["ABORT"], 1)

            with (
                patch.object(
                    runner, "analyze_graph_safe",
                    side_effect=AssertionError("resume recomputed a graph"),
                ) as analyze,
                patch.object(runner, "initialize_native_thread", return_value=None),
            ):
                resumed, resumed_infra = runner.run_campaign(
                    graphs,
                    configuration=configuration,
                    workers=2,
                    checkpoint_every=1,
                    progress_every=99,
                    retry_infra_errors=False,
                    checkpoint_dir=root / "checkpoints",
                    report_path=root / "resumed.json",
                    decisions_path=root / "resumed.tsv",
                )
            analyze.assert_not_called()
            self.assertFalse(resumed_infra)
            self.assertEqual(resumed["runtime"]["newly_completed_graphs"], 0)
            self.assertEqual(
                (root / "first.tsv").read_text(),
                (root / "resumed.tsv").read_text(),
            )

    def test_unapproved_optional_report_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fake.json"
            path.write_text(json.dumps({
                "decisions": {
                    "strict_H_rejected": [1, 2],
                    "strict_H_survivors": [],
                }
            }))
            with self.assertRaisesRegex(ValueError, "unapproved strict-H"):
                runner.apply_k7_strict_h(
                    {1: "K7", 2: "K7"}, {1, 2}, path
                )


if __name__ == "__main__":
    unittest.main()
