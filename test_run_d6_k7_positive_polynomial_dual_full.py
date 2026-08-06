#!/usr/bin/env python3
"""Focused controls for the production K7 positive-dual wrapper."""

from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_k7_positive_dual_selection as builder
import run_d6_k7_positive_polynomial_dual_full as runner
import verify_d6_k7_positive_polynomial_dual as exact
import verify_d6_k7_positive_polynomial_dual_full as archive_verifier


ROOT = Path(__file__).resolve().parent


class PositiveDualFullRunnerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.input_path = ROOT / ".runs/d6_k7_rank_survivors.json"
        cls.selection_path = ROOT / "d6_k7_positive_dual_selection.json"

    def test_current_residue_rebuilds_exactly(self) -> None:
        recorded = json.loads(self.selection_path.read_text(encoding="utf-8"))
        rebuilt = builder.build_selection(self.input_path)
        self.assertEqual(recorded, rebuilt)
        self.assertEqual(recorded["selected"], 12_941)
        self.assertEqual(
            recorded["selected_indices_sha256"],
            "a07bafe688c45ff132921519320a19687d2917a79a6c9009bda51a19a6b7f290",
        )
        increments = {
            item["layer"]: item["incremental_rejections"]
            for item in recorded["layer_accounting"]
        }
        self.assertEqual(increments, builder.EXPECTED_INCREMENTAL_COUNTS)

    def test_load_selection_checks_explicit_manifest_pin(self) -> None:
        digest = runner.sha256(self.selection_path)
        graphs, provenance = runner.load_selection(
            self.input_path, self.selection_path, digest
        )
        self.assertEqual(len(graphs), 12_941)
        self.assertEqual(provenance["selection_report_sha256"], digest)
        with self.assertRaises(ValueError):
            runner.load_selection(self.input_path, self.selection_path, "0" * 64)

    def test_atomic_gzip_is_byte_deterministic(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT / ".runs") as directory:
            path = Path(directory) / "archive.gz"
            runner.atomic_deterministic_gzip(path, b"proof\narchive\n")
            first = path.read_bytes()
            runner.atomic_deterministic_gzip(path, b"proof\narchive\n")
            second = path.read_bytes()
            self.assertEqual(first, second)
            self.assertEqual(gzip.decompress(second), b"proof\narchive\n")

    def test_complete_checkpoint_resumes_without_graph_work(self) -> None:
        graph = dict(json.loads(self.input_path.read_text(encoding="utf-8"))[
            "graphs"
        ][0])
        graph["ordinal"] = 0
        result = {
            "ordinal": 0,
            "index": int(graph["index"]),
            "status": "COMPLETE",
            "analysis": {
                "index": int(graph["index"]),
                "applicable": True,
                "seeds": 1,
                "covers": 0,
                "enhanced_passing_covers": 0,
                "strict_h_passing_covers": 0,
                "dual_passing_covers": 0,
                "strict_h_passing_cliques": 0,
                "dual_failing_cliques": 0,
                "strict_h_rejected": False,
                "dual_rejected": False,
                "marginal_dual_rejected": False,
                "dual_failure_witnesses": [],
            },
            "elapsed_seconds": 0.125,
            "error_type": None,
            "error_message": None,
            "error_traceback": None,
        }
        configuration = {"schema": 1, "multiplier_degree": 1, "test": True}
        config_hash = runner.stable_hash(configuration)
        with tempfile.TemporaryDirectory(dir=ROOT / ".runs") as directory:
            base = Path(directory)
            checkpoint_dir = base / "checkpoints"
            checkpoint = runner.result_path(checkpoint_dir, graph)
            runner.atomic_json(
                checkpoint, runner.checkpoint_value(config_hash, result)
            )
            arguments = {
                "graphs": [graph],
                "configuration": configuration,
                "workers": 11,
                "max_inflight": 22,
                "retry_infra_errors": False,
                "checkpoint_dir": checkpoint_dir,
                "decisions_archive": base / "decisions.tsv.gz",
                "certificate_archive": base / "certificates.jsonl.gz",
                "checkpoint_copy": base / "checkpoint.json",
                "report_path": base / "report.json",
                "pid_file": base / "pid.json",
            }
            report, has_infra = runner.run_campaign(**arguments)
            self.assertFalse(has_infra)
            self.assertEqual(report["runtime"]["newly_completed_graphs"], 0)
            decisions_hash = hashlib.sha256(
                arguments["decisions_archive"].read_bytes()
            ).hexdigest()
            report2, has_infra2 = runner.run_campaign(**arguments)
            self.assertFalse(has_infra2)
            self.assertEqual(report2["runtime"]["newly_completed_graphs"], 0)
            self.assertEqual(
                decisions_hash,
                hashlib.sha256(arguments["decisions_archive"].read_bytes()).hexdigest(),
            )

    def test_new_campaign_refuses_preexisting_output(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT / ".runs") as directory:
            base = Path(directory)
            output = base / "decisions.tsv.gz"
            output.write_bytes(b"do not overwrite")
            with self.assertRaisesRegex(RuntimeError, "refusing to overwrite"):
                runner.run_campaign(
                    [], configuration={"schema": 1, "multiplier_degree": 1},
                    workers=11, max_inflight=22, retry_infra_errors=False,
                    checkpoint_dir=base / "checkpoints",
                    decisions_archive=output,
                    certificate_archive=base / "certificates.jsonl.gz",
                    checkpoint_copy=base / "checkpoint.json",
                    report_path=base / "report.json",
                    pid_file=base / "pid.json",
                )

    def test_known_rejection_round_trips_through_compact_archives(self) -> None:
        graphs = json.loads(self.input_path.read_text(encoding="utf-8"))["graphs"]
        graph = dict(next(
            item for item in graphs if int(item["index"]) == 3_649_646
        ))
        graph["ordinal"] = 0
        result = runner.evaluate_graph(graph, 1)
        runner.validate_result(result, graph)
        self.assertEqual(result["status"], "COMPLETE")
        self.assertTrue(result["analysis"]["marginal_dual_rejected"])

        rows = archive_verifier.parse_decisions(runner.render_decisions([result]))
        certificates = archive_verifier.parse_certificates(
            runner.render_certificates([result])
        )
        witnesses = certificates[graph["index"]]["dual_failure_witnesses"]
        rows[0]["analysis"]["dual_failure_witnesses"] = witnesses
        for witness in witnesses:
            archive_verifier.verify_witness(graph, witness)
        exact.verify_graph(graph, rows[0]["analysis"])


if __name__ == "__main__":
    unittest.main()
