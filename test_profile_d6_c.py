#!/usr/bin/env python3
"""Run exact positive and negative controls through the C profiler itself.

The production build has ``N=19``.  ``profile_d6.c`` deliberately permits a
compile-time ``N`` override so the same implementation can also consume the
smaller, known-realizable controls without padding them with fictitious
vertices.  No filter logic is replaced by this harness.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_d6_reference_filters import (
    cardinality_cover_failure,
    ineligible_disjoint_edge_failure,
    k7_hall_failure,
    lower_bound_18_graph,
    reflected_facet_graph,
    tight_cover_matching_failure,
)


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "profile_d6.c"


class CProfilerControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary = tempfile.TemporaryDirectory(prefix="profile-d6-controls-")
        cls.tmp = Path(cls._temporary.name)
        cls.binaries: dict[int, Path] = {}

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    @classmethod
    def profiler(cls, n: int) -> Path:
        if n not in cls.binaries:
            binary = cls.tmp / f"profile_d6_n{n}"
            subprocess.run(
                [
                    "cc",
                    "-O3",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-pthread",
                    f"-DN={n}",
                    "-DEXPECTED=1u",
                    "-o",
                    str(binary),
                    str(SOURCE),
                ],
                check=True,
            )
            cls.binaries[n] = binary
        return cls.binaries[n]

    @classmethod
    def profile(cls, adjacency: list[int]) -> dict:
        n = len(adjacency)
        corpus = cls.tmp / f"control_n{n}.txt"
        corpus.write_text(
            f"{n} " + " ".join(map(str, adjacency)) + "\n",
            encoding="ascii",
        )
        process = subprocess.run(
            [str(cls.profiler(n)), str(corpus), "-", "1", "1"],
            check=True,
            text=True,
            capture_output=True,
        )
        return json.loads(process.stdout)["populations"]["all"]

    def test_known_realizable_controls_pass(self) -> None:
        controls = [
            reflected_facet_graph(1),
            reflected_facet_graph(2),
            lower_bound_18_graph(),
        ]
        for adjacency in controls:
            with self.subTest(n=len(adjacency)):
                result = self.profile(adjacency)
                self.assertEqual(result["cumulative_exact_rejections"], 0)

    def test_ineligible_disjoint_edge_rejected(self) -> None:
        result = self.profile(ineligible_disjoint_edge_failure())
        self.assertEqual(
            result["K7_disjoint_edge_bounded_cover_rejections"], 1
        )
        self.assertEqual(result["cumulative_previous_exact_rejections"], 0)

    def test_cover_bound_seven_is_sharp_control(self) -> None:
        result = self.profile(cardinality_cover_failure())
        self.assertEqual(
            result["K7_disjoint_edge_bounded_cover_rejections"], 1
        )
        self.assertEqual(result["K7_clique_Hall_rejections"], 0)
        self.assertEqual(result["cumulative_previous_exact_rejections"], 0)

    def test_k7_hall_negative_control(self) -> None:
        result = self.profile(k7_hall_failure())
        self.assertEqual(result["K7_clique_Hall_rejections"], 1)

    def test_tight_cover_matching_negative_control(self) -> None:
        result = self.profile(tight_cover_matching_failure())
        self.assertEqual(result["K7_tight_cover_matching_rejections"], 1)

    def test_malformed_kill_log_is_rejected(self) -> None:
        adjacency = reflected_facet_graph(1)
        n = len(adjacency)
        corpus = self.tmp / "malformed_log_control.txt"
        kill_log = self.tmp / "malformed_kill_log.txt"
        corpus.write_text(
            f"{n} " + " ".join(map(str, adjacency)) + "\n",
            encoding="ascii",
        )
        kill_log.write_text("0\n", encoding="ascii")
        process = subprocess.run(
            [str(self.profiler(n)), str(corpus), str(kill_log), "1", "1"],
            text=True,
            capture_output=True,
        )
        self.assertEqual(process.returncode, 2)
        self.assertIn("malformed kill-log line 1", process.stderr)

    def test_trailing_corpus_record_is_rejected(self) -> None:
        adjacency = reflected_facet_graph(1)
        n = len(adjacency)
        record = f"{n} " + " ".join(map(str, adjacency)) + "\n"
        corpus = self.tmp / "trailing_corpus_control.txt"
        corpus.write_text(record + record, encoding="ascii")
        process = subprocess.run(
            [str(self.profiler(n)), str(corpus), "-"],
            text=True,
            capture_output=True,
        )
        self.assertEqual(process.returncode, 2)
        self.assertIn("trailing data after candidate 0", process.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
