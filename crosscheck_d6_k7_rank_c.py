#!/usr/bin/env python3
"""Compile and cross-check the standalone C K7 evaluator on all 512 samples.

The Python report is the independent decision oracle for both the original
support/subspace-K/component-B union and the enhanced clique, complement-
degree, rational basis-kernel, and saturating-mask rules.  The C evaluator
checks only the now-proved cap-three covers; the Python report reaches the
same endpoint by retaining the larger covers and applying the proved direct
cap transitions.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "profile_d6_k7_rank.c"
SAMPLE = ROOT / "d6_k7_rank_sample.json"
REFERENCE_REPORT = ROOT / "d6_k7_rank_report.json"


def selected(rows: list[dict[str, str]], column: str) -> set[int]:
    return {int(row["index"]) for row in rows if int(row[column])}


def main() -> None:
    sample = json.loads(SAMPLE.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE_REPORT.read_text(encoding="utf-8"))
    workers = max(1, os.cpu_count() or 1)
    with tempfile.TemporaryDirectory(prefix="d6-k7-rank-c-") as directory:
        temporary = Path(directory)
        binary = temporary / "profile_d6_k7_rank"
        input_path = temporary / "sample.txt"
        decisions_path = temporary / "decisions.tsv"
        input_path.write_text(
            "".join(
                f"{graph['index']} 19 "
                f"{' '.join(map(str, graph['adjacency']))}\n"
                for graph in sample["graphs"]
            ),
            encoding="ascii",
        )
        subprocess.run(
            [
                "cc",
                "-O3",
                "-Wall",
                "-Wextra",
                "-Werror",
                "-pthread",
                "-o",
                str(binary),
                str(SOURCE),
            ],
            check=True,
            cwd=ROOT,
        )
        completed = subprocess.run(
            [
                str(binary),
                str(input_path),
                str(workers),
                str(decisions_path),
            ],
            check=True,
            cwd=ROOT,
            env={**os.environ, "D6_ALLOW_BACKGROUND_TEST_ONLY": "1"},
            text=True,
            capture_output=True,
        )
        c_report = json.loads(completed.stdout)
        with decisions_path.open(newline="", encoding="ascii") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))

    if len(rows) != len(sample["graphs"]):
        raise AssertionError(f"C emitted {len(rows)} decision rows")
    expected_order = [graph["index"] for graph in sample["graphs"]]
    observed_order = [int(row["index"]) for row in rows]
    if observed_order != expected_order:
        raise AssertionError("C decisions changed sample order")
    if any(int(row["internal_error"]) for row in rows):
        raise AssertionError("C exact-pattern invariant failed")

    expected = reference["individual_graph_decisions"]
    python_joint = set(expected["joint_existential_rejected"])
    python_cap3 = set(expected["cap4_to3_rejected"])
    python_survivors = set(expected["survivors"])
    c_base_joint = selected(rows, "base_joint_rejected")
    c_cap3 = selected(rows, "cap3_rejected")
    c_base_survivors = set(observed_order) - c_base_joint
    if c_base_joint != python_joint:
        raise AssertionError(
            "base joint mismatch: "
            f"C-only={sorted(c_base_joint-python_joint)[:12]}, "
            f"Python-only={sorted(python_joint-c_base_joint)[:12]}"
        )
    if c_cap3 != python_cap3:
        raise AssertionError(
            f"cap-3 mismatch: C={sorted(c_cap3)}, "
            f"Python={sorted(python_cap3)}"
        )
    if c_base_survivors != python_survivors:
        raise AssertionError("base survivor set differs from Python")

    expected_histogram = {
        key: int(value)
        for key, value in reference["cover_size_histogram"].items()
        if int(key) <= 3
    }
    if c_report["cover_size_histogram"] != expected_histogram:
        raise AssertionError(
            f"cover histogram mismatch: {c_report['cover_size_histogram']} "
            f"!= {expected_histogram}"
        )
    for key in ("graphs", "applicable_K7_graphs", "seeds_checked"):
        reference_key = "applicable_K7_graphs" if key == "applicable_K7_graphs" else key
        if c_report[key] != reference[reference_key]:
            raise AssertionError(f"{key} mismatch")
    expected_covers = sum(expected_histogram.values())
    if c_report["eligible_cap3_covers_checked"] != expected_covers:
        raise AssertionError("cap-3 cover count mismatch")
    expected_support_failures = int(
        reference["sequential_marginal_cover_failures"]["support_after_cap3"]
    )
    if c_report["cover_failures"]["support"] != expected_support_failures:
        raise AssertionError("support-CSP cover decisions differ")

    c_full_joint = selected(rows, "joint_rejected")
    if not c_base_joint <= c_full_joint:
        raise AssertionError("cheap exact rules removed a base rejection")
    graph_decisions = reference["individual_graph_decisions"]
    if "enhanced_joint_existential_rejected" in graph_decisions:
        if c_full_joint != set(
            graph_decisions["enhanced_joint_existential_rejected"]
        ):
            raise AssertionError("enhanced joint graph decisions differ")
    print(
        "K7 C/Python cross-check: PASS; "
        f"{len(rows)} graphs, {c_report['seeds_checked']} seeds, "
        f"{expected_covers} cap-3 covers; "
        f"base rejected {len(c_base_joint)}, "
        f"enhanced exact {len(c_full_joint)}; "
        f"C kernel {c_report['wall_seconds']:.3f}s/{workers} workers"
    )


if __name__ == "__main__":
    main()
