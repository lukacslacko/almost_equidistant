#!/usr/bin/env python3
"""Compare C-profiler decisions with the independent Python reference sample.

Build the current profiler first, for example:

    cc -O3 -Wall -Wextra -Werror -pthread -o /tmp/profile_d6 profile_d6.c
    python3 crosscheck_d6_reference_c.py /tmp/profile_d6 \
        d6_reference_sample.json

The profiler's optional one-graph limit makes this intentionally simple: each
sample graph is run in its own process, so aggregate agreement cannot conceal
two decisions swapped between graphs.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

from d6_reference_filters import evaluate


C_KEYS = {
    "K7_clique_Hall": "K7_clique_Hall_rejections",
    "K7_disjoint_edge_bounded_cover": (
        "K7_disjoint_edge_bounded_cover_rejections"
    ),
    "K7_tight_cover_matching": "K7_tight_cover_matching_rejections",
    "K6_clique_Hall": "K6_clique_Hall_rejections",
}


def c_decisions(profiler: Path, adj: list[int], directory: Path) -> tuple[dict, int]:
    """Run the C profiler on one graph and return decisions and clique number."""

    corpus = directory / "one_graph.txt"
    kill_log = directory / "empty_kill_log.txt"
    corpus.write_text("19 " + " ".join(map(str, adj)) + "\n", encoding="ascii")
    kill_log.write_text("", encoding="ascii")
    process = subprocess.run(
        [str(profiler), str(corpus), str(kill_log), "1", "1"],
        check=True,
        text=True,
        capture_output=True,
    )
    population = json.loads(process.stdout)["populations"]["all"]
    decisions = {
        name: bool(population[key]) for name, key in C_KEYS.items()
    }
    clique_counts = population["clique_number"]
    clique_number = 7 if clique_counts["7"] else 6 if clique_counts["6"] else 5
    return decisions, clique_number


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiler", type=Path)
    parser.add_argument("sample", type=Path)
    args = parser.parse_args()
    profiler = args.profiler.resolve()
    if not profiler.is_file():
        raise SystemExit(f"profiler does not exist: {profiler}")
    with args.sample.open(encoding="utf-8") as stream:
        sample = json.load(stream)

    per_stratum: dict[str, dict[str, int]] = {}
    with tempfile.TemporaryDirectory(prefix="d6-reference-crosscheck-") as tmp:
        directory = Path(tmp)
        for graph in sample["graphs"]:
            python_result = evaluate(graph["adjacency"])
            expected = {
                name: decision.rejected for name, decision in python_result.items()
            }
            observed, omega = c_decisions(profiler, graph["adjacency"], directory)
            expected_omega = 7 if graph["clique_class"] == "K7" else 6
            if observed != expected or omega != expected_omega:
                raise SystemExit(
                    "C/Python mismatch at original corpus index "
                    f"{graph['index']}: expected decisions {expected}, "
                    f"observed {observed}, expected omega {expected_omega}, "
                    f"observed {omega}"
                )
            counts = per_stratum.setdefault(
                graph["stratum"], {name: 0 for name in C_KEYS}
            )
            for name, rejected in observed.items():
                counts[name] += int(rejected)

    print(
        json.dumps(
            {
                "result": "PASS",
                "graphs_compared_individually": len(sample["graphs"]),
                "rejections_by_stratum": per_stratum,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
