#!/usr/bin/env python3
"""Build a tiny deterministic current-residue sample for tetrad screening."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rank-input", type=Path,
        default=Path("/tmp/d6_k7_rank_full_residue.json"),
    )
    parser.add_argument(
        "--selection", type=Path,
        default=Path("d6_k7_positive_dual_selection.json"),
    )
    parser.add_argument(
        "--decisions", type=Path,
        default=Path("d6_k7_positive_dual_full_decisions.tsv.gz"),
    )
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_rankone_tetrad_sample.json"),
    )
    args = parser.parse_args()
    if args.count <= 0:
        raise ValueError("sample count must be positive")
    rank_payload = json.loads(args.rank_input.read_text(encoding="utf-8"))
    graph_by_index = {
        int(graph["index"]): graph for graph in rank_payload["graphs"]
    }
    selection = json.loads(args.selection.read_text(encoding="utf-8"))
    selected = {int(value) for value in selection["selected_indices"]}
    candidates = []
    with gzip.open(args.decisions, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            if (
                row["status"] == "SURVIVOR"
                and index in selected
                and int(row["dual_passing_covers"]) == 1
                and int(row["strict_h_passing_cliques"]) == 0
            ):
                candidates.append((int(row["covers"]), index))
    candidates.sort()
    chosen = [index for _covers, index in candidates[:args.count]]
    if len(chosen) != args.count:
        raise ValueError("not enough qualifying current-residue graphs")
    if any(index not in graph_by_index for index in chosen):
        raise ValueError("selected index absent from rank input")
    output = {
        "schema": 1,
        "kind": "d6_k7_rankone_tetrad_current_residue_sample",
        "selection_rule": (
            "first count records sorted by (all-cover count,index) among "
            "full degree-one SURVIVOR rows with exactly one dual-passing "
            "cover and zero strict-H-passing saturating cliques"
        ),
        "count": args.count,
        "indices": chosen,
        "source": {
            "builder": {
                "path": Path(__file__).name,
                "sha256": sha256(Path(__file__)),
            },
            "rank_input": {
                "path": str(args.rank_input), "sha256": sha256(args.rank_input),
            },
            "positive_dual_selection": {
                "path": str(args.selection), "sha256": sha256(args.selection),
            },
            "positive_dual_decisions": {
                "path": str(args.decisions), "sha256": sha256(args.decisions),
            },
        },
        "graphs": [graph_by_index[index] for index in chosen],
    }
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
