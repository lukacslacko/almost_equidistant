#!/usr/bin/env python3
"""Extract the exact 17,764-graph K7 rank-survivor corpus.

The full K7 JSON contains all 113,136 graphs entering the frozen rank layer.
This utility joins it to the independently archived Python decision TSV and
writes only the rows labelled ``SURVIVOR``.  Both inputs are hash-pinned.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path


EXPECTED_FULL_SHA256 = (
    "7097ddd5f333326bd4a4ab1c30379b054182542ca4cd3820f3cc3bc4308ec6da"
)
EXPECTED_DECISIONS_SHA256 = (
    "d331016042c14ba412a29e42b1f2ee06ee101a7d4376067b64e99e249677f10f"
)
EXPECTED_FULL_GRAPHS = 113_136
EXPECTED_SURVIVORS = 17_764


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full", type=Path, default=Path("/tmp/d6_k7_rank_full_residue.json")
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=Path("d6_k7_rank_python_full_decisions.tsv.gz"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    args = parser.parse_args()

    full_hash = sha256(args.full)
    decisions_hash = sha256(args.decisions)
    if full_hash != EXPECTED_FULL_SHA256:
        raise SystemExit(f"unexpected full K7 JSON SHA-256 {full_hash}")
    if decisions_hash != EXPECTED_DECISIONS_SHA256:
        raise SystemExit(f"unexpected decision archive SHA-256 {decisions_hash}")

    survivor_indices: set[int] = set()
    decision_rows = 0
    with gzip.open(args.decisions, "rt", encoding="ascii", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames is None or not {"index", "decision"} <= set(
            reader.fieldnames
        ):
            raise SystemExit("decision archive has an unexpected header")
        for row in reader:
            decision_rows += 1
            index = int(row["index"])
            decision = row["decision"]
            if decision not in {"REJECTED", "SURVIVOR"}:
                raise SystemExit(f"unexpected decision {decision!r} at {index}")
            if decision == "SURVIVOR":
                if index in survivor_indices:
                    raise SystemExit(f"duplicate survivor index {index}")
                survivor_indices.add(index)
    if (
        decision_rows != EXPECTED_FULL_GRAPHS
        or len(survivor_indices) != EXPECTED_SURVIVORS
    ):
        raise SystemExit(
            f"unexpected decision counts {decision_rows}/{len(survivor_indices)}"
        )

    full = json.loads(args.full.read_text(encoding="utf-8"))
    graphs = full.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != EXPECTED_FULL_GRAPHS:
        raise SystemExit("full K7 JSON has an unexpected graph population")
    selected = [
        graph for graph in graphs if int(graph["index"]) in survivor_indices
    ]
    selected_indices = {int(graph["index"]) for graph in selected}
    if selected_indices != survivor_indices or len(selected) != EXPECTED_SURVIVORS:
        raise SystemExit("full K7 JSON and decision archive do not join exactly")

    output = {
        "schema": 1,
        "description": "Exact K7 enhanced-rank survivors for later exact layers.",
        "sources": {
            "full_K7_JSON": str(args.full),
            "full_K7_JSON_sha256": full_hash,
            "decision_archive": str(args.decisions),
            "decision_archive_sha256": decisions_hash,
            "decision_rows": decision_rows,
        },
        "graphs": selected,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(output, stream, separators=(",", ":"), sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, args.output)
    print(
        f"wrote {len(selected)} K7 rank survivors to {args.output}; "
        f"SHA-256 {sha256(args.output)}"
    )


if __name__ == "__main__":
    main()
