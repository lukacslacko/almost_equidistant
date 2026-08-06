#!/usr/bin/env python3
"""Build the hash-pinned 12,839-graph input for the full K7 tetrad pass."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


EXPECTED_RANK_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_RANK_GRAPHS = 17_764
EXPECTED_PRIOR_SELECTION_SHA256 = (
    "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0"
)
EXPECTED_PRIOR_REPORT_SHA256 = (
    "c2de7e06bbe4f3d163a4f747d42721f7f60da985b867ce74669664431738e345"
)
EXPECTED_PRIOR_DECISIONS_SHA256 = (
    "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
)
EXPECTED_PRIOR_CERTIFICATES_SHA256 = (
    "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
)
EXPECTED_PRIOR_GRAPHS = 12_941
EXPECTED_SURVIVORS = 12_839


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def require_hash(path: Path, expected: str) -> str:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(f"hash mismatch for {path}: {observed} != {expected}")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rank-input", type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    parser.add_argument(
        "--prior-selection", type=Path,
        default=Path("d6_k7_positive_dual_selection.json"),
    )
    parser.add_argument(
        "--prior-report", type=Path,
        default=Path("d6_k7_positive_dual_full_report.json"),
    )
    parser.add_argument(
        "--prior-decisions", type=Path,
        default=Path("d6_k7_positive_dual_full_decisions.tsv.gz"),
    )
    parser.add_argument(
        "--prior-certificates", type=Path,
        default=Path("d6_k7_positive_dual_full_certificates.jsonl.gz"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_selection.json"),
    )
    args = parser.parse_args()
    hashes = {
        "rank_input": require_hash(
            args.rank_input, EXPECTED_RANK_INPUT_SHA256
        ),
        "prior_selection": require_hash(
            args.prior_selection, EXPECTED_PRIOR_SELECTION_SHA256
        ),
        "prior_report": require_hash(
            args.prior_report, EXPECTED_PRIOR_REPORT_SHA256
        ),
        "prior_decisions": require_hash(
            args.prior_decisions, EXPECTED_PRIOR_DECISIONS_SHA256
        ),
        "prior_certificates": require_hash(
            args.prior_certificates, EXPECTED_PRIOR_CERTIFICATES_SHA256
        ),
    }
    rank_payload = json.loads(args.rank_input.read_text(encoding="utf-8"))
    if len(rank_payload.get("graphs", [])) != EXPECTED_RANK_GRAPHS:
        raise ValueError("rank input graph count mismatch")
    rank_indices = {int(graph["index"]) for graph in rank_payload["graphs"]}
    if len(rank_indices) != EXPECTED_RANK_GRAPHS:
        raise ValueError("rank input repeats an index")
    prior_selection = json.loads(
        args.prior_selection.read_text(encoding="utf-8")
    )
    prior_indices = [int(value) for value in prior_selection["selected_indices"]]
    if len(prior_indices) != EXPECTED_PRIOR_GRAPHS:
        raise ValueError("prior selection count mismatch")
    prior_report = json.loads(args.prior_report.read_text(encoding="utf-8"))
    summary = prior_report["summary"]
    if not (
        summary["graphs"] == EXPECTED_PRIOR_GRAPHS
        and summary["complete"] == EXPECTED_PRIOR_GRAPHS
        and summary["infra_errors"] == 0
        and summary["survivors"] == EXPECTED_SURVIVORS
        and prior_report["artifacts"]["decisions_archive_sha256"]
        == hashes["prior_decisions"]
        and prior_report["artifacts"]["certificate_archive_sha256"]
        == hashes["prior_certificates"]
    ):
        raise ValueError("prior full report is incomplete or inconsistent")
    rows = []
    with gzip.open(
        args.prior_decisions, "rt", encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if [int(row["index"]) for row in rows] != prior_indices:
        raise ValueError("prior decision order differs from prior selection")
    if any(row["status"] == "INFRA_ERROR" for row in rows):
        raise ValueError("prior decisions contain an infrastructure error")
    survivors = [
        int(row["index"]) for row in rows if row["status"] == "SURVIVOR"
    ]
    if len(survivors) != EXPECTED_SURVIVORS:
        raise ValueError("prior survivor count mismatch")
    if any(index not in rank_indices for index in survivors):
        raise ValueError("prior survivor absent from rank input")
    source = Path(__file__)
    output = {
        "schema": 1,
        "kind": "d6_k7_rankone_tetrad_full_selection",
        "rank_input": {
            "path": str(args.rank_input),
            "sha256": hashes["rank_input"],
            "graphs": EXPECTED_RANK_GRAPHS,
        },
        "prior_degree_one": {
            "selection": {
                "path": str(args.prior_selection),
                "sha256": hashes["prior_selection"],
            },
            "report": {
                "path": str(args.prior_report),
                "sha256": hashes["prior_report"],
            },
            "decisions": {
                "path": str(args.prior_decisions),
                "sha256": hashes["prior_decisions"],
            },
            "certificates": {
                "path": str(args.prior_certificates),
                "sha256": hashes["prior_certificates"],
            },
            "graphs": EXPECTED_PRIOR_GRAPHS,
            "rejected": EXPECTED_PRIOR_GRAPHS - EXPECTED_SURVIVORS,
            "survivors": EXPECTED_SURVIVORS,
        },
        "selected_indices": survivors,
        "selected_indices_sha256": stable_hash(survivors),
        "source": {"path": source.name, "sha256": sha256(source)},
    }
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
