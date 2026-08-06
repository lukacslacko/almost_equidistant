#!/usr/bin/env python3
"""Extract all 1,106 survivors of the pure K6 two-light-ray filter.

The exact indices come from the full profiler's decision TSV.  This extractor
hashes that entire sidecar and the complete 3,971,787-graph corpus, verifies
the committed aggregate profile, and writes the original adjacency rows for
all K6-only graphs with ``K6_two_light_ray=0``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from d6_k6_lorentz_reference import find_clique_mask, validate_graph


EXPECTED_GRAPHS = 3_971_787
EXPECTED_RESIDUE = 1_106
EXPECTED_CORPUS_SHA256 = (
    "12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5"
)
EXPECTED_DECISIONS_SHA256 = (
    "c2e39f7b366a67c2ea117531e2cd5788a0f40eca32e605601bec59152aecea75"
)
EXPECTED_PROFILE_SHA256 = (
    "77d507fcc55c4e34a44af2de1ce7c16d0743476ea09d08ce0215e65c80ab17d9"
)
EXPECTED_HEADER = (
    "index\tomega\tlink_mask\treflection\tdefect_csp\tK7_Hall\tK7_cover\t"
    "K7_tight_cover\tK6_Hall\tK6_two_light_ray\tK7_Hall_seed\tK7_Hall_aux\t"
    "K7_cover_seed\tK7_cover_aux\tK7_tight_cover_seed\t"
    "K7_tight_cover_aux\tK6_Hall_seed\tK6_Hall_aux\t"
    "K6_two_light_ray_seed\tK6_two_light_ray_aux"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def scan_decisions(path: Path) -> tuple[str, list[int], int]:
    """Hash the full sidecar and return the exact pure-filter survivors."""

    digest = hashlib.sha256()
    survivors = []
    rows = 0
    with path.open("rb") as stream:
        raw_header = stream.readline()
        digest.update(raw_header)
        header = raw_header.decode("ascii").rstrip("\r\n")
        if header != EXPECTED_HEADER:
            raise ValueError("unexpected decisions TSV header")
        for expected_index, raw in enumerate(stream):
            digest.update(raw)
            fields = raw.split()
            if len(fields) != 20:
                raise ValueError(f"malformed decisions row {expected_index}")
            index = int(fields[0])
            if index != expected_index:
                raise ValueError(
                    f"decision index {index}, expected {expected_index}"
                )
            omega = int(fields[1])
            k6_hall = int(fields[8])
            pure_lorentz = int(fields[9])
            if omega == 6 and pure_lorentz == 0:
                if k6_hall:
                    raise ValueError(f"K6 Hall rejected survivor {index}")
                survivors.append(index)
            rows += 1
    return digest.hexdigest(), survivors, rows


def scan_corpus(
    path: Path, wanted: set[int]
) -> tuple[str, dict[int, list[int]], int]:
    """Hash every corpus byte and parse only the 1,106 wanted records."""

    digest = hashlib.sha256()
    selected = {}
    rows = 0
    with path.open("rb") as stream:
        for index, raw in enumerate(stream):
            digest.update(raw)
            if index in wanted:
                fields = raw.split()
                if len(fields) != 20 or fields[0] != b"19":
                    raise ValueError(f"malformed corpus row {index}")
                adj = [int(field) for field in fields[1:]]
                validate_graph(adj, require_alpha_two=True)
                if find_clique_mask(adj, 7):
                    raise ValueError(f"residue row {index} contains K7")
                if not find_clique_mask(adj, 6):
                    raise ValueError(f"residue row {index} has no K6")
                selected[index] = adj
            rows = index + 1
    return digest.hexdigest(), selected, rows


def check_profile(path: Path) -> str:
    digest = sha256(path)
    if digest != EXPECTED_PROFILE_SHA256:
        raise ValueError(f"unexpected profile SHA-256 {digest}")
    with path.open(encoding="utf-8") as stream:
        profile = json.load(stream)
    deferred = profile["populations"]["deferred"]
    if deferred["K6_odd_component_two_light_ray_rejections"] != 174_713:
        raise ValueError("unexpected pure K6 Lorentz rejection count")
    if deferred["residue_by_clique_number"]["6"] != EXPECTED_RESIDUE:
        raise ValueError("unexpected K6 residue count")
    return digest


def build_residue(corpus: Path, decisions: Path, profile: Path) -> dict:
    decision_hash, indices, decision_rows = scan_decisions(decisions)
    if decision_hash != EXPECTED_DECISIONS_SHA256:
        raise ValueError(f"unexpected decisions SHA-256 {decision_hash}")
    if decision_rows != EXPECTED_GRAPHS:
        raise ValueError(f"decisions have {decision_rows} graph rows")
    if len(indices) != EXPECTED_RESIDUE:
        raise ValueError(f"decisions contain {len(indices)} K6 survivors")
    corpus_hash, selected, corpus_rows = scan_corpus(corpus, set(indices))
    if corpus_hash != EXPECTED_CORPUS_SHA256:
        raise ValueError(f"unexpected corpus SHA-256 {corpus_hash}")
    if corpus_rows != EXPECTED_GRAPHS:
        raise ValueError(f"corpus has {corpus_rows} graph rows")
    if len(selected) != EXPECTED_RESIDUE:
        raise ValueError(f"loaded only {len(selected)} residue graphs")
    profile_hash = check_profile(profile)

    return {
        "schema": 1,
        "description": (
            "Complete 1,106-graph K6-only residue after the exact pure "
            "odd-component/two-light-ray filter."
        ),
        "sources": {
            "corpus": {
                "file": corpus.name,
                "sha256": corpus_hash,
                "graphs": corpus_rows,
                "tracked": False,
            },
            "pure_filter_decisions": {
                "file": str(decisions),
                "sha256": decision_hash,
                "graph_rows": decision_rows,
                "tracked": False,
            },
            "pure_filter_profile": {
                "file": profile.name,
                "sha256": profile_hash,
                "K6_rejections": 174_713,
                "K6_residue": EXPECTED_RESIDUE,
            },
        },
        "selection": {
            "rule": "omega == 6 and K6_two_light_ray == 0",
            "complete_not_sampled": True,
            "indices_in_increasing_corpus_order": True,
        },
        "preconditions": {
            "alpha_at_most_two_checked": True,
            "clique_number_exactly_six_checked": True,
            "candidate_nonedges": "unconstrained; allowed support may be unused",
        },
        "graphs": [
            {"index": index, "adjacency": selected[index]} for index in indices
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("aeq_d6_n19.txt"))
    parser.add_argument(
        "--decisions",
        type=Path,
        default=Path("/tmp/d6_profile_stream_decisions.tsv"),
    )
    parser.add_argument("--profile", type=Path, default=Path("d6_profile.json"))
    parser.add_argument(
        "--output", type=Path, default=Path("d6_k6_support_residue.json")
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    residue = build_residue(args.corpus, args.decisions, args.profile)
    rendered = json.dumps(residue, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} is not reproducible")
        print(f"PASS: {args.output} is reproducible")
    else:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}: {len(residue['graphs'])} complete survivors")


if __name__ == "__main__":
    main()
