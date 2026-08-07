#!/usr/bin/env python3
"""Extract a deterministic sample of the pre-Lorentz K6-only residue.

The selection is the first 1,024 K7-free graphs in a fixed 32,768-index
SplitMix64 pool.  The extractor hashes the complete corpus and old kill log,
checks that every selected K6-only graph is deferred, and binds the preserved
input profile saying all 175,819 K6-only graphs reached this Lorentz filter.
Indices are zero based.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from d6_k6_lorentz_reference import find_clique_mask, validate_graph


EXPECTED_GRAPHS = 3_971_787
EXPECTED_K6_RESIDUE = 175_819
EXPECTED_CORPUS_SHA256 = (
    "12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5"
)
EXPECTED_KILL_LOG_SHA256 = (
    "11c790b12d3543f4476f1eb4c221203f806f1c434f06360b7d25346cace5e305"
)
EXPECTED_PROFILE_SHA256 = (
    "2229bb5cdea73cb7fd37669747af743e4f673f88aaf2620e37b1d5513e949fde"
)
SAMPLE_SIZE = 1_024
POOL_SIZE = 32_768
SPLITMIX64_SEED = 0x6C06_0DDC_0A57_2026
MASK64 = (1 << 64) - 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_pool() -> list[int]:
    """Return a platform-independent sequence of unique corpus indices."""

    state = SPLITMIX64_SEED
    answer = []
    seen = set()
    while len(answer) < POOL_SIZE:
        state = (state + 0x9E3779B97F4A7C15) & MASK64
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
        z ^= z >> 31
        index = z % EXPECTED_GRAPHS
        if index not in seen:
            seen.add(index)
            answer.append(index)
    return answer


def hash_and_load_candidates(
    path: Path, wanted: set[int]
) -> tuple[str, dict[int, list[int]], int]:
    """Hash every byte but parse adjacency only at selected pool indices."""

    digest = hashlib.sha256()
    selected: dict[int, list[int]] = {}
    count = 0
    with path.open("rb") as stream:
        for index, raw in enumerate(stream):
            digest.update(raw)
            if index in wanted:
                fields = raw.split()
                if len(fields) != 20 or fields[0] != b"19":
                    raise ValueError(f"malformed candidate line {index}")
                adj = [int(field) for field in fields[1:]]
                validate_graph(adj)
                selected[index] = adj
            count = index + 1
    return digest.hexdigest(), selected, count


def hash_and_load_killed(path: Path, wanted: set[int]) -> tuple[str, set[int]]:
    """Hash the full old kill log and retain membership for the pool."""

    digest = hashlib.sha256()
    killed = set()
    with path.open("rb") as stream:
        for line_number, raw in enumerate(stream, 1):
            digest.update(raw)
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(f"malformed kill-log line {line_number}")
            index = int(fields[0])
            if index in wanted:
                killed.add(index)
    return digest.hexdigest(), killed


def check_profile(path: Path) -> tuple[str, dict]:
    profile_hash = sha256(path)
    if profile_hash != EXPECTED_PROFILE_SHA256:
        raise ValueError(f"unexpected profile SHA-256 {profile_hash}")
    with path.open(encoding="utf-8") as stream:
        profile = json.load(stream)
    deferred = profile["populations"]["deferred"]
    if deferred["clique_number"]["6"] != EXPECTED_K6_RESIDUE:
        raise ValueError("profile K6-only deferred count changed")
    if deferred["residue_by_clique_number"]["6"] != EXPECTED_K6_RESIDUE:
        raise ValueError("profile K6-only exact residue count changed")
    if deferred["K6_clique_Hall_rejections"] != 0:
        raise ValueError("profile unexpectedly rejects a K6-only graph")
    return profile_hash, profile


def build_sample(corpus: Path, kill_log: Path, profile_path: Path) -> dict:
    pool = selection_pool()
    pool_set = set(pool)
    corpus_hash, loaded, corpus_count = hash_and_load_candidates(corpus, pool_set)
    kill_hash, killed = hash_and_load_killed(kill_log, pool_set)
    profile_hash, _ = check_profile(profile_path)
    if corpus_hash != EXPECTED_CORPUS_SHA256:
        raise ValueError(f"unexpected corpus SHA-256 {corpus_hash}")
    if kill_hash != EXPECTED_KILL_LOG_SHA256:
        raise ValueError(f"unexpected kill-log SHA-256 {kill_hash}")
    if corpus_count != EXPECTED_GRAPHS:
        raise ValueError(f"corpus has {corpus_count}, expected {EXPECTED_GRAPHS}")
    if len(loaded) != POOL_SIZE:
        raise ValueError(f"loaded {len(loaded)} of {POOL_SIZE} pool graphs")

    k6_only = []
    pool_k6 = 0
    pool_k7 = 0
    for priority, index in enumerate(pool):
        adj = loaded[index]
        if find_clique_mask(adj, 7):
            pool_k7 += 1
            continue
        if not find_clique_mask(adj, 6):
            raise ValueError(f"pool graph {index} has clique number below six")
        pool_k6 += 1
        if index in killed:
            raise ValueError(
                f"K6-only pool graph {index} occurs in the old kill log"
            )
        if len(k6_only) < SAMPLE_SIZE:
            k6_only.append(
                {
                    "index": index,
                    "selection_priority": priority,
                    "clique_class": "K6_only",
                    "old_status": "deferred",
                    "input_exact_status": "pre_Lorentz_residue",
                    "adjacency": adj,
                }
            )
    if len(k6_only) != SAMPLE_SIZE:
        raise ValueError(
            f"pool supplies only {len(k6_only)} K6-only graphs; need {SAMPLE_SIZE}"
        )

    return {
        "schema": 1,
        "description": (
            "Deterministic uniform-priority sample of 1,024 graphs from the "
            "175,819-graph K6-only exact residue entering the Lorentz filter."
        ),
        "sources": {
            "corpus": {
                "file": corpus.name,
                "sha256": corpus_hash,
                "graphs": corpus_count,
                "tracked": False,
            },
            "old_kill_log": {
                "file": kill_log.name,
                "sha256": kill_hash,
                "tracked": False,
            },
            "previous_exact_profile": {
                "file": profile_path.name,
                "sha256": profile_hash,
                "K6_only_input_residue": EXPECTED_K6_RESIDUE,
            },
        },
        "selection": {
            "algorithm": (
                "first 1,024 K7-free graphs by SplitMix64 pool priority; "
                "outputs reduced modulo 3,971,787; duplicates skipped"
            ),
            "seed_hex": f"0x{SPLITMIX64_SEED:016x}",
            "pool_size": POOL_SIZE,
            "sample_size": SAMPLE_SIZE,
            "pool_classification": {"K6_only": pool_k6, "K7": pool_k7},
            "largest_selected_priority": k6_only[-1]["selection_priority"],
        },
        "preconditions": {
            "n": 19,
            "alpha_at_most_two_checked_per_graph": True,
            "required_K6_present": True,
            "required_K7_absent": True,
            "candidate_nonedges": "unconstrained; allowed defects may be zero",
        },
        "graphs": k6_only,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("aeq_d6_n19.txt"))
    parser.add_argument(
        "--kill-log", type=Path, default=Path("killed_d6_n19.log")
    )
    parser.add_argument(
        "--profile",
        type=Path,
        default=Path("d6_k6_lorentz_input_profile.json"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("d6_k6_lorentz_sample.json")
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify byte identity with --output instead of rewriting it",
    )
    args = parser.parse_args()
    sample = build_sample(args.corpus, args.kill_log, args.profile)
    rendered = json.dumps(sample, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} is not reproducible")
        print(f"PASS: {args.output} is reproducible")
    else:
        args.output.write_text(rendered, encoding="utf-8")
        selected = sample["graphs"]
        print(
            f"wrote {args.output}: {len(selected)} K6-only residue graphs; "
            f"max priority {selected[-1]['selection_priority']}; pool "
            f"{sample['selection']['pool_classification']}"
        )


if __name__ == "__main__":
    main()
