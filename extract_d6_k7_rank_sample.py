#!/usr/bin/env python3
"""Extract a deterministic sample of the current exact K7 residue.

Selection takes the first 512 qualifying graphs in a fixed 32,768-index
SplitMix64 pool.  A qualifying graph is absent from the old interval kill log,
has profiler clique number seven, and passes every exact K7 flag preceding the
new support/rank layer.  The extractor hashes and validates all four source
artifacts and verifies the full 113,136-graph input-residue count.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from d6_k7_rank_reference import clique_masks, validate_graph


EXPECTED_GRAPHS = 3_971_787
EXPECTED_K7_RESIDUE = 113_136
EXPECTED_CORPUS_SHA256 = (
    "12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5"
)
EXPECTED_KILL_LOG_SHA256 = (
    "11c790b12d3543f4476f1eb4c221203f806f1c434f06360b7d25346cace5e305"
)
EXPECTED_PROFILE_SHA256 = (
    "77d507fcc55c4e34a44af2de1ce7c16d0743476ea09d08ce0215e65c80ab17d9"
)
EXPECTED_DECISIONS_SHA256 = (
    "c2e39f7b366a67c2ea117531e2cd5788a0f40eca32e605601bec59152aecea75"
)
SAMPLE_SIZE = 512
POOL_SIZE = 32_768
SPLITMIX64_SEED = 0xD607_5A77_C0DE_2026
MASK64 = (1 << 64) - 1
EXACT_K7_COLUMNS = (
    "link_mask",
    "reflection",
    "defect_csp",
    "K7_Hall",
    "K7_cover",
    "K7_tight_cover",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_pool() -> list[int]:
    """Return a platform-independent sequence of unique corpus indices."""

    state = SPLITMIX64_SEED
    answer: list[int] = []
    seen: set[int] = set()
    while len(answer) < POOL_SIZE:
        state = (state + 0x9E3779B97F4A7C15) & MASK64
        value = state
        value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
        value ^= value >> 31
        index = value % EXPECTED_GRAPHS
        if index not in seen:
            seen.add(index)
            answer.append(index)
    return answer


def hash_and_load_killed(path: Path) -> tuple[str, bytearray, int]:
    """Hash and load exact old-kill membership for all corpus indices."""

    digest = hashlib.sha256()
    killed = bytearray(EXPECTED_GRAPHS)
    unique = 0
    with path.open("rb") as stream:
        for line_number, raw in enumerate(stream, 1):
            digest.update(raw)
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(f"malformed kill-log line {line_number}")
            index = int(fields[0])
            if not 0 <= index < EXPECTED_GRAPHS:
                raise ValueError(f"bad kill-log index {index}")
            if not killed[index]:
                killed[index] = 1
                unique += 1
    return digest.hexdigest(), killed, unique


def hash_and_classify_decisions(
    path: Path, killed: bytearray, pool: set[int]
) -> tuple[str, set[int], int, dict[str, int]]:
    """Verify all decisions and retain current-residue membership in the pool."""

    digest = hashlib.sha256()
    selected_pool: set[int] = set()
    residue_count = 0
    classification = {
        "current_K7_residue": 0,
        "old_killed": 0,
        "K6_only": 0,
        "K7_exact_rejected": 0,
    }
    with path.open("rb") as stream:
        raw_header = stream.readline()
        digest.update(raw_header)
        header = raw_header.decode("ascii").rstrip("\n").split("\t")
        position = {name: i for i, name in enumerate(header)}
        required = {"index", "omega", *EXACT_K7_COLUMNS}
        missing = required - position.keys()
        if missing:
            raise ValueError(f"decision header missing {sorted(missing)}")
        count = 0
        for expected_index, raw in enumerate(stream):
            digest.update(raw)
            fields = raw.decode("ascii").rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(f"bad decision row {expected_index}")
            index = int(fields[position["index"]])
            if index != expected_index:
                raise ValueError(
                    f"decision index {index}, expected {expected_index}"
                )
            omega = int(fields[position["omega"]])
            exact_clear = all(
                int(fields[position[name]]) == 0 for name in EXACT_K7_COLUMNS
            )
            qualifies = omega == 7 and exact_clear and not killed[index]
            if qualifies:
                residue_count += 1
                if index in pool:
                    selected_pool.add(index)
            if index in pool:
                if killed[index]:
                    classification["old_killed"] += 1
                elif omega == 6:
                    classification["K6_only"] += 1
                elif qualifies:
                    classification["current_K7_residue"] += 1
                else:
                    classification["K7_exact_rejected"] += 1
            count = expected_index + 1
    if count != EXPECTED_GRAPHS:
        raise ValueError(f"decisions have {count}, expected {EXPECTED_GRAPHS}")
    return digest.hexdigest(), selected_pool, residue_count, classification


def hash_and_load_candidates(
    path: Path, wanted: set[int]
) -> tuple[str, dict[int, list[int]], int]:
    """Hash every corpus byte and parse only selected sample candidates."""

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
                adjacency = [int(field) for field in fields[1:]]
                validate_graph(adjacency)
                if next(clique_masks(adjacency, 7), None) is None:
                    raise ValueError(f"selected graph {index} has no K7")
                selected[index] = adjacency
            count = index + 1
    return digest.hexdigest(), selected, count


def check_profile(path: Path) -> tuple[str, dict]:
    profile_hash = sha256(path)
    if profile_hash != EXPECTED_PROFILE_SHA256:
        raise ValueError(f"unexpected profile SHA-256 {profile_hash}")
    with path.open(encoding="utf-8") as stream:
        profile = json.load(stream)
    residue = profile["populations"]["deferred"]["residue_by_clique_number"]
    if residue["7"] != EXPECTED_K7_RESIDUE:
        raise ValueError(f"profile K7 residue is {residue['7']}")
    return profile_hash, profile


def build_sample(
    corpus: Path, kill_log: Path, profile_path: Path, decisions: Path
) -> dict:
    pool = selection_pool()
    pool_set = set(pool)
    kill_hash, killed, killed_count = hash_and_load_killed(kill_log)
    decisions_hash, residue_pool, residue_count, classification = (
        hash_and_classify_decisions(decisions, killed, pool_set)
    )
    profile_hash, _ = check_profile(profile_path)
    if kill_hash != EXPECTED_KILL_LOG_SHA256:
        raise ValueError(f"unexpected kill-log SHA-256 {kill_hash}")
    if decisions_hash != EXPECTED_DECISIONS_SHA256:
        raise ValueError(f"unexpected decisions SHA-256 {decisions_hash}")
    if residue_count != EXPECTED_K7_RESIDUE:
        raise ValueError(
            f"decisions give K7 residue {residue_count}, expected "
            f"{EXPECTED_K7_RESIDUE}"
        )

    chosen: list[tuple[int, int]] = []
    for priority, index in enumerate(pool):
        if index in residue_pool and len(chosen) < SAMPLE_SIZE:
            chosen.append((priority, index))
    if len(chosen) != SAMPLE_SIZE:
        raise ValueError(
            f"pool supplies {len(chosen)} current K7 residue graphs"
        )
    wanted = {index for _, index in chosen}
    corpus_hash, loaded, corpus_count = hash_and_load_candidates(corpus, wanted)
    if corpus_hash != EXPECTED_CORPUS_SHA256:
        raise ValueError(f"unexpected corpus SHA-256 {corpus_hash}")
    if corpus_count != EXPECTED_GRAPHS:
        raise ValueError(f"corpus has {corpus_count}, expected {EXPECTED_GRAPHS}")
    if len(loaded) != SAMPLE_SIZE:
        raise ValueError(f"loaded {len(loaded)} of {SAMPLE_SIZE} graphs")

    graphs = [
        {
            "index": index,
            "selection_priority": priority,
            "clique_class": "K7",
            "old_status": "deferred",
            "stratum": "K7_current_exact_residue",
            "input_exact_status": "pre_support_rank_residue",
            "adjacency": loaded[index],
        }
        for priority, index in chosen
    ]
    return {
        "schema": 1,
        "description": (
            "Deterministic uniform-priority sample of 512 graphs from the "
            "113,136-graph current exact K7 residue entering the existential "
            "support/B-rank reference."
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
                "unique_indices": killed_count,
                "tracked": False,
            },
            "current_exact_profile": {
                "file": profile_path.name,
                "sha256": profile_hash,
                "K7_input_residue": residue_count,
            },
            "current_exact_decisions": {
                "file": decisions.name,
                "sha256": decisions_hash,
                "tracked": False,
                "required_zero_columns": list(EXACT_K7_COLUMNS),
            },
        },
        "selection": {
            "algorithm": (
                "first 512 current-K7-residue graphs by SplitMix64 pool "
                "priority; outputs reduced modulo 3,971,787; duplicates skipped"
            ),
            "seed_hex": f"0x{SPLITMIX64_SEED:016x}",
            "pool_size": POOL_SIZE,
            "sample_size": SAMPLE_SIZE,
            "pool_classification": classification,
            "largest_selected_priority": graphs[-1]["selection_priority"],
        },
        "preconditions": {
            "n": 19,
            "alpha_at_most_two_checked_per_graph": True,
            "required_K7_present": True,
            "old_kill_absent": True,
            "all_prior_exact_K7_flags_clear": True,
            "candidate_nonedges": "unconstrained; allowed defects may be zero",
        },
        "graphs": graphs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("aeq_d6_n19.txt"))
    parser.add_argument(
        "--kill-log", type=Path, default=Path("killed_d6_n19.log")
    )
    parser.add_argument(
        "--profile", type=Path, default=Path("d6_k7_rank_input_profile.json")
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=Path("/tmp/d6_profile_stream_decisions.tsv"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("d6_k7_rank_sample.json")
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    sample = build_sample(
        args.corpus, args.kill_log, args.profile, args.decisions
    )
    rendered = json.dumps(sample, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} is not reproducible")
        print(f"PASS: {args.output} is reproducible")
    else:
        args.output.write_text(rendered, encoding="utf-8")
        print(
            f"wrote {args.output}: {len(sample['graphs'])} current K7 residue "
            f"graphs; max priority "
            f"{sample['selection']['largest_selected_priority']}; pool "
            f"{sample['selection']['pool_classification']}"
        )


if __name__ == "__main__":
    main()
