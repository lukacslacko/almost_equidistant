#!/usr/bin/env python3
"""Extract the deterministic real-graph sample for d6_reference_filters.py.

The extractor hashes and scans the complete local corpus, but performs clique
classification only on a deterministic 2,048-index SplitMix64 pool.  It then
keeps 16 graphs from each nonempty campaign stratum:

* K7 and old-certified;
* K7 and old-deferred;
* K6-only and old-deferred.

It additionally includes four explicitly named deferred witnesses for the K7
tight-cover matching rule.  Those form a separate targeted stratum and are
not presented as random sample members.

The old campaign contains no K6-only certified graph, so that empty stratum is
recorded explicitly rather than fabricated.  Candidate indices are zero based.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from d6_reference_filters import find_clique_mask, validate_graph


EXPECTED_GRAPHS = 3_971_787
EXPECTED_CORPUS_SHA256 = (
    "12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5"
)
EXPECTED_KILL_LOG_SHA256 = (
    "11c790b12d3543f4476f1eb4c221203f806f1c434f06360b7d25346cace5e305"
)
POOL_SIZE = 2_048
PER_STRATUM = 16
SPLITMIX64_SEED = 0xD619_AE05_6C0D_E55E
MASK64 = (1 << 64) - 1
TARGETED_TIGHT_INDICES = (2_751_055, 2_887_126, 3_368_019, 3_936_501)
TARGETED_STRATUM = "K7_deferred_tight_targeted"


def selection_pool() -> list[int]:
    """Return a platform-independent pseudorandom sequence of unique indices."""

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


def hash_and_load_kill_membership(path: Path, wanted: set[int]) -> tuple[str, set[int]]:
    """Hash the full kill log while retaining membership only for ``wanted``."""

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


def hash_and_load_candidates(path: Path, wanted: set[int]) -> tuple[str, dict[int, list[int]], int]:
    """Hash the full corpus and parse just the selected graph lines."""

    digest = hashlib.sha256()
    selected = {}
    count = 0
    with path.open("rb") as stream:
        for index, raw in enumerate(stream):
            digest.update(raw)
            if index in wanted:
                fields = raw.split()
                if len(fields) != 20 or fields[0] != b"19":
                    raise ValueError(f"malformed candidate line at index {index}")
                adj = [int(field) for field in fields[1:]]
                validate_graph(adj)
                selected[index] = adj
            count = index + 1
    return digest.hexdigest(), selected, count


def build_sample(corpus: Path, kill_log: Path) -> dict:
    """Build the deterministic JSON-serializable sample object."""

    pool = selection_pool()
    pool_set = set(pool)
    selected_indices = pool_set | set(TARGETED_TIGHT_INDICES)
    kill_hash, killed = hash_and_load_kill_membership(kill_log, selected_indices)
    corpus_hash, selected, corpus_count = hash_and_load_candidates(
        corpus, selected_indices
    )
    if corpus_count != EXPECTED_GRAPHS:
        raise ValueError(
            f"corpus has {corpus_count} graphs, expected {EXPECTED_GRAPHS}"
        )
    if corpus_hash != EXPECTED_CORPUS_SHA256:
        raise ValueError(f"unexpected corpus SHA-256 {corpus_hash}")
    if kill_hash != EXPECTED_KILL_LOG_SHA256:
        raise ValueError(f"unexpected kill-log SHA-256 {kill_hash}")
    if len(selected) != len(selected_indices):
        raise ValueError(
            f"loaded {len(selected)} of {len(selected_indices)} selected graphs"
        )

    grouped: dict[str, list[dict]] = {
        "K7_certified": [],
        "K7_deferred": [],
        "K6_only_certified": [],
        "K6_only_deferred": [],
    }
    full = (1 << 19) - 1
    for priority, index in enumerate(pool):
        adj = selected[index]
        clique_class = "K7" if find_clique_mask(adj, full, 7) else "K6_only"
        status = "certified" if index in killed else "deferred"
        stratum = f"{clique_class}_{status}"
        grouped[stratum].append(
            {
                "index": index,
                "selection_priority": priority,
                "status": status,
                "clique_class": clique_class,
                "stratum": stratum,
                "adjacency": adj,
            }
        )

    random_strata = ("K7_certified", "K7_deferred", "K6_only_deferred")
    for stratum in random_strata:
        if len(grouped[stratum]) < PER_STRATUM:
            raise ValueError(
                f"selection pool has only {len(grouped[stratum])} graphs in "
                f"{stratum}; need {PER_STRATUM}"
            )
    if grouped["K6_only_certified"]:
        raise ValueError(
            "found a K6-only old-certified graph, contrary to the full profile"
        )

    graphs = []
    for stratum in random_strata:
        graphs.extend(grouped[stratum][:PER_STRATUM])
    for targeted_order, index in enumerate(TARGETED_TIGHT_INDICES):
        adj = selected[index]
        if index in killed:
            raise ValueError(f"targeted tight witness {index} is old-certified")
        if not find_clique_mask(adj, full, 7):
            raise ValueError(f"targeted tight witness {index} has no K7")
        graphs.append(
            {
                "index": index,
                "selection_priority": None,
                "targeted_order": targeted_order,
                "target_reason": "full-profiler K7 tight-cover rejection witness",
                "status": "deferred",
                "clique_class": "K7",
                "stratum": TARGETED_STRATUM,
                "adjacency": adj,
            }
        )
    return {
        "schema": 2,
        "description": (
            "Fixed plain-Python cross-check sample for the dimension-six "
            "K7/K6 exact graph filters, with a separate targeted tight-cover "
            "witness stratum."
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
        },
        "selection": {
            "algorithm": "SplitMix64 outputs modulo 3971787, duplicates skipped",
            "seed_hex": f"0x{SPLITMIX64_SEED:016x}",
            "pool_size": POOL_SIZE,
            "graphs_per_nonempty_stratum": PER_STRATUM,
            "ordering": "increasing SplitMix64 output priority within each stratum",
        },
        "targeted_selection": {
            "stratum": TARGETED_STRATUM,
            "indices": list(TARGETED_TIGHT_INDICES),
            "reason": (
                "Deterministic full-profiler witnesses included to cross-check "
                "tight-cover rejections, not random coverage."
            ),
            "expected_overlap_partition": {
                "also_old_cover_rejected": list(TARGETED_TIGHT_INDICES[:2]),
                "incremental_over_old_cover": list(TARGETED_TIGHT_INDICES[2:]),
            },
        },
        "selection_pool_stratum_counts": {
            stratum: len(items) for stratum, items in grouped.items()
        },
        "stratum_counts": {
            "K7_certified": PER_STRATUM,
            "K7_deferred": PER_STRATUM,
            "K6_only_certified": 0,
            "K6_only_deferred": PER_STRATUM,
            TARGETED_STRATUM: len(TARGETED_TIGHT_INDICES),
        },
        "empty_stratum_note": (
            "The completed full profile found every old-certified graph had a "
            "K7 and every K6-only graph was deferred."
        ),
        "graphs": graphs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("aeq_d6_n19.txt"))
    parser.add_argument(
        "--kill-log", type=Path, default=Path("killed_d6_n19.log")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("d6_reference_sample.json")
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that --output is byte-identical instead of rewriting it",
    )
    args = parser.parse_args()
    sample = build_sample(args.corpus, args.kill_log)
    rendered = json.dumps(sample, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} is not reproducible")
        print(f"PASS: {args.output} is reproducible")
    else:
        args.output.write_text(rendered, encoding="utf-8")
        print(
            f"wrote {args.output}: {len(sample['graphs'])} graphs; pool strata "
            f"{sample['selection_pool_stratum_counts']}"
        )


if __name__ == "__main__":
    main()
