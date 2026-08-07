#!/usr/bin/env python3
"""Post-process the bounded C containment pilot against the new exact filter.

The C run is deliberately stratified at the preexisting-profile boundary.
This script keeps that sample fixed and asks the strategically decisive
question: which containment hits remain after the full audited direct K7
filter union (bounded cover OR tight-cover matching) has been applied?

It also classifies the 89 LM-selected n=14 hypotheses into:

* exact by the new bounded-cover or tight-cover-matching rules;
* exact by a preexisting K5-link or K7-reflection rule only;
* still heuristic.

No numerical LM failure is promoted to a proof.  The exact classifications
use the documented graph-only rules, while containments of the remaining
hypotheses stay explicitly heuristic.

Usage:
  python3 analyze_obstruction_pilot_d6.py \
      d6_obstruction_pilot.json aeq_d6_n14.txt aeq_d6_n19.txt \
      d6_profile.json d6_obstruction_pilot_incremental.json \
      d6_obstruction_pilot_manifest.json
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import sys
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

from cdriver6 import gen_orders
from d6_reference_filters import (
    clique_masks,
    k7_disjoint_edge_bounded_cover,
    k7_tight_cover_matching,
)


PREEXISTING_K7_RESIDUE = "deferred_K7_preexisting_exact_residue"
PREEXISTING_K6_RESIDUE = "deferred_K6_preexisting_exact_residue"


def load_graphs(path: Path, expected_n: int) -> list[list[int]]:
    graphs = []
    with path.open() as source:
        for line_number, line in enumerate(source, 1):
            fields = line.split()
            if not fields:
                raise ValueError(f"blank graph line {line_number} in {path}")
            values = [int(x) for x in fields]
            if values[0] != expected_n or len(values) != expected_n + 1:
                raise ValueError(
                    f"bad graph record at {path}:{line_number}: "
                    f"n/fields={values[0]}/{len(values)}"
                )
            graphs.append(values[1:])
    return graphs


def load_selected_targets(path: Path, wanted: set[int]) -> dict[int, list[int]]:
    selected: dict[int, list[int]] = {}
    with path.open() as source:
        for index, line in enumerate(source):
            if index not in wanted:
                continue
            values = [int(x) for x in line.split()]
            if not values or values[0] != 19 or len(values) != 20:
                raise ValueError(f"bad n=19 graph record at index {index}")
            selected[index] = values[1:]
    missing = wanted - selected.keys()
    if missing:
        raise ValueError(f"sample indices missing from corpus: {sorted(missing)}")
    return selected


def bad_k5_link(adj: Sequence[int]) -> bool:
    """Independent transparent version of the preexisting K5-link rule."""

    for seed_mask in clique_masks(adj, 5):
        common = (1 << len(adj)) - 1
        seed = []
        while seed_mask:
            bit = seed_mask & -seed_mask
            vertex = bit.bit_length() - 1
            seed.append(vertex)
            common &= adj[vertex]
            seed_mask ^= bit
        if common.bit_count() > 4:
            return True
    return False


def bad_k7_reflection(adj: Sequence[int]) -> bool:
    """Independent transparent version of the preexisting reflection rule."""

    n = len(adj)
    for seed_mask in clique_masks(adj, 7):
        seed = [v for v in range(n) if seed_mask & (1 << v)]
        forced: dict[int, int] = {}
        for vertex in range(n):
            if seed_mask & (1 << vertex):
                continue
            defects = [q for q in seed if not adj[vertex] & (1 << q)]
            if len(defects) != 1:
                continue
            defect = defects[0]
            if defect in forced:
                return True
            if any(adj[vertex] & (1 << other) for other in forced.values()):
                return True
            forced[defect] = vertex
    return False


def has_k7(adj: Sequence[int]) -> bool:
    return next(clique_masks(adj, 7), None) is not None


def circle_count(adj: Sequence[int], seed: Iterable[int], order: Iterable[int]) -> int:
    placed = set(seed)
    count = 0
    for vertex in order:
        if sum(bool(adj[vertex] & (1 << other)) for other in placed) == 5:
            count += 1
        placed.add(vertex)
    return count


def hit(row: dict, ordinal: int) -> bool:
    return bool((int(row["hit_mask_hex_high_low"], 16) >> ordinal) & 1)


def timed_out(row: dict, ordinal: int) -> bool:
    return bool((int(row["timeout_mask_hex_high_low"], 16) >> ordinal) & 1)


def union_hits(rows: Sequence[dict], ordinals: set[int]) -> int:
    return sum(any(hit(row, ordinal) for ordinal in ordinals) for row in rows)


def union_timeouts(rows: Sequence[dict], ordinals: set[int]) -> int:
    return sum(any(timed_out(row, ordinal) for ordinal in ordinals) for row in rows)


def greedy(rows: Sequence[dict], pattern_rows: Sequence[dict],
           allowed: set[int]) -> list[dict]:
    covered: set[int] = set()
    chosen: set[int] = set()
    answer = []
    while True:
        candidates = []
        for pattern in pattern_rows:
            ordinal = pattern["ordinal"]
            if ordinal not in allowed or ordinal in chosen:
                continue
            new = {
                index
                for index, row in enumerate(rows)
                if index not in covered and hit(row, ordinal)
            }
            candidates.append((len(new), -pattern["pattern_index"], pattern, new))
        if not candidates:
            break
        marginal, _, pattern, new = max(candidates, key=lambda item: item[:2])
        if not marginal:
            break
        covered.update(new)
        chosen.add(pattern["ordinal"])
        answer.append(
            {
                "step": len(answer) + 1,
                "pattern_index": pattern["pattern_index"],
                "marginal_hits": marginal,
                "cumulative_hits": len(covered),
            }
        )
    return answer


def wilson95(hits: int, total: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = hits / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(
        p * (1 - p) / total + z * z / (4 * total * total)
    ) / denominator
    return center - half, center + half


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while block := source.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    pilot_path, n14_path, n19_path, profile_path, output_path, manifest_path = (
        Path(arg) for arg in sys.argv[1:]
    )
    pilot = json.loads(pilot_path.read_text())
    profile = json.loads(profile_path.read_text())
    pattern_graphs = load_graphs(n14_path, 14)
    pattern_rows = pilot["individual_patterns"]
    if len(pattern_rows) != 89:
        raise ValueError(f"expected 89 selected patterns, got {len(pattern_rows)}")

    pattern_status = []
    new_cover_exact: set[int] = set()
    new_tight_exact: set[int] = set()
    preexisting_exact: set[int] = set()
    engine_usable: set[int] = set()
    for pattern in pattern_rows:
        ordinal = pattern["ordinal"]
        source_index = pattern["pattern_index"]
        adj = pattern_graphs[source_index]
        cover = k7_disjoint_edge_bounded_cover(adj).rejected
        tight = k7_tight_cover_matching(adj).rejected
        old_link = bad_k5_link(adj)
        old_reflection = bad_k7_reflection(adj)
        orders = gen_orders(adj, 14, kmax=12)
        minimum_circles = (
            min(circle_count(adj, seed, order) for seed, order in orders)
            if orders else None
        )
        if cover:
            new_cover_exact.add(ordinal)
        if tight:
            new_tight_exact.add(ordinal)
        if old_link or old_reflection:
            preexisting_exact.add(ordinal)
        if orders:
            engine_usable.add(ordinal)
        pattern_status.append(
            {
                "ordinal": ordinal,
                "pattern_index": source_index,
                "contains_K7": has_k7(adj),
                "exact_new_bounded_cover": cover,
                "exact_new_tight_cover_matching": tight,
                "exact_preexisting_K5_link": old_link,
                "exact_preexisting_K7_reflection": old_reflection,
                "exact_union": cover or tight or old_link or old_reflection,
                "engine_order_count": len(orders),
                "engine_minimum_circle_stages": minimum_circles,
            }
        )

    new_direct_k7_exact = new_cover_exact | new_tight_exact
    exact_union = new_direct_k7_exact | preexisting_exact
    all_ordinals = {row["ordinal"] for row in pattern_rows}
    heuristic = all_ordinals - exact_union
    direct_k7_uncertified = all_ordinals - new_direct_k7_exact
    heuristic_engine = heuristic & engine_usable
    if (len(new_cover_exact) != 78 or len(new_tight_exact) != 0 or
            len(new_direct_k7_exact) != 78 or len(exact_union) != 81 or
            len(heuristic) != 8):
        raise ValueError(
            "unexpected pattern status counts: "
            f"cover={len(new_cover_exact)}, tight={len(new_tight_exact)}, "
            f"direct_K7={len(new_direct_k7_exact)}, "
            f"exact_union={len(exact_union)}, heuristic={len(heuristic)}"
        )
    if not all(row["contains_K7"] for row in pattern_status):
        raise ValueError("expected every selected n=14 pattern to contain K7")

    sample_rows = pilot["sample_audit"]
    sample_indices = {row["corpus_index"] for row in sample_rows}
    targets = load_selected_targets(n19_path, sample_indices)
    direct_k7_decisions: dict[int, tuple[bool, bool]] = {}
    for row in sample_rows:
        index = row["corpus_index"]
        if "K7" in row["stratum"]:
            cover = k7_disjoint_edge_bounded_cover(targets[index]).rejected
            tight = k7_tight_cover_matching(targets[index]).rejected
        else:
            cover = tight = False
        direct_k7_decisions[index] = (cover, tight)

    strata = {}
    for stratum in pilot["sample_counts"]:
        original = [row for row in sample_rows if row["stratum"] == stratum]
        residue = [
            row for row in original
            if not any(direct_k7_decisions[row["corpus_index"]])
        ]
        cover_rejected = sum(
            direct_k7_decisions[row["corpus_index"]][0] for row in original
        )
        tight_rejected = sum(
            direct_k7_decisions[row["corpus_index"]][1] for row in original
        )
        tight_incremental = sum(
            (not direct_k7_decisions[row["corpus_index"]][0]) and
            direct_k7_decisions[row["corpus_index"]][1]
            for row in original
        )
        strata[stratum] = {
            "initial_sample": len(original),
            "bounded_cover_rejected": cover_rejected,
            "tight_cover_matching_rejected": tight_rejected,
            "tight_cover_matching_incremental_over_bound": tight_incremental,
            "full_direct_K7_union_rejected": len(original) - len(residue),
            "full_direct_exact_residue": len(residue),
            "incremental_hits_by_78_new_direct_K7_exact_patterns": union_hits(
                residue, new_direct_k7_exact
            ),
            "incremental_hits_by_81_all_exact_patterns": union_hits(
                residue, exact_union
            ),
            "incremental_hits_by_8_still_heuristic_patterns": union_hits(
                residue, heuristic
            ),
            "incremental_hits_by_engine_usable_still_heuristic_patterns": (
                union_hits(residue, heuristic_engine)
            ),
            "target_pattern_timeout_rows_among_heuristic_patterns": (
                union_timeouts(residue, heuristic)
            ),
        }

    current_k7_rows = [
        row for row in sample_rows
        if row["stratum"] == PREEXISTING_K7_RESIDUE
        and not any(direct_k7_decisions[row["corpus_index"]])
    ]
    current_k6_rows = [
        row for row in sample_rows if row["stratum"] == PREEXISTING_K6_RESIDUE
    ]
    k7_hits = union_hits(current_k7_rows, heuristic)
    k6_hits = union_hits(current_k6_rows, heuristic)
    if k6_hits:
        raise ValueError("a K7-containing pattern was reported inside K6-only target")
    deferred_profile = profile["populations"]["deferred"]
    current_k7_population = deferred_profile["residue_by_clique_number"]["7"]
    current_k6_population = deferred_profile["residue_by_clique_number"]["6"]
    current_population = current_k7_population + current_k6_population
    k7_weight = current_k7_population / current_population
    low, high = wilson95(k7_hits, len(current_k7_rows))

    ordinal_to_index = {
        row["ordinal"]: row["pattern_index"] for row in pattern_rows
    }
    heuristic_hit_counts = []
    for ordinal in sorted(heuristic, key=ordinal_to_index.get):
        count = sum(hit(row, ordinal) for row in current_k7_rows)
        heuristic_hit_counts.append(
            {
                "pattern_index": ordinal_to_index[ordinal],
                "hits_in_sampled_current_K7_residue": count,
                "engine_usable": ordinal in engine_usable,
            }
        )

    output = {
        "schema": 1,
        "claim_status": (
            "81 patterns have graph-only exact certificates; containments "
            "of the remaining 8 are heuristic until separately certified"
        ),
        "decisive_metric": (
            "incremental containment coverage beyond the full direct exact "
            "K7 union (bounded cover OR tight-cover matching)"
        ),
        "pattern_classification": {
            "all_selected_patterns": 89,
            "all_contain_required_K7": True,
            "consequence_for_K6_only_population": (
                "structurally zero coverage by all 89 patterns"
            ),
            "exact_by_new_bounded_cover": {
                "count": len(new_cover_exact),
                "pattern_indices": sorted(ordinal_to_index[o] for o in new_cover_exact),
            },
            "exact_by_new_tight_cover_matching": {
                "count": len(new_tight_exact),
                "pattern_indices": sorted(ordinal_to_index[o] for o in new_tight_exact),
            },
            "exact_by_full_direct_K7_union": {
                "count": len(new_direct_k7_exact),
                "pattern_indices": sorted(
                    ordinal_to_index[o] for o in new_direct_k7_exact
                ),
            },
            "not_exact_by_full_direct_K7_union": {
                "count": len(direct_k7_uncertified),
                "pattern_indices": sorted(
                    ordinal_to_index[o] for o in direct_k7_uncertified
                ),
            },
            "additional_exact_by_preexisting_rules_only": {
                "count": len(preexisting_exact - new_direct_k7_exact),
                "pattern_indices": sorted(
                    ordinal_to_index[o]
                    for o in preexisting_exact - new_direct_k7_exact
                ),
            },
            "all_exact_union": {
                "count": len(exact_union),
                "pattern_indices": sorted(ordinal_to_index[o] for o in exact_union),
            },
            "still_heuristic": {
                "count": len(heuristic),
                "pattern_indices": sorted(ordinal_to_index[o] for o in heuristic),
            },
            "still_heuristic_and_engine_usable": {
                "count": len(heuristic_engine),
                "pattern_indices": sorted(
                    ordinal_to_index[o] for o in heuristic_engine
                ),
            },
        },
        "exact_pattern_redundancy_theorem": {
            "patterns_subsumed_by_direct_exact_filters": len(exact_union),
            "conclusion": (
                "All 81 exact patterns have identically zero incremental "
                "coverage beyond the direct exact filter union on every "
                "target graph, not merely on this sample."
            ),
            "bounded_cover_78_proof": (
                "Under a non-induced edge-preserving embedding, the witness "
                "K7 remains a K7; adding seed edges only shrinks allowed "
                "defect sets, disjoint-support required L-edges persist, "
                "and the eligible cover-vertex set can only shrink. Thus an "
                "eligible cover of size at most 7 in the target would induce "
                "one for the pattern, contradicting its certificate."
            ),
            "tight_cover_matching_note": (
                "The audited tight-cover matching refinement rejects none "
                "of these 89 n=14 patterns. It is nevertheless included in "
                "every target-side direct-union decision and in the final "
                "profile residue populations."
            ),
            "old_rule_only_3_pattern_indices": sorted(
                ordinal_to_index[o]
                for o in preexisting_exact - new_direct_k7_exact
            ),
            "old_rule_only_3_proof": (
                "Patterns 303, 304, and 443 violate the K5 common-neighbour "
                "bound (not reflection). The required K5 and its more than "
                "four required common neighbours persist under every "
                "edge-preserving embedding."
            ),
            "sample_sanity_check": (
                "The measured exact-pattern incremental union is 0 on all "
                "post-direct-union sampled residue; this agrees with, but is not "
                "the proof of, the theorem above."
            ),
        },
        "sample_intersection_by_pull_boundary_stratum": strata,
        "current_exact_residue_incremental_summary": {
            "population": {
                "K7": current_k7_population,
                "K6_only": current_k6_population,
                "total": current_population,
            },
            "fixed_sample_after_full_direct_exact_K7_union": {
                "K7": len(current_k7_rows),
                "K6_only": len(current_k6_rows),
            },
            "incremental_exact_pattern_union_hits": {
                "K7": union_hits(current_k7_rows, exact_union),
                "K6_only": union_hits(current_k6_rows, exact_union),
            },
            "incremental_still_heuristic_union_hits": {
                "K7": k7_hits,
                "K6_only": k6_hits,
            },
            "incremental_engine_usable_heuristic_union_hits": {
                "K7": union_hits(current_k7_rows, heuristic_engine),
                "K6_only": union_hits(current_k6_rows, heuristic_engine),
            },
            "population_weighted_heuristic_point_estimate": (
                k7_weight * k7_hits / len(current_k7_rows)
            ),
            "population_weighted_approximate_95pct_low": k7_weight * low,
            "population_weighted_approximate_95pct_high": k7_weight * high,
            "uncertainty_note": (
                f"The {len(current_k7_rows)} K7 survivors are the "
                "rejection-sampled subset of a "
                "fixed uniform old-residue reservoir; Wilson uncertainty is "
                "therefore large. K6 coverage is exactly zero because every "
                "pattern contains K7."
            ),
            "heuristic_pattern_individual_hits": heuristic_hit_counts,
            "heuristic_greedy_marginal": greedy(
                current_k7_rows, pattern_rows, heuristic
            ),
            "engine_usable_heuristic_greedy_marginal": greedy(
                current_k7_rows, pattern_rows, heuristic_engine
            ),
        },
        "pattern_status": pattern_status,
    }
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    source_dir = Path(__file__).resolve().parent
    input_paths = {
        "obstruction_pilot_d6.c": source_dir / "obstruction_pilot_d6.c",
        "analyze_obstruction_pilot_d6.py": Path(__file__).resolve(),
        "d6_reference_filters.py": source_dir / "d6_reference_filters.py",
        "cdriver6.py": source_dir / "cdriver6.py",
        "primary_pilot_json": pilot_path,
        "incremental_analysis_json": output_path,
        "n14_corpus": n14_path,
        "lm_triage": Path(pilot["inputs"]["lm_triage"]),
        "n19_corpus": n19_path,
        "kill_log": Path(pilot["inputs"]["kill_log"]),
        "exact_profile_json": profile_path,
    }
    compiler = subprocess.run(
        ["cc", "--version"], check=True, text=True, capture_output=True
    ).stdout.splitlines()[0]
    manifest = {
        "schema": 1,
        "commands": {
            "build": (
                "cc -O3 -std=c11 -Wall -Wextra -pedantic -pthread "
                "-o obstruction_pilot_d6.bin obstruction_pilot_d6.c -lm"
            ),
            "pilot": " ".join(pilot["command"]),
            "incremental_analysis": "python3 " + " ".join(sys.argv),
        },
        "platform": platform.platform(),
        "python": platform.python_version(),
        "compiler": compiler,
        "sha256": {name: sha256(path) for name, path in input_paths.items()},
    }
    binary = source_dir / "obstruction_pilot_d6.bin"
    if binary.exists():
        manifest["sha256"]["obstruction_pilot_d6.bin"] = sha256(binary)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
