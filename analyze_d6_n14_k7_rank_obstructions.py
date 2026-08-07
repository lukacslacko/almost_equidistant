#!/usr/bin/env python3
"""Certify LM-selected n=14 patterns with dimension-independent K7 theory.

The level-19 theorem |Z|<=3 is intentionally NOT used here.  For every K7
seed this program enumerates all eligible covers through size seven and uses
only:

* the orthonormal zero-factor support CSP (including square row rules);
* rank(K) <= min(7-|Z|, nu_N, nu_T-|Z|);
* positive-definite required-clique and perpendicular-neighbourhood bounds;
* rank(B) <= rank(K)+1;
* ordinary zero forcing and the exact component/inertia nullity bound for B.

The input patterns have alpha at most two.  A rejection is monotone under
edge addition, so it certifies the required-edge pattern under non-induced
containment; candidate nonedges are never prescribed to be non-unit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

from d6_k7_rank_reference import (
    SupportSolver,
    ZeroForcingSolver,
    clique_masks,
    complement_graph,
    component_inertia_nullity_caps,
    eligible_covers,
    induced_graph,
    matching_size,
    seed_instance,
    validate_graph,
)


EXPECTED_PATTERNS = 89
EXPECTED_N14 = 1052
RULES = (
    "bounded_cover",
    "support",
    "K_rank",
    "K_clique",
    "F_degree",
    "ordinary_B",
    "component_B",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_n14(path: Path) -> list[list[int]]:
    graphs = []
    with path.open(encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if len(fields) != 15 or fields[0] != "14":
                raise ValueError(f"malformed n14 record at line {line_number}")
            adjacency = [int(value) for value in fields[1:]]
            validate_graph(adjacency)
            graphs.append(adjacency)
    if len(graphs) != EXPECTED_N14:
        raise ValueError(f"n14 corpus has {len(graphs)} graphs")
    return graphs


def rank_upper_bounds(
    n_size: int, nu_n: int, nu_t: int, z_size: int
) -> tuple[int, int]:
    """Dimension-independent cross-orthogonal K/B rank upper bounds."""

    k_upper = min(n_size, nu_n, max(0, nu_t - z_size), 7 - z_size)
    return k_upper, min(n_size, k_upper + 1, 8 - z_size)


def clique_number(adjacency: Sequence[int]) -> int:
    """Exact maximum-clique size; here the graph has at most seven vertices."""

    best = 0
    for mask in range(1 << len(adjacency)):
        size = mask.bit_count()
        if size <= best:
            continue
        if all(
            not (mask & ~(1 << vertex) & ~adjacency[vertex])
            for vertex in range(len(adjacency))
            if mask & (1 << vertex)
        ):
            best = size
    return best


def analyze_seed(
    adjacency: Sequence[int],
    seed_mask: int,
    support_solver: SupportSolver,
    zero_forcing: ZeroForcingSolver,
) -> dict:
    seed, outside, defects, ladj, eligible = seed_instance(adjacency, seed_mask)
    covers = eligible_covers(ladj, eligible, cap=7)
    nu_t = matching_size(defects)
    counts: Counter[str] = Counter()
    sizes: Counter[str] = Counter()
    passing = {rule: 0 for rule in RULES[1:]}
    joint_passes = 0
    first_joint_failure = None
    for zmask in covers:
        zlocal = [i for i in range(len(outside)) if zmask & (1 << i)]
        nlocal = [i for i in range(len(outside)) if not (zmask & (1 << i))]
        sizes[str(len(zlocal))] += 1
        support_failed = support_solver.solve(defects[i] for i in zlocal) is None
        ndefects = [defects[i] for i in nlocal]
        nu_n = matching_size(ndefects)
        selected = [outside[i] for i in nlocal]
        graph_n = induced_graph(adjacency, selected)

        # This equality also checks the optional-zero/cover semantics.
        for i in range(len(nlocal)):
            for j in range(i):
                assert bool(graph_n[i] & (1 << j)) == bool(
                    ndefects[i] & ndefects[j]
                )

        k_zf = zero_forcing.solve(graph_n)
        k_upper, b_upper = rank_upper_bounds(
            len(nlocal), nu_n, nu_t, len(zlocal)
        )
        assert len(nlocal) - k_zf.number <= min(len(nlocal), nu_n)
        k_failed = len(nlocal) - k_zf.number > k_upper

        fgraph = complement_graph(graph_n)
        k_clique_failed = clique_number(graph_n) > k_upper
        f_degree_failed = bool(nlocal) and max(
            row.bit_count() for row in fgraph
        ) > k_upper - 1
        b_zf = zero_forcing.solve(fgraph)
        ordinary_b_failed = len(nlocal) - b_zf.number > b_upper
        component_caps, _ = component_inertia_nullity_caps(
            fgraph, zero_forcing
        )
        component_b_failed = (
            len(nlocal) - component_caps["exact"] > b_upper
        )
        assert not ordinary_b_failed or component_b_failed
        failures = {
            "support": support_failed,
            "K_rank": k_failed,
            "K_clique": k_clique_failed,
            "F_degree": f_degree_failed,
            "ordinary_B": ordinary_b_failed,
            "component_B": component_b_failed,
        }
        counts.update({rule: int(failed) for rule, failed in failures.items()})
        for rule, failed in failures.items():
            passing[rule] += int(not failed)
        joint_failed = (
            support_failed
            or k_failed
            or k_clique_failed
            or f_degree_failed
            or component_b_failed
        )
        joint_passes += int(not joint_failed)
        if joint_failed and first_joint_failure is None:
            first_joint_failure = {
                "Z": [outside[i] for i in zlocal],
                "failure_flags": failures,
                "nu_N": nu_n,
                "nu_T": nu_t,
                "K_zero_forcing": k_zf.number,
                "B_zero_forcing": b_zf.number,
                "component_nullity_cap": component_caps["exact"],
            }

    old_rejected = not covers
    return {
        "seed": seed,
        "eligible_covers": len(covers),
        "cover_size_histogram": dict(sizes),
        "cover_failures": dict(counts),
        "bounded_cover_rejected": old_rejected,
        "support_only_rejected": bool(covers) and passing["support"] == 0,
        "K_rank_only_rejected": bool(covers) and passing["K_rank"] == 0,
        "K_clique_only_rejected": (
            bool(covers) and passing["K_clique"] == 0
        ),
        "F_degree_only_rejected": (
            bool(covers) and passing["F_degree"] == 0
        ),
        "ordinary_B_only_rejected": (
            bool(covers) and passing["ordinary_B"] == 0
        ),
        "component_B_only_rejected": (
            bool(covers) and passing["component_B"] == 0
        ),
        "joint_rejected": bool(covers) and joint_passes == 0,
        "joint_passes": joint_passes,
        "first_joint_failure": first_joint_failure,
    }


def analyze_pattern(payload: tuple[dict, list[int]]) -> dict:
    pattern, adjacency = payload
    support_solver = SupportSolver()
    zero_forcing = ZeroForcingSolver()
    started = time.perf_counter()
    seeds = [
        analyze_seed(adjacency, seed, support_solver, zero_forcing)
        for seed in clique_masks(adjacency, 7)
    ]
    if not seeds:
        raise ValueError(f"selected pattern {pattern['pattern_index']} has no K7")
    decisions = {
        "bounded_cover": any(seed["bounded_cover_rejected"] for seed in seeds),
        "support": any(seed["support_only_rejected"] for seed in seeds),
        "K_rank": any(seed["K_rank_only_rejected"] for seed in seeds),
        "K_clique": any(seed["K_clique_only_rejected"] for seed in seeds),
        "F_degree": any(seed["F_degree_only_rejected"] for seed in seeds),
        "ordinary_B": any(seed["ordinary_B_only_rejected"] for seed in seeds),
        "component_B": any(
            seed["component_B_only_rejected"] for seed in seeds
        ),
        "joint": any(seed["joint_rejected"] for seed in seeds),
    }
    first = {
        rule: next(
            (
                seed["seed"]
                for seed in seeds
                if seed[
                    "bounded_cover_rejected"
                    if rule == "bounded_cover"
                    else f"{rule}_only_rejected"
                ]
            ),
            None,
        )
        for rule in RULES
    }
    cover_sizes: Counter[str] = Counter()
    cover_failures: Counter[str] = Counter()
    for seed in seeds:
        cover_sizes.update(seed["cover_size_histogram"])
        cover_failures.update(seed["cover_failures"])
    return {
        "ordinal": pattern["ordinal"],
        "pattern_index": pattern["pattern_index"],
        "edges": pattern["edges"],
        "lm_distinct_residual": pattern["lm_distinct_residual"],
        "seeds_checked": len(seeds),
        "eligible_covers_checked": sum(seed["eligible_covers"] for seed in seeds),
        "cover_size_histogram": dict(cover_sizes),
        "cover_failures": dict(cover_failures),
        "decisions": decisions,
        "first_rejecting_seed": {
            rule: seed for rule, seed in first.items() if seed is not None
        },
        "elapsed_seconds": time.perf_counter() - started,
        "support_CSP_nodes": support_solver.nodes,
        "zero_forcing_initial_sets_checked": zero_forcing.initial_sets_checked,
    }


def hit(sample_row: dict, ordinal: int) -> bool:
    return bool((int(sample_row["hit_mask_hex_high_low"], 16) >> ordinal) & 1)


def coverage(pilot: dict, exact_ordinals: set[int]) -> dict:
    rows = pilot["sample_audit"]
    by_stratum = {}
    for stratum, total in pilot["sample_counts"].items():
        selected = [row for row in rows if row["stratum"] == stratum]
        hits = sum(
            any(hit(row, ordinal) for ordinal in exact_ordinals)
            for row in selected
        )
        assert len(selected) == total
        by_stratum[stratum] = {"hits": hits, "sample": total}
    union = sum(
        any(hit(row, ordinal) for ordinal in exact_ordinals) for row in rows
    )
    return {
        "known_fixed_sample_union_hits": union,
        "known_fixed_sample_size": len(rows),
        "by_stratum": by_stratum,
        "interpretation": (
            "Exact on these fixed containment hits, but the deliberately "
            "stratified sample is not an unweighted population estimate."
        ),
    }


def build_report(
    pilot: dict,
    prior: dict,
    n14: list[list[int]],
    workers: int,
) -> dict:
    selected = pilot["individual_patterns"]
    if len(selected) != EXPECTED_PATTERNS:
        raise ValueError(f"pilot has {len(selected)} selected patterns")
    payloads = [(row, n14[row["pattern_index"]]) for row in selected]
    started = time.perf_counter()
    if workers == 1:
        results = [analyze_pattern(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(analyze_pattern, payloads, chunksize=1))
    wall = time.perf_counter() - started

    decision_sets = {
        rule: {
            result["ordinal"] for result in results
            if result["decisions"][rule]
        }
        for rule in (*RULES, "joint")
    }
    exact_new = decision_sets["bounded_cover"] | decision_sets["joint"]
    prior_status = {row["ordinal"]: row for row in prior["pattern_status"]}
    prior_exact = {
        ordinal for ordinal, row in prior_status.items() if row["exact_union"]
    }
    newly_exact = exact_new - prior_exact
    combined_exact = exact_new | prior_exact
    all_ordinals = {row["ordinal"] for row in selected}
    still_heuristic = all_ordinals - combined_exact

    individual = []
    pilot_by_ordinal = {row["ordinal"]: row for row in selected}
    result_by_ordinal = {row["ordinal"]: row for row in results}
    for ordinal in sorted(combined_exact):
        source = pilot_by_ordinal[ordinal]
        result = result_by_ordinal[ordinal]
        individual.append(
            {
                "ordinal": ordinal,
                "pattern_index": source["pattern_index"],
                "new_K7_decisions": result["decisions"],
                "previously_exact": ordinal in prior_exact,
                "newly_exact": ordinal in newly_exact,
                "known_containment_hits": source["hits"],
                "known_hits_by_stratum": source["hits_by_stratum"],
            }
        )

    cover_sizes: Counter[str] = Counter()
    cover_failures: Counter[str] = Counter()
    for result in results:
        cover_sizes.update(result["cover_size_histogram"])
        cover_failures.update(result["cover_failures"])
    by_index = {row["pattern_index"]: row for row in individual}
    return {
        "schema": 1,
        "claim_status": (
            "Exact finite graph/algebra certificates for listed patterns; "
            "LM residuals are retained only as historical triage data."
        ),
        "n14_scope": {
            "outside_vertices_per_K7": 7,
            "covers_enumerated": "every eligible cover with size at most 7",
            "level19_cap3_used": False,
            "candidate_nonedge_semantics": (
                "optional unit pairs; exactness follows by monotonicity over "
                "all alpha-at-most-two edge supergraphs"
            ),
        },
        "workers": workers,
        "patterns_checked": len(results),
        "seeds_checked": sum(row["seeds_checked"] for row in results),
        "eligible_covers_checked": sum(
            row["eligible_covers_checked"] for row in results
        ),
        "cover_size_histogram": dict(
            sorted(cover_sizes.items(), key=lambda item: int(item[0]))
        ),
        "individual_cover_failures": dict(cover_failures),
        "decision_counts": {
            rule: len(ordinals) for rule, ordinals in decision_sets.items()
        },
        "decision_pattern_indices": {
            rule: sorted(
                pilot_by_ordinal[ordinal]["pattern_index"]
                for ordinal in ordinals
            )
            for rule, ordinals in decision_sets.items()
        },
        "classification": {
            "new_K7_exact_count": len(exact_new),
            "new_K7_exact_pattern_indices": sorted(
                pilot_by_ordinal[o]["pattern_index"] for o in exact_new
            ),
            "newly_exact_over_previous_union_count": len(newly_exact),
            "newly_exact_over_previous_union_pattern_indices": sorted(
                pilot_by_ordinal[o]["pattern_index"] for o in newly_exact
            ),
            "combined_exact_count": len(combined_exact),
            "still_heuristic_count": len(still_heuristic),
            "still_heuristic_pattern_indices": sorted(
                pilot_by_ordinal[o]["pattern_index"] for o in still_heuristic
            ),
        },
        "frequent_patterns": {
            str(index): by_index.get(index, {
                "pattern_index": index,
                "exact": False,
                "known_containment_hits": next(
                    row["hits"] for row in selected
                    if row["pattern_index"] == index
                ),
            })
            for index in (569, 571)
        },
        "known_containment_coverage": coverage(pilot, combined_exact),
        "monotonicity_certificate": {
            "defect_masks": "edge additions only shrink every D_x",
            "covers": (
                "eligible target covers remain eligible pattern covers; "
                "pattern L edges persist"
            ),
            "support": "shrinking support domains cannot repair infeasibility",
            "rank": (
                "on target N, upgraded pattern nonedges cannot remain in N; "
                "hence K/F patterns agree while nu_N and nu_T only decrease"
            ),
            "positive_definite_submatrices": (
                "target and pattern K/F graphs agree on N while the target "
                "rank upper bound can only decrease; clique and degree "
                "violations therefore persist"
            ),
            "consequence": (
                "a rejected pattern has no realizable alpha-at-most-two edge "
                "supergraph and is valid for non-induced containment"
            ),
        },
        "runtime": {
            "wall_seconds": wall,
            "sum_pattern_seconds": sum(row["elapsed_seconds"] for row in results),
            "maximum_pattern_seconds": max(row["elapsed_seconds"] for row in results),
            "support_CSP_nodes": sum(row["support_CSP_nodes"] for row in results),
            "zero_forcing_initial_sets_checked": sum(
                row["zero_forcing_initial_sets_checked"] for row in results
            ),
        },
        "exact_patterns": individual,
        "per_pattern_reference": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, default=Path("d6_obstruction_pilot.json"))
    parser.add_argument(
        "--prior",
        type=Path,
        default=Path("d6_obstruction_pilot_incremental.json"),
    )
    parser.add_argument("--n14", type=Path, default=Path("aeq_d6_n14.txt"))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--output", type=Path, default=Path("d6_n14_k7_rank_obstructions.json")
    )
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("workers must be positive")
    pilot = json.loads(args.pilot.read_text(encoding="utf-8"))
    prior = json.loads(args.prior.read_text(encoding="utf-8"))
    report = build_report(pilot, prior, load_n14(args.n14), args.workers)
    report["inputs"] = {
        "pilot": {"file": args.pilot.name, "sha256": sha256(args.pilot)},
        "prior": {"file": args.prior.name, "sha256": sha256(args.prior)},
        "n14": {"file": args.n14.name, "sha256": sha256(args.n14)},
        "shared_reference": {
            "file": "d6_k7_rank_reference.py",
            "sha256": sha256(Path("d6_k7_rank_reference.py")),
        },
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"wrote {args.output}: new K7 exact {report['classification']['new_K7_exact_count']}, "
        f"new over prior {report['classification']['newly_exact_over_previous_union_count']}, "
        f"still heuristic {report['classification']['still_heuristic_count']}"
    )


if __name__ == "__main__":
    main()
