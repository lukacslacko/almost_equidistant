#!/usr/bin/env python3
"""Exact K6 arbitrary-subset rank Hall on the pinned 831-graph residue.

For a positive-sign nonbipartite connected side span F, every selected subset
S is allowed as one alternative block.  Components H of F[S] are mutually
Gram-orthogonal and contribute the maximum of ordinary zero-forcing,
positive-sign inertia, and the PSD--Z |H|-1 bound when H is bipartite.
Exactly one subset is selected per original connected span; subset ranks from
the same span are never added.  The zero-coverage one-vertex extension from
the discovery probe is deliberately omitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import d6_k6_fused_side_rank as frozen
import d6_k6_psd_z_hereditary as parent
import d6_k6_psd_zmatrix as rank_parent
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    build_instance,
    clique_masks,
    find_clique_mask,
    support_matching,
    validate_graph,
    vertices,
)
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-arbitrary-subset-hall-v1"
CERTIFICATE_SCHEMA = "d6-k6-arbitrary-subset-hall-certificates-v1"
CHECKPOINT_SCHEMA = "d6-k6-arbitrary-subset-hall-checkpoint-v1"
EXPECTED_INPUT = 831
EXPECTED_INPUT_SHA256 = (
    "2bcdad095c6bd3038a4bb1d117f2faadfa98113d814633351ec9e79ff553a44b"
)
EXPECTED = {
    "d6_k6_psd_z_hereditary.py": (
        "7b869895cb9b4f7f3bbf1a1d1b11ec71507b6355c8ec0b7d7f8eaec2a51b180d"
    ),
    "d6_k6_psd_z_hereditary_report.json": (
        "202a844d6505d3c68c9a0bea8a5d82a983a7e44b96cb9f3c711d80c47f6c00c1"
    ),
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
}


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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def components_in_mask(rows: Sequence[int], selected: int) -> tuple[int, ...]:
    remaining = selected
    answer = []
    while remaining:
        root = remaining & -remaining
        reached = root
        frontier = root
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            new = rows[vertex] & remaining & ~reached
            reached |= new
            frontier |= new
        answer.append(reached)
        remaining &= ~reached
    return tuple(answer)


def induced_rows(rows: Sequence[int], selected: int) -> tuple[int, ...]:
    chosen = vertices(selected)
    positions = {vertex: local for local, vertex in enumerate(chosen)}
    return tuple(
        sum(
            1 << positions[other]
            for other in chosen
            if rows[vertex] & (1 << other)
        )
        for vertex in chosen
    )


def is_bipartite(rows: Sequence[int]) -> bool:
    colors: dict[int, int] = {}
    for root in range(len(rows)):
        if root in colors:
            continue
        colors[root] = 0
        stack = [root]
        while stack:
            vertex = stack.pop()
            for other in vertices(rows[vertex]):
                if other in colors:
                    if colors[other] == colors[vertex]:
                        return False
                else:
                    colors[other] = 1 - colors[vertex]
                    stack.append(other)
    return True


class PositiveSubsetRank:
    def __init__(self, inertia: InertiaCache, zero_forcing: ZeroForcingSolver):
        self.inertia = inertia
        self.zero_forcing = zero_forcing
        self.values: dict[tuple[tuple[int, ...], int], int] = {}
        self.hits = 0

    def solve(self, full_rows: Sequence[int], selected: int) -> int:
        key = tuple(full_rows), selected
        if key in self.values:
            self.hits += 1
            return self.values[key]
        rank = 0
        for component in components_in_mask(full_rows, selected):
            local = induced_rows(full_rows, component)
            size = len(local)
            pattern = tuple(row | (1 << i) for i, row in enumerate(local))
            _, negative, _ = self.inertia.solve(pattern)
            zf_rank = size - self.zero_forcing.solve(local).number
            inertia_rank = size - negative
            psd_z_rank = size - 1 if is_bipartite(local) else 0
            rank += max(zf_rank, inertia_rank, psd_z_rank)
        self.values[key] = rank
        return rank


@dataclass(frozen=True)
class GroupStats:
    raw_options: int = 0
    compressed_options: int = 0
    arbitrary_groups: int = 0
    arbitrary_proper_options: int = 0


def allowed_union(masks: Sequence[int], selected: int) -> int:
    answer = 0
    for vertex in vertices(selected):
        answer |= masks[vertex]
    return answer


def compress_choices(
    label: str,
    raw: Sequence[parent.Choice],
    applicable: bool,
) -> parent.Group:
    best: dict[int, parent.Choice] = {}
    for choice in raw:
        old = best.get(choice.coordinates)
        if old is None or choice.rank > old.rank:
            best[choice.coordinates] = choice
    return parent.Group(
        label,
        tuple(
            best[mask]
            for mask in sorted(best, key=lambda item: (item.bit_count(), item))
        ),
        applicable,
    )


def make_side_groups(
    instance: K6LorentzInstance,
    label: str,
    side,
    subset_rank: PositiveSubsetRank,
) -> tuple[tuple[parent.Group, ...], GroupStats]:
    local_by_absolute = {
        absolute: local for local, absolute in enumerate(instance.outside)
    }
    groups = []
    stats = GroupStats()
    raw_options = compressed_options = arbitrary_groups = proper_options = 0
    for number, detail in enumerate(side.components):
        absolute = tuple(detail["vertices"])
        masks = tuple(
            instance.defects[local_by_absolute[vertex]] for vertex in absolute
        )
        size = len(absolute)
        applicable = bool(detail["psd_zmatrix_applicable"])
        full_rank = int(detail["componentwise_fused_rank_lower"])
        if applicable:
            group = parent.component_group(
                f"{label}:F{number}", absolute, masks, full_rank, True
            )
            raw_count = 1 << size
        else:
            raw = [parent.Choice(0, 0, (), "empty")]
            for selected in range(1, 1 << size):
                full = selected == (1 << size) - 1
                rank = subset_rank.solve(detail["adjacency_rows"], selected)
                if full and rank != full_rank:
                    raise AssertionError(
                        "arbitrary full-subset rank differs from frozen bound"
                    )
                raw.append(parent.Choice(
                    rank,
                    allowed_union(masks, selected),
                    tuple(absolute[i] for i in vertices(selected)),
                    "full" if full else "arbitrary_positive_nonbip",
                ))
            group = compress_choices(f"{label}:F{number}", raw, False)
            raw_count = len(raw)
            arbitrary_groups += 1
            proper_options += len(raw) - 2
        groups.append(group)
        raw_options += raw_count
        compressed_options += len(group.choices)
    stats = GroupStats(
        raw_options, compressed_options, arbitrary_groups, proper_options
    )
    return tuple(groups), stats


def orientation_result(
    instance: K6LorentzInstance,
    z0: int,
    side_a,
    side_b,
    subset_rank: PositiveSubsetRank,
) -> tuple[parent.HallResult, GroupStats]:
    groups_a, stats_a = make_side_groups(instance, "A", side_a, subset_rank)
    groups_b, stats_b = make_side_groups(instance, "B", side_b, subset_rank)
    z0_groups = tuple(
        parent.singleton_group(
            f"Z0:{instance.outside[local]}",
            instance.outside[local],
            instance.defects[local],
        )
        for local in vertices(z0)
    )
    stats = GroupStats(
        stats_a.raw_options + stats_b.raw_options + 2 * len(z0_groups),
        stats_a.compressed_options + stats_b.compressed_options
        + 2 * len(z0_groups),
        stats_a.arbitrary_groups + stats_b.arbitrary_groups,
        stats_a.arbitrary_proper_options + stats_b.arbitrary_proper_options,
    )
    return parent.check_groups(groups_a + groups_b + z0_groups), stats


def check_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
    subset_rank: PositiveSubsetRank,
) -> dict:
    graph_a, absolute_a = induced_required_graph(adj, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adj, instance, component.side_b)
    orientations = []
    transitions = raw_options = compressed_options = 0
    arbitrary_groups = proper_options = 0
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        side_a = rank_parent.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = rank_parent.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        result, stats = orientation_result(
            instance, z0, side_a, side_b, subset_rank
        )
        transitions += result.transitions
        raw_options += stats.raw_options
        compressed_options += stats.compressed_options
        arbitrary_groups += stats.arbitrary_groups
        proper_options += stats.arbitrary_proper_options
        orientations.append({
            "case": name,
            "passed": result.passed,
            "first_failure": result.first_failure,
        })
    light = z0 | component.component
    light_dimension = light.bit_count()
    light_matching = (
        light_dimension <= COORDINATES
        and support_matching(light, instance.defects) is not None
    )
    passed = any(row["passed"] for row in orientations) or light_matching
    return {
        "passed": passed,
        "transitions": transitions,
        "raw_options": raw_options,
        "compressed_options": compressed_options,
        "arbitrary_groups": arbitrary_groups,
        "arbitrary_proper_options": proper_options,
        "certificate": None if passed else {
            "component": frozen.absolute_vertices(instance, component.component),
            "orientations": orientations,
            "lightlike": {
                "orthonormal_vectors": light_dimension,
                "dimension_passed": light_dimension <= COORDINATES,
                "allowed_mask_matching_passed": light_matching,
                "passed": light_matching,
            },
        },
    }


COUNTERS = (
    "z0_considered",
    "z0_matchable",
    "z0_arbitrary_passed",
    "z0_arbitrary_failed",
    "bipartite_components_checked",
    "bipartite_components_failed",
    "generic_orientations_checked",
    "arbitrary_positive_nonbip_groups",
    "arbitrary_proper_options",
    "raw_group_options",
    "compressed_group_options",
    "dominance_dp_transitions",
)


def empty_counts() -> dict[str, int]:
    return {name: 0 for name in COUNTERS}


def solve_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
    subset_rank: PositiveSubsetRank,
) -> dict:
    counts = empty_counts()
    rows = []
    for z0 in rank_parent.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        matchable = support_matching(z0, instance.defects) is not None
        if not matchable:
            rows.append({
                "Z0": frozen.absolute_vertices(instance, z0),
                "matchable": False,
                "first_failed_component": None,
            })
            continue
        counts["z0_matchable"] += 1
        failure = None
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                continue
            result = check_component(
                adj, instance, z0, component, inertia, zero_forcing, subset_rank
            )
            counts["bipartite_components_checked"] += 1
            counts["generic_orientations_checked"] += 2
            counts["arbitrary_positive_nonbip_groups"] += result["arbitrary_groups"]
            counts["arbitrary_proper_options"] += result["arbitrary_proper_options"]
            counts["raw_group_options"] += result["raw_options"]
            counts["compressed_group_options"] += result["compressed_options"]
            counts["dominance_dp_transitions"] += result["transitions"]
            if not result["passed"]:
                counts["bipartite_components_failed"] += 1
                failure = result["certificate"]
                break
        if failure is None:
            counts["z0_arbitrary_passed"] += 1
            return {"feasible": True, "counts": counts, "certificate_rows": None}
        counts["z0_arbitrary_failed"] += 1
        rows.append({
            "Z0": frozen.absolute_vertices(instance, z0),
            "matchable": True,
            "first_failed_component": failure,
        })
    return {"feasible": False, "counts": counts, "certificate_rows": rows}


def evaluate_graph(adj: Sequence[int]) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("arbitrary-subset production target must be K6-only")
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    subset_rank = PositiveSubsetRank(inertia, zero_forcing)
    totals = empty_counts()
    seeds_checked = 0
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        result = solve_seed(
            adj, instance, inertia, zero_forcing, subset_rank
        )
        for name in COUNTERS:
            totals[name] += result["counts"][name]
        if not result["feasible"]:
            return {
                "rejected": True,
                "seeds_checked": seeds_checked,
                "impossible_seeds": 1,
                **totals,
                "subset_rank_cache_entries": len(subset_rank.values),
                "subset_rank_cache_hits": subset_rank.hits,
                "first_impossible_seed": list(instance.seed),
                "certificate": {
                    "seed": list(instance.seed),
                    "outside": list(instance.outside),
                    "eligible_Z0": frozen.absolute_vertices(
                        instance, instance.eligible_z0_mask
                    ),
                    "allowed_defects": {
                        str(instance.outside[local]): [
                            instance.seed[coordinate]
                            for coordinate in vertices(mask)
                        ]
                        for local, mask in enumerate(instance.defects)
                    },
                    "choices": result["certificate_rows"],
                },
            }
    return {
        "rejected": False,
        "seeds_checked": seeds_checked,
        "impossible_seeds": 0,
        **totals,
        "subset_rank_cache_entries": len(subset_rank.values),
        "subset_rank_cache_hits": subset_rank.hits,
        "first_impossible_seed": None,
        "certificate": None,
    }


def evaluate_record(record: dict) -> dict:
    return {"index": record["index"], "decision": evaluate_graph(record["adjacency"])}


def synthetic_controls() -> dict:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    rank = PositiveSubsetRank(inertia, zero_forcing)
    c5 = (18, 5, 10, 20, 9)
    selected_p3 = 7
    if rank.solve(c5, selected_p3) != 2 or rank.solve(c5, 31) != 3:
        raise AssertionError("positive C5 subset-rank controls failed")
    # Full C5 uses four allowed coordinates and rank three; the nonadjacent
    # pair {0,2} has rank two but only one allowed coordinate.
    masks = (1, 2, 1, 4, 8)
    raw = [parent.Choice(0, 0, (), "empty")]
    for selected in range(1, 32):
        raw.append(parent.Choice(
            rank.solve(c5, selected),
            allowed_union(masks, selected),
            tuple(vertices(selected)),
            "full" if selected == 31 else "arbitrary_positive_nonbip",
        ))
    arbitrary = compress_choices("C5", raw, False)
    frozen_full = compress_choices(
        "C5",
        (
            parent.Choice(0, 0, (), "empty"),
            parent.Choice(3, 15, tuple(range(5)), "full"),
        ),
        False,
    )
    if not parent.check_groups((frozen_full,)).passed:
        raise AssertionError("frozen full-only C5 control failed")
    if parent.check_groups((arbitrary,)).passed:
        raise AssertionError("arbitrary-subset C5 Hall control should fail")
    overlap = parent.component_group(
        "overlap", (0, 1, 2), (1, 2, 4), 2, True
    )
    if not parent.check_groups((overlap,)).passed:
        raise AssertionError("one-span overlap control falsely added alternatives")
    return {
        "C5_selected_P3_rank": 2,
        "C5_full_rank": 3,
        "frozen_full_only_passed": True,
        "arbitrary_subset_failed": True,
        "one_span_overlap_alternatives_not_added": True,
        "one_vertex_extension_omitted": True,
    }


def verify_inputs() -> tuple[dict[str, str], list[dict], list[int]]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"arbitrary-subset dependency hash mismatch: {observed}")
    _, parent_records, _ = parent.verify_inputs()
    report = json.loads(
        (ROOT / "d6_k6_psd_z_hereditary_report.json").read_text()
    )
    indices = [
        row["index"] for row in report["graph_results"]
        if not row["decision"]["rejected"]
    ]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("pinned arbitrary-subset input boundary changed")
    by_index = {record["index"]: record for record in parent_records}
    return observed, [by_index[index] for index in indices], indices


def checkpoint_payload(source: str, indices: Sequence[int], completed: list[dict]) -> dict:
    return {
        "schema": CHECKPOINT_SCHEMA,
        "production_source_sha256": source,
        "input_indices_sha256": stable_hash(indices),
        "completed": completed,
    }


def run(
    output: Path,
    certificates: Path,
    checkpoint: Path,
    workers: int,
    checkpoint_every: int,
    resume: bool,
) -> dict:
    observed, records, indices = verify_inputs()
    source = sha256(Path(__file__).resolve())
    completed: list[dict] = []
    if resume and checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        if (
            saved.get("schema") != CHECKPOINT_SCHEMA
            or saved.get("production_source_sha256") != source
            or saved.get("input_indices_sha256") != stable_hash(indices)
        ):
            raise ValueError("arbitrary-subset checkpoint boundary mismatch")
        completed = saved["completed"]
        if [row["index"] for row in completed] != indices[:len(completed)]:
            raise ValueError("checkpoint is not an ordered input prefix")
    elif not resume:
        atomic_json(checkpoint, checkpoint_payload(source, indices, completed))
    started = time.perf_counter()
    for offset in range(len(completed), len(records), checkpoint_every):
        batch = records[offset:offset + checkpoint_every]
        if workers == 1:
            fresh = [evaluate_record(record) for record in batch]
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                fresh = list(pool.map(evaluate_record, batch, chunksize=1))
        completed.extend(fresh)
        atomic_json(checkpoint, checkpoint_payload(source, indices, completed))
    wall = time.perf_counter() - started
    rejected = [row for row in completed if row["decision"]["rejected"]]
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "status": "COMPLETE",
        "input_indices_sha256": stable_hash(indices),
        "production_source_sha256": source,
        "rejected_indices": [row["index"] for row in rejected],
        "certificates": [
            {"index": row["index"], "certificate": row["decision"]["certificate"]}
            for row in rejected
        ],
    }
    atomic_json(certificates, archive)
    graph_results = []
    for row in completed:
        decision = dict(row["decision"])
        decision.pop("certificate")
        graph_results.append({"index": row["index"], "decision": decision})
    numeric = ("seeds_checked", "impossible_seeds", *COUNTERS,
               "subset_rank_cache_entries", "subset_rank_cache_hits")
    totals = {
        name: sum(row["decision"][name] for row in graph_results)
        for name in numeric
    }
    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point control failed")
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact arbitrary-subset componentwise ZF/inertia/bipartite-PSD-Z "
            "rank Hall inside positive nonbipartite K6 side spans."
        ),
        "input_graphs": len(indices),
        "input_indices_sha256": stable_hash(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": [row["index"] for row in rejected],
        **totals,
        "graph_results": graph_results,
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "bytes": certificates.stat().st_size,
            "certificates": len(rejected),
        },
        "checkpoint": {
            "path": checkpoint.name,
            "sha256": sha256(checkpoint),
            "completed": len(completed),
        },
        "positive_control": {"passed": True, "K6_seeds": positive["seeds_checked"]},
        "synthetic_controls": synthetic_controls(),
        "sources": {**observed, Path(__file__).name: source},
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "span_quantifier": "one selected subset per original connected span",
            "one_vertex_extension": "deliberately omitted (zero pilot coverage)",
            "arithmetic": "exact integer/rational only",
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds_this_invocation": wall,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "d6_k6_arbitrary_subset_hall_report.json")
    parser.add_argument("--certificates", type=Path, default=ROOT / "d6_k6_arbitrary_subset_hall_certificates.json")
    parser.add_argument("--checkpoint", type=Path, default=ROOT / "d6_k6_arbitrary_subset_hall_checkpoint.json")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--checkpoint-every", type=int, default=50)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.checkpoint_every < 1:
        parser.error("workers and checkpoint-every must be positive")
    report = run(
        args.output.resolve(), args.certificates.resolve(), args.checkpoint.resolve(),
        args.workers, args.checkpoint_every, args.resume,
    )
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "graphs_surviving": report["graphs_surviving"],
        "rejected_indices": report["rejected_indices"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
