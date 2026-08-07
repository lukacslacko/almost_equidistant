#!/usr/bin/env python3
"""Exact pilot: arbitrary-subset rank Hall in positive nonbipartite K6 spans.

The frozen hereditary PSD--Z layer permits only the empty or full choice from
a positive-sign nonbipartite connected side span.  This pilot permits every
selected subset S.  It decomposes F[S] into connected induced components and
sums, across their mutually orthogonal Gram blocks, the componentwise maximum
of ordinary zero-forcing, positive-sign inertia for I+Adj, and the PSD--Z
|H|-1 bound when H is bipartite.

If a bipartite connected H in F[S] has a vertex v outside S for which
F[H union {v}] is connected and bipartite, then Gram(H) is a proper principal
submatrix of an irreducible PSD Z-matrix after signature switching.  It is
therefore positive definite and contributes |H|.

Exactly one subset is chosen from each original connected side span.  Subset
ranks from one span are alternatives and are never added to each other.
Allowed defect masks remain upper bounds and candidate nonedges may be unit.
All decisions are integer/exact; floating point is used only for timing.
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
import d6_k6_psd_z_hereditary as hereditary
import d6_k6_psd_zmatrix as psd_parent
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
SCHEMA = "d6-k6-arbitrary-subset-hall-pilot-v1"
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
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


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
    position = {vertex: local for local, vertex in enumerate(chosen)}
    return tuple(
        sum(
            1 << position[other]
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


def has_one_vertex_bipartite_extension(
    full_rows: Sequence[int], selected: int, component: int
) -> bool:
    """Whether H has a connected bipartite induced extension H+v outside S."""

    outside = ((1 << len(full_rows)) - 1) & ~selected
    while outside:
        bit = outside & -outside
        outside ^= bit
        vertex = bit.bit_length() - 1
        if not full_rows[vertex] & component:
            continue
        if is_bipartite(induced_rows(full_rows, component | bit)):
            return True
    return False


@dataclass(frozen=True)
class RankDetail:
    rank: int
    components: tuple[dict, ...]
    extension_upgrades: int


class SubsetRankCache:
    def __init__(self, inertia: InertiaCache, zero_forcing: ZeroForcingSolver):
        self.inertia = inertia
        self.zero_forcing = zero_forcing
        self.values: dict[tuple[tuple[int, ...], int, bool], RankDetail] = {}
        self.hits = 0

    def positive_rank(
        self, full_rows: Sequence[int], selected: int, use_extension: bool
    ) -> RankDetail:
        key = tuple(full_rows), selected, use_extension
        if key in self.values:
            self.hits += 1
            return self.values[key]
        total = 0
        details = []
        upgrades = 0
        for component in components_in_mask(full_rows, selected):
            local = induced_rows(full_rows, component)
            size = len(local)
            pattern = tuple(row | (1 << i) for i, row in enumerate(local))
            positive, negative, zero = self.inertia.solve(pattern)
            zf_number = self.zero_forcing.solve(local).number
            zf_rank = size - zf_number
            inertia_rank = size - negative
            bipartite = is_bipartite(local)
            psd_z_rank = size - 1 if bipartite else 0
            extension = bool(
                use_extension
                and bipartite
                and has_one_vertex_bipartite_extension(
                    full_rows, selected, component
                )
            )
            rank = max(
                zf_rank,
                inertia_rank,
                psd_z_rank,
                size if extension else 0,
            )
            upgrades += extension and rank > max(zf_rank, inertia_rank, psd_z_rank)
            total += rank
            details.append({
                "vertices": vertices(component),
                "size": size,
                "bipartite": bipartite,
                "ordinary_zero_forcing_number": zf_number,
                "ordinary_zero_forcing_rank_lower": zf_rank,
                "inertia_I_plus_Adj": [positive, negative, zero],
                "positive_sign_inertia_rank_lower": inertia_rank,
                "psd_zmatrix_rank_lower": psd_z_rank,
                "one_vertex_bipartite_extension": extension,
                "rank_lower": rank,
            })
        result = RankDetail(total, tuple(details), upgrades)
        self.values[key] = result
        return result


@dataclass(frozen=True)
class Choice:
    rank: int
    coordinates: int
    selected: tuple[int, ...]
    kind: str
    extension_upgrades: int = 0


@dataclass(frozen=True)
class Group:
    label: str
    choices: tuple[Choice, ...]
    raw_choices: tuple[Choice, ...]
    arbitrary: bool


@dataclass(frozen=True)
class HallResult:
    passed: bool
    transitions: int
    first_failure: dict | None


def compress_group(
    label: str, raw: Sequence[Choice], arbitrary: bool
) -> Group:
    best: dict[int, Choice] = {}
    for choice in raw:
        old = best.get(choice.coordinates)
        if old is None or choice.rank > old.rank:
            best[choice.coordinates] = choice
    choices = tuple(
        best[mask]
        for mask in sorted(best, key=lambda item: (item.bit_count(), item))
    )
    return Group(label, choices, tuple(raw), arbitrary)


def allowed_union(masks: Sequence[int], selected: int) -> int:
    answer = 0
    for vertex in vertices(selected):
        answer |= masks[vertex]
    return answer


def make_side_groups(
    instance: K6LorentzInstance,
    label: str,
    side,
    rank_cache: SubsetRankCache,
    use_extension: bool,
) -> tuple[Group, ...]:
    local_by_absolute = {
        absolute: local for local, absolute in enumerate(instance.outside)
    }
    answer = []
    for number, detail in enumerate(side.components):
        absolute = tuple(detail["vertices"])
        masks = tuple(
            instance.defects[local_by_absolute[vertex]] for vertex in absolute
        )
        size = len(absolute)
        applicable = bool(detail["psd_zmatrix_applicable"])
        full_rank = int(detail["componentwise_fused_rank_lower"])
        raw = [Choice(0, 0, (), "empty")]
        for selected in range(1, 1 << size):
            full = selected == (1 << size) - 1
            if applicable:
                rank = full_rank if full else selected.bit_count()
                upgrades = 0
                kind = "full" if full else "hereditary_proper"
            else:
                rank_detail = rank_cache.positive_rank(
                    detail["adjacency_rows"], selected, use_extension
                )
                rank = rank_detail.rank
                upgrades = rank_detail.extension_upgrades
                kind = "full" if full else "arbitrary_positive_nonbip"
                if full and rank != full_rank:
                    raise AssertionError(
                        "arbitrary full-subset rank differs from frozen bound"
                    )
            raw.append(Choice(
                rank,
                allowed_union(masks, selected),
                tuple(absolute[i] for i in vertices(selected)),
                kind,
                upgrades,
            ))
        answer.append(compress_group(
            f"{label}:F{number}", raw, not applicable
        ))
    return tuple(answer)


def z0_groups(instance: K6LorentzInstance, z0: int) -> tuple[Group, ...]:
    return tuple(
        compress_group(
            f"Z0:{instance.outside[local]}",
            (
                Choice(0, 0, (), "empty"),
                Choice(
                    1,
                    instance.defects[local],
                    (instance.outside[local],),
                    "Z0",
                ),
            ),
            False,
        )
        for local in vertices(z0)
    )


def check_groups(groups: Sequence[Group]) -> HallResult:
    states: dict[int, tuple[int, tuple[tuple[str, Choice], ...]]] = {0: (0, ())}
    transitions = 0
    for group in groups:
        new: dict[int, tuple[int, tuple[tuple[str, Choice], ...]]] = {}
        for coordinates, (rank, witness) in states.items():
            for choice in group.choices:
                transitions += 1
                union = coordinates | choice.coordinates
                total = rank + choice.rank
                selected = witness
                if choice.selected:
                    selected += ((group.label, choice),)
                if total > union.bit_count():
                    return HallResult(False, transitions, {
                        "rank_lower": total,
                        "coordinate_capacity": union.bit_count(),
                        "allowed_coordinates": vertices(union),
                        "selected_groups": [
                            {
                                "group": name,
                                "kind": item.kind,
                                "rank_lower": item.rank,
                                "selected_vertices": list(item.selected),
                                "allowed_coordinates": vertices(item.coordinates),
                                "extension_upgrades": item.extension_upgrades,
                            }
                            for name, item in selected
                        ],
                    })
                old = new.get(union)
                if old is None or total > old[0]:
                    new[union] = (total, selected)
        states = new
    return HallResult(True, transitions, None)


def brute_check_groups(groups: Sequence[Group]) -> bool:
    """Independent raw Cartesian-product check, with no union dominance."""

    def visit(position: int, rank: int, coordinates: int) -> bool:
        if position == len(groups):
            return True
        for choice in groups[position].raw_choices:
            union = coordinates | choice.coordinates
            total = rank + choice.rank
            if total > union.bit_count():
                return False
            if not visit(position + 1, total, union):
                return False
        return True

    return visit(0, 0, 0)


@dataclass
class CrossCheck:
    systems: int = 0


def orientation_result(
    instance: K6LorentzInstance,
    z0: int,
    side_a,
    side_b,
    rank_cache: SubsetRankCache,
    use_extension: bool,
    cross_check: CrossCheck | None,
) -> tuple[HallResult, tuple[Group, ...]]:
    groups = make_side_groups(
        instance, "A", side_a, rank_cache, use_extension
    ) + make_side_groups(
        instance, "B", side_b, rank_cache, use_extension
    ) + z0_groups(instance, z0)
    result = check_groups(groups)
    if cross_check is not None:
        brute = brute_check_groups(groups)
        cross_check.systems += 1
        if brute != result.passed:
            raise AssertionError("dominance DP and raw Cartesian Hall differ")
    return result, groups


def check_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
    rank_cache: SubsetRankCache,
    use_extension: bool,
    cross_check: CrossCheck | None,
) -> dict:
    graph_a, absolute_a = induced_required_graph(adj, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adj, instance, component.side_b)
    orientations = []
    transitions = 0
    raw_options = 0
    compressed_options = 0
    arbitrary_groups = 0
    arbitrary_options = 0
    extension_options = 0
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        side_a = psd_parent.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = psd_parent.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        new, groups = orientation_result(
            instance,
            z0,
            side_a,
            side_b,
            rank_cache,
            use_extension,
            cross_check,
        )
        old, _ = hereditary.orientation_result(instance, z0, side_a, side_b)
        if not old.passed and new.passed:
            raise AssertionError("arbitrary-subset Hall revived frozen failure")
        transitions += new.transitions
        raw_options += sum(len(group.raw_choices) for group in groups)
        compressed_options += sum(len(group.choices) for group in groups)
        arbitrary_groups += sum(group.arbitrary for group in groups)
        arbitrary_options += sum(
            len(group.raw_choices) - 2 for group in groups if group.arbitrary
        )
        extension_options += sum(
            choice.extension_upgrades > 0
            for group in groups
            for choice in group.raw_choices
            if group.arbitrary
        )
        orientations.append({
            "case": name,
            "frozen_passed": old.passed,
            "arbitrary_subset_passed": new.passed,
            "strict_new_failure": old.passed and not new.passed,
            "first_failure": new.first_failure,
        })
    light = z0 | component.component
    light_dimension = light.bit_count()
    light_passed = (
        light_dimension <= COORDINATES
        and support_matching(light, instance.defects) is not None
    )
    old_passed = any(row["frozen_passed"] for row in orientations) or light_passed
    new_passed = any(row["arbitrary_subset_passed"] for row in orientations) or light_passed
    if not old_passed and new_passed:
        raise AssertionError("arbitrary-subset component revived frozen failure")
    return {
        "passed": new_passed,
        "frozen_passed": old_passed,
        "strict_new_failure": old_passed and not new_passed,
        "orientations": orientations,
        "lightlike_passed": light_passed,
        "component": frozen.absolute_vertices(instance, component.component),
        "transitions": transitions,
        "raw_options": raw_options,
        "compressed_options": compressed_options,
        "arbitrary_groups": arbitrary_groups,
        "arbitrary_options": arbitrary_options,
        "extension_options": extension_options,
    }


COUNTERS = (
    "z0_considered",
    "z0_matchable",
    "z0_arbitrary_passed",
    "z0_arbitrary_failed",
    "components_checked",
    "components_failed",
    "components_strict_new_failed",
    "orientations_checked",
    "orientations_strict_new_failed",
    "arbitrary_groups",
    "arbitrary_proper_options",
    "extension_upgraded_options",
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
    rank_cache: SubsetRankCache,
    use_extension: bool,
    cross_check: CrossCheck | None,
) -> dict:
    counts = empty_counts()
    choices = []
    for z0 in psd_parent.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            choices.append({
                "Z0": frozen.absolute_vertices(instance, z0),
                "matchable": False,
                "first_failure": None,
            })
            continue
        counts["z0_matchable"] += 1
        failure = None
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                continue
            result = check_component(
                adj,
                instance,
                z0,
                component,
                inertia,
                zero_forcing,
                rank_cache,
                use_extension,
                cross_check,
            )
            counts["components_checked"] += 1
            counts["components_failed"] += not result["passed"]
            counts["components_strict_new_failed"] += result["strict_new_failure"]
            counts["orientations_checked"] += 2
            counts["orientations_strict_new_failed"] += sum(
                row["strict_new_failure"] for row in result["orientations"]
            )
            counts["arbitrary_groups"] += result["arbitrary_groups"]
            counts["arbitrary_proper_options"] += result["arbitrary_options"]
            counts["extension_upgraded_options"] += result["extension_options"]
            counts["raw_group_options"] += result["raw_options"]
            counts["compressed_group_options"] += result["compressed_options"]
            counts["dominance_dp_transitions"] += result["transitions"]
            if not result["passed"]:
                failure = result
                break
        if failure is None:
            counts["z0_arbitrary_passed"] += 1
            return {"feasible": True, "counts": counts, "choices": None}
        counts["z0_arbitrary_failed"] += 1
        choices.append({
            "Z0": frozen.absolute_vertices(instance, z0),
            "matchable": True,
            "first_failure": failure,
        })
    return {"feasible": False, "counts": counts, "choices": choices}


def add_counts(target: dict[str, int], source: dict[str, int]) -> None:
    for name in COUNTERS:
        target[name] += source[name]


def evaluate_graph(
    adj: Sequence[int],
    use_extension: bool = True,
    cross_check: CrossCheck | None = None,
) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("arbitrary-subset K6 target contains K7")
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    rank_cache = SubsetRankCache(inertia, zero_forcing)
    totals = empty_counts()
    seeds_checked = 0
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        result = solve_seed(
            adj,
            instance,
            inertia,
            zero_forcing,
            rank_cache,
            use_extension,
            cross_check,
        )
        add_counts(totals, result["counts"])
        if not result["feasible"]:
            return {
                "rejected": True,
                "seeds_checked": seeds_checked,
                "impossible_seeds": 1,
                **totals,
                "subset_rank_cache_entries": len(rank_cache.values),
                "subset_rank_cache_hits": rank_cache.hits,
                "first_impossible_seed": list(instance.seed),
                "certificate": {
                    "seed": list(instance.seed),
                    "choices": result["choices"],
                },
            }
    return {
        "rejected": False,
        "seeds_checked": seeds_checked,
        "impossible_seeds": 0,
        **totals,
        "subset_rank_cache_entries": len(rank_cache.values),
        "subset_rank_cache_hits": rank_cache.hits,
        "first_impossible_seed": None,
        "certificate": None,
    }


def evaluate_record(record: dict) -> dict:
    return {
        "index": record["index"],
        "decision": evaluate_graph(record["adjacency"]),
    }


def evaluate_record_without_extension(record: dict) -> dict:
    return {
        "index": record["index"],
        "decision": evaluate_graph(record["adjacency"], use_extension=False),
    }


def synthetic_controls() -> dict:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    cache = SubsetRankCache(inertia, zero_forcing)
    c5 = (18, 5, 10, 20, 9)
    path3 = (1 << 0) | (1 << 1) | (1 << 2)
    without = cache.positive_rank(c5, path3, False)
    with_extension = cache.positive_rank(c5, path3, True)
    if without.rank != 2 or with_extension.rank != 3:
        raise AssertionError("C5/P3 one-vertex extension control failed")

    absolute = (0, 1, 2, 3, 4)
    masks = (3, 1, 3, 4, 8)
    detail = {
        "vertices": list(absolute),
        "adjacency_rows": list(c5),
        "componentwise_fused_rank_lower": cache.positive_rank(
            c5, (1 << 5) - 1, True
        ).rank,
        "psd_zmatrix_applicable": False,
    }
    fake_side = type("Side", (), {"components": (detail,)})()
    fake_instance = type("Instance", (), {
        "outside": absolute,
        "defects": masks,
    })()
    group_without = make_side_groups(
        fake_instance, "A", fake_side, cache, False
    )[0]
    group_with = make_side_groups(
        fake_instance, "A", fake_side, cache, True
    )[0]
    old = compress_group(
        "old",
        (
            Choice(0, 0, (), "empty"),
            Choice(detail["componentwise_fused_rank_lower"], 15, absolute, "full"),
        ),
        False,
    )
    if not check_groups((old,)).passed:
        raise AssertionError("synthetic frozen full-only control failed")
    if not check_groups((group_without,)).passed:
        raise AssertionError("synthetic no-extension arbitrary control failed")
    if check_groups((group_with,)).passed:
        raise AssertionError("synthetic extension Hall control should fail")
    if brute_check_groups((group_with,)):
        raise AssertionError("synthetic brute Hall missed extension failure")
    return {
        "C5_selected_P3_rank_without_extension": without.rank,
        "C5_selected_P3_rank_with_extension": with_extension.rank,
        "frozen_full_only_passed": True,
        "arbitrary_without_extension_passed": True,
        "arbitrary_with_extension_failed": True,
        "dominance_and_brute_agree": True,
    }


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"arbitrary-subset dependency mismatch: {observed}")
    _, parent_records, _ = hereditary.verify_inputs()
    report = json.loads(
        (ROOT / "d6_k6_psd_z_hereditary_report.json").read_text()
    )
    indices = [
        row["index"] for row in report["graph_results"]
        if not row["decision"]["rejected"]
    ]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("pinned 831 arbitrary-subset boundary changed")
    by_index = {record["index"]: record for record in parent_records}
    return [by_index[index] for index in indices], indices, observed


def cross_check_sample(
    records: Sequence[dict], count: int, use_extension: bool
) -> dict:
    checked = CrossCheck()
    graph_results = []
    started = time.perf_counter()
    for record in records[:count]:
        result = evaluate_graph(
            record["adjacency"],
            use_extension=use_extension,
            cross_check=checked,
        )
        graph_results.append({"index": record["index"], "rejected": result["rejected"]})
    return {
        "graphs": len(graph_results),
        "systems": checked.systems,
        "all_dominance_decisions_equal_raw_cartesian": True,
        "graph_results": graph_results,
        "wall_seconds": time.perf_counter() - started,
    }


def run(
    output: Path,
    workers: int,
    cross_check_graphs: int,
    use_extension: bool,
) -> dict:
    controls = synthetic_controls()
    records, indices, observed = load_input()
    cross = cross_check_sample(records, cross_check_graphs, use_extension)
    evaluator = evaluate_record if use_extension else evaluate_record_without_extension
    started = time.perf_counter()
    if workers == 1:
        graph_results = [evaluator(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            graph_results = list(pool.map(evaluator, records, chunksize=1))
    wall = time.perf_counter() - started
    rejected = [row for row in graph_results if row["decision"]["rejected"]]
    numeric = sorted({
        key
        for row in graph_results
        for key, value in row["decision"].items()
        if isinstance(value, int) and not isinstance(value, bool)
    })
    totals = {
        key: sum(row["decision"].get(key, 0) for row in graph_results)
        for key in numeric
    }
    positive = evaluate_graph(
        lower_bound_18_graph(), use_extension=use_extension
    )
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point control failed")
    report = {
        "schema": SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact pilot of arbitrary-subset rank Hall inside positive-sign "
            "nonbipartite connected K6 side spans."
        ),
        "one_vertex_bipartite_extension_enabled": use_extension,
        "input_graphs": len(indices),
        "input_indices_sha256": stable_hash(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": [row["index"] for row in rejected],
        **totals,
        "graph_results": graph_results,
        "synthetic_controls": controls,
        "brute_force_cross_check": cross,
        "positive_control": {
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
        },
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "span_quantifier": "one selected subset per original connected span",
            "arithmetic": "exact integer decisions only",
            "status": "pilot, not production theorem credit",
        },
        "sources": {**observed, Path(__file__).name: sha256(Path(__file__).resolve())},
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": wall,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_arbitrary_subset_hall_pilot_report.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--cross-check-graphs", type=int, default=8)
    parser.add_argument("--no-extension", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.cross_check_graphs < 1:
        parser.error("workers and cross-check-graphs must be positive")
    report = run(
        args.output.resolve(),
        args.workers,
        args.cross_check_graphs,
        not args.no_extension,
    )
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "graphs_surviving": report["graphs_surviving"],
        "rejected_indices": report["rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
