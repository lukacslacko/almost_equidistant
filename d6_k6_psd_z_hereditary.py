#!/usr/bin/env python3
"""Exact hereditary PSD--Z support-Hall layer on the pinned K6 residue.

Every proper principal submatrix of a connected irreducible PSD Z-matrix is
positive definite.  Thus, inside each PSD--Z-applicable connected side-graph
block, every proper subset of defect vectors is independent.  For one
generic Lorentz orientation we select at most one subset from each connected
side span, plus singleton Z0 blocks, and apply all coordinate-support Hall
inequalities.  The lightlike Lorentz case remains a separate alternative.

All decisions are exact integer computations.  Candidate nonedges are
unconstrained and allowed defect masks are upper bounds only.
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
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k6_fused_side_rank as frozen
import d6_k6_psd_zmatrix as parent
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
REPORT_SCHEMA = "d6-k6-psd-z-hereditary-v1"
CERTIFICATE_SCHEMA = "d6-k6-psd-z-hereditary-certificates-v1"
CHECKPOINT_SCHEMA = "d6-k6-psd-z-hereditary-checkpoint-v1"
EXPECTED_INPUT = 861
EXPECTED_INPUT_SHA256 = (
    "09ebce17d2b72fa6514fc6a8a938376626d373161d2e1b4dd8d31c5e00c18db5"
)
EXPECTED = {
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "d6_k6_psd_zmatrix_report.json": (
        "abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30"
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


@dataclass(frozen=True)
class Choice:
    rank: int
    coordinates: int
    selected: tuple[int, ...]
    kind: str


@dataclass(frozen=True)
class Group:
    label: str
    choices: tuple[Choice, ...]
    applicable: bool


@dataclass(frozen=True)
class HallResult:
    passed: bool
    transitions: int
    first_failure: dict | None


def component_group(
    label: str,
    absolute_vertices: Sequence[int],
    coordinate_masks: Sequence[int],
    full_rank: int,
    applicable: bool,
) -> Group:
    """Construct alternatives for one span; alternatives never add together."""

    size = len(absolute_vertices)
    if len(coordinate_masks) != size or not 0 <= full_rank <= size:
        raise ValueError("invalid hereditary component group")
    if applicable and full_rank < size - 1:
        raise ValueError("applicable PSD Z component lost its n-1 rank bound")
    raw = [Choice(0, 0, (), "empty")]
    subsets = range(1, 1 << size) if applicable else ((1 << size) - 1,)
    for encoded in subsets:
        positions = tuple(i for i in range(size) if encoded & (1 << i))
        selected = tuple(absolute_vertices[i] for i in positions)
        coordinates = 0
        for i in positions:
            coordinates |= coordinate_masks[i]
        full = encoded == (1 << size) - 1
        raw.append(Choice(
            full_rank if full else len(positions),
            coordinates,
            selected,
            "full" if full else "proper",
        ))
    # Same coordinate union: larger rank dominates for every future group.
    best: dict[int, Choice] = {}
    for choice in raw:
        old = best.get(choice.coordinates)
        if old is None or choice.rank > old.rank:
            best[choice.coordinates] = choice
    return Group(
        label,
        tuple(best[mask] for mask in sorted(best, key=lambda x: (x.bit_count(), x))),
        applicable,
    )


def singleton_group(label: str, absolute: int, coordinates: int) -> Group:
    return Group(
        label,
        (Choice(0, 0, (), "empty"), Choice(1, coordinates, (absolute,), "Z0")),
        True,
    )


def check_groups(groups: Sequence[Group]) -> HallResult:
    """Exact 64-state dominance DP over one choice from each orthogonal span."""

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
                            }
                            for name, item in selected
                        ],
                    })
                old = new.get(union)
                if old is None or total > old[0]:
                    new[union] = (total, selected)
        states = new
    return HallResult(True, transitions, None)


def side_groups(
    instance: K6LorentzInstance, label: str, side
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
        answer.append(component_group(
            f"{label}:F{number}",
            absolute,
            masks,
            detail["componentwise_fused_rank_lower"],
            detail["psd_zmatrix_applicable"],
        ))
    return tuple(answer)


def orientation_result(
    instance: K6LorentzInstance, z0: int, side_a, side_b
) -> tuple[HallResult, tuple[Group, ...]]:
    groups = side_groups(instance, "A", side_a) + side_groups(instance, "B", side_b)
    groups += tuple(
        singleton_group(
            f"Z0:{instance.outside[local]}",
            instance.outside[local],
            instance.defects[local],
        )
        for local in vertices(z0)
    )
    return check_groups(groups), groups


def check_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> dict:
    graph_a, absolute_a = induced_required_graph(adj, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adj, instance, component.side_b)
    orientations = []
    transitions = 0
    applicable_groups = 0
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        side_a = parent.side_rank_lower(graph_a, absolute_a, sign_a, inertia, zero_forcing)
        side_b = parent.side_rank_lower(graph_b, absolute_b, sign_b, inertia, zero_forcing)
        result, groups = orientation_result(instance, z0, side_a, side_b)
        transitions += result.transitions
        applicable_groups += sum(group.applicable for group in groups if not group.label.startswith("Z0:"))
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
        "applicable_groups": applicable_groups,
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
    "z0_hereditary_passed",
    "z0_hereditary_failed",
    "bipartite_components_checked",
    "bipartite_components_failed",
    "generic_orientations_checked",
    "psd_zmatrix_applicable_groups",
    "dominance_dp_transitions",
)


def empty_counts() -> dict[str, int]:
    return {name: 0 for name in COUNTERS}


def solve_seed(adj: Sequence[int], instance: K6LorentzInstance) -> dict:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    counts = empty_counts()
    certificate_rows = []
    for z0 in parent.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        matchable = support_matching(z0, instance.defects) is not None
        if not matchable:
            certificate_rows.append({
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
                adj, instance, z0, component, inertia, zero_forcing
            )
            counts["bipartite_components_checked"] += 1
            counts["generic_orientations_checked"] += 2
            counts["psd_zmatrix_applicable_groups"] += result["applicable_groups"]
            counts["dominance_dp_transitions"] += result["transitions"]
            if not result["passed"]:
                counts["bipartite_components_failed"] += 1
                failure = result["certificate"]
                break
        if failure is None:
            counts["z0_hereditary_passed"] += 1
            return {"feasible": True, "counts": counts, "certificate_rows": None}
        counts["z0_hereditary_failed"] += 1
        certificate_rows.append({
            "Z0": frozen.absolute_vertices(instance, z0),
            "matchable": True,
            "first_failed_component": failure,
        })
    return {"feasible": False, "counts": counts, "certificate_rows": certificate_rows}


def evaluate_graph(adj: Sequence[int]) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("hereditary production target must be K6-only")
    totals = empty_counts()
    seeds_checked = 0
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        result = solve_seed(adj, instance)
        for name in COUNTERS:
            totals[name] += result["counts"][name]
        if not result["feasible"]:
            return {
                "rejected": True,
                "seeds_checked": seeds_checked,
                "impossible_seeds": 1,
                **totals,
                "first_impossible_seed": list(instance.seed),
                "certificate": {
                    "seed": list(instance.seed),
                    "outside": list(instance.outside),
                    "eligible_Z0": frozen.absolute_vertices(instance, instance.eligible_z0_mask),
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
        "first_impossible_seed": None,
        "certificate": None,
    }


def evaluate_record(record: dict) -> dict:
    return {"index": record["index"], "decision": evaluate_graph(record["adjacency"])}


def synthetic_controls() -> dict:
    # One group with two overlapping proper subsets.  Each subset is an
    # alternative choice, so the system passes; adding their ranks (4 > 3)
    # would be the deliberately forbidden false inference.
    overlap = component_group(
        "overlap",
        (0, 1, 2),
        (1, 2, 4),
        2,
        True,
    )
    overlap_result = check_groups((overlap,))
    if not overlap_result.passed:
        raise AssertionError("one-span overlapping-subset control falsely failed")
    separate = check_groups((overlap, overlap))
    if separate.passed:
        raise AssertionError("distinct-span additive control should fail")
    strict = component_group(
        "strict",
        (0, 1, 2),
        (1, 1, 2),
        2,
        True,
    )
    strict_result = check_groups((strict,))
    if strict_result.passed:
        raise AssertionError("proper-subset Hall control should fail")
    return {
        "one_span_overlapping_subsets_passed": overlap_result.passed,
        "forbidden_naive_overlap_rank": 4,
        "forbidden_naive_overlap_capacity": 3,
        "two_distinct_spans_failed": not separate.passed,
        "proper_subset_strict_failure": strict_result.first_failure,
    }


def verify_inputs() -> tuple[dict[str, str], list[dict], list[int]]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"hereditary dependency hash mismatch: {observed}")
    _, records, _ = parent.verify_inputs()
    report = json.loads((ROOT / "d6_k6_psd_zmatrix_report.json").read_text())
    indices = [
        item["index"] for item in report["graph_results"]
        if not item["decision"]["rejected"]
    ]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("pinned hereditary input boundary changed")
    by_index = {record["index"]: record for record in records}
    return observed, [by_index[index] for index in indices], indices


def checkpoint_payload(source_hash: str, indices: Sequence[int], completed: list[dict]) -> dict:
    return {
        "schema": CHECKPOINT_SCHEMA,
        "production_source_sha256": source_hash,
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
    source_hash = sha256(Path(__file__).resolve())
    completed: list[dict] = []
    if resume and checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        if (
            saved.get("schema") != CHECKPOINT_SCHEMA
            or saved.get("production_source_sha256") != source_hash
            or saved.get("input_indices_sha256") != stable_hash(indices)
        ):
            raise ValueError("checkpoint source/input boundary mismatch")
        completed = saved["completed"]
        if [row["index"] for row in completed] != indices[:len(completed)]:
            raise ValueError("checkpoint is not an ordered input prefix")
    elif not resume:
        atomic_json(checkpoint, checkpoint_payload(source_hash, indices, completed))

    started = time.perf_counter()
    for offset in range(len(completed), len(records), checkpoint_every):
        batch = records[offset:offset + checkpoint_every]
        if workers == 1:
            fresh = [evaluate_record(record) for record in batch]
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                fresh = list(pool.map(evaluate_record, batch, chunksize=1))
        completed.extend(fresh)
        atomic_json(checkpoint, checkpoint_payload(source_hash, indices, completed))
    wall = time.perf_counter() - started

    rejected = [row for row in completed if row["decision"]["rejected"]]
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "status": "COMPLETE",
        "input_indices_sha256": stable_hash(indices),
        "production_source_sha256": source_hash,
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
    numeric = (
        "seeds_checked", "impossible_seeds", *COUNTERS,
    )
    totals = {
        name: sum(row["decision"][name] for row in graph_results) for name in numeric
    }
    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point control failed")
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": "Exact hereditary principal-submatrix PSD-Z support Hall.",
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
        "positive_control": {
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
        },
        "synthetic_controls": synthetic_controls(),
        "sources": {**observed, Path(__file__).name: source_hash},
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "span_quantifier": "one selected subset per connected orthogonal span",
            "lightlike": "separate alternative",
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
    parser.add_argument("--output", type=Path, default=ROOT / "d6_k6_psd_z_hereditary_report.json")
    parser.add_argument("--certificates", type=Path, default=ROOT / "d6_k6_psd_z_hereditary_certificates.json")
    parser.add_argument("--checkpoint", type=Path, default=ROOT / "d6_k6_psd_z_hereditary_checkpoint.json")
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
