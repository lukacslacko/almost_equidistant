#!/usr/bin/env python3
"""Read-only exact pilot for promoting an empty K6 defect to an actual K7.

Fix a required K6 ``Q`` in a putative realization.  If an outside point ``x``
has actual K6 defect vector zero, then all six pairs ``xq`` are unit.  Adding
those six (possibly omitted) required edges produces a sound supergraph with
the required K7 ``Q union {x}``.  This pilot subjects that *particular* K7
seed to generic, archive-free exact K7 graph/rank, labeled-support,
propagated-link, and sparse-value filters.

An eliminated augmentation proves only ``u_Q(x) != 0``.  The original graph
is never called rejected: the all-nonempty branch and every surviving
one-/two-empty branch remain outside this pilot.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import d6_k6_arbitrary_subset_hall as parent
import d6_k7_propagated_link_caps as link_caps
import d6_k7_rank_reference as k7
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation
from d6_k6_lorentz_reference import (
    COORDINATES,
    add_edge,
    clique_masks as k6_clique_masks,
    validate_graph as validate_k6_graph,
    vertices,
)


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_arbitrary_subset_hall_report.json"
EXPECTED_PARENT_REPORT_SHA256 = (
    "bcc4baecf28a95443a5b13e87a796cbba0a97cd42fa76573119bfd81ef3735fc"
)
EXPECTED_RESIDUE = 822
EXPECTED_RESIDUE_SHA256 = (
    "cf94aac41eba28904ee6254f9c50996e73d38dcc15316fb6a411b29ff5cd05d8"
)
SAMPLE_STRATA = 4
SAMPLE_PER_STRATUM = 5
STATUS_GRAPH = "GENERIC_GRAPH_RANK_ELIMINATED"
STATUS_SUPPORT = "LABELED_SUPPORT_ELIMINATED"
STATUS_LINK = "PROPAGATED_LINK_CAP_ELIMINATED"
STATUS_SPARSE = "SPARSE_VALUE_ELIMINATED"
STATUS_SURVIVOR = "SURVIVOR"
STATUSES = (STATUS_GRAPH, STATUS_SUPPORT, STATUS_LINK, STATUS_SPARSE, STATUS_SURVIVOR)
QOS_CLASS_USER_INITIATED = 0x19


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


def worker_init() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    if sys.platform != "darwin":
        return
    libc = ctypes.CDLL(None)
    setter = libc.pthread_set_qos_class_self_np
    setter.argtypes = (ctypes.c_uint, ctypes.c_int)
    setter.restype = ctypes.c_int
    error = int(setter(QOS_CLASS_USER_INITIATED, 0))
    if error:
        raise OSError(error, "pthread_set_qos_class_self_np failed")


def load_residue() -> tuple[list[dict], list[int], dict]:
    if sha256(PARENT_REPORT) != EXPECTED_PARENT_REPORT_SHA256:
        raise ValueError("K6 arbitrary-subset parent report hash changed")
    _, records, ordered_831 = parent.verify_inputs()
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    if report.get("status") != "COMPLETE" or report.get("input_graphs") != 831:
        raise ValueError("K6 arbitrary-subset parent is not complete")
    rejected = {int(index) for index in report["rejected_indices"]}
    residue_indices = [index for index in ordered_831 if index not in rejected]
    if (
        len(residue_indices) != EXPECTED_RESIDUE
        or stable_hash(residue_indices) != EXPECTED_RESIDUE_SHA256
    ):
        raise ValueError("ordered 822-graph K6 residue changed")
    by_index = {int(record["index"]): record for record in records}
    residue = [by_index[index] for index in residue_indices]
    return residue, residue_indices, report


def count_k6(adjacency: Sequence[int]) -> int:
    return sum(1 for _ in k6_clique_masks(adjacency, COORDINATES))


def stratified_selection(records: Sequence[dict]) -> list[dict]:
    """Choose five stable-hash records from each K6-count quartile."""

    scored = []
    for record in records:
        count = count_k6(record["adjacency"])
        priority = hashlib.sha256(str(record["index"]).encode("ascii")).digest()
        scored.append((count, priority, record))
    scored.sort(key=lambda row: (row[0], row[1]))
    selected = []
    for stratum in range(SAMPLE_STRATA):
        lo = len(scored) * stratum // SAMPLE_STRATA
        hi = len(scored) * (stratum + 1) // SAMPLE_STRATA
        block = sorted(scored[lo:hi], key=lambda row: row[1])
        for count, _, record in block[:SAMPLE_PER_STRATUM]:
            selected.append({**record, "k6_count": count, "stratum": stratum})
    if len(selected) != SAMPLE_STRATA * SAMPLE_PER_STRATUM:
        raise AssertionError("stratified selection size changed")
    return selected


def augment_with_virtual_k7(
    adjacency: Sequence[int], seed: Sequence[int], apex: int
) -> tuple[tuple[int, ...], int]:
    """Add exactly the branch-forced apex-to-K6 edges."""

    if len(seed) != 6 or apex in seed:
        raise ValueError("virtual promotion needs a K6 and an outside apex")
    augmented = list(map(int, adjacency))
    added = 0
    for q in seed:
        if not (augmented[apex] & (1 << q)):
            add_edge(augmented, apex, q)
            added += 1
    validate_k6_graph(augmented)
    seed_mask = sum(1 << vertex for vertex in (*seed, apex))
    if next(k7.clique_masks(augmented, 7, ), None) is None:
        raise AssertionError("promotion failed to create a K7")
    # Check the intended seed directly; the graph may contain other new K7s.
    intended = list(vertices(seed_mask))
    if len(intended) != 7 or any(
        not (augmented[u] & (1 << v)) for u, v in combinations(intended, 2)
    ):
        raise AssertionError("intended virtual K7 is not a clique")
    return tuple(augmented), added


def baseline_failure_names(analysis: k7.CoverAnalysis) -> tuple[str, ...]:
    tests = (
        ("direct_cap", analysis.direct_cap_failure is not None),
        ("zero_factor_support", analysis.support_failed),
        ("subspace_K", analysis.subspace_k_failed),
        ("component_B", analysis.component_b_failed),
        ("PD_clique", analysis.pd_clique_failed),
        ("perpendicular_degree", analysis.perpendicular_degree_failed),
        ("basis_kernel", analysis.basis_kernel_failed),
        ("saturating_mask", analysis.saturating_mask_failed),
    )
    return tuple(name for name, failed in tests if failed)


def analyze_virtual_seed(
    adjacency: Sequence[int],
    seed: Sequence[int],
    apex: int,
    support_solver: k7.SupportSolver,
    zero_forcing: k7.ZeroForcingSolver,
    clique_solver: k7.CliqueStructureSolver,
) -> dict:
    """Apply the exact archive-free conjunction to one promoted K7 seed."""

    augmented, added_edges = augment_with_virtual_k7(adjacency, seed, apex)
    k7.validate_graph(augmented)
    seed_mask = sum(1 << vertex for vertex in (*seed, apex))
    _, outside, defects, ladj, eligible = k7.seed_instance(augmented, seed_mask)
    all_covers = k7.eligible_covers(ladj, eligible, cap=7)
    covers = [cover for cover in all_covers if cover.bit_count() <= 3]
    total_term_rank = k7.matching_size(defects)

    baseline_failures: Counter[str] = Counter()
    propagation_failures: Counter[str] = Counter()
    sparse_failures: Counter[str] = Counter()
    link_failure_sizes: Counter[str] = Counter()
    counters: Counter[str] = Counter(
        eligible_covers=len(all_covers),
        small_covers=len(covers),
        direct_cap_covers=len(all_covers) - len(covers),
    )
    survivor_witness = None

    for zmask in covers:
        counters["covers_analyzed"] += 1
        baseline = k7.analyze_cover(
            augmented,
            outside,
            defects,
            zmask,
            support_solver,
            zero_forcing,
            total_term_rank,
            clique_solver,
        )
        if baseline.enhanced_joint_failed:
            names = baseline_failure_names(baseline)
            if not names:
                raise AssertionError("enhanced baseline failed without a reason")
            baseline_failures.update(names)
            continue
        counters["baseline_passing_covers"] += 1

        zvertices = tuple(k7.bits(zmask))
        nvertices = tuple(
            vertex for vertex in range(len(outside)) if not (zmask & (1 << vertex))
        )
        graph_n = k7.induced_graph(
            augmented, [outside[vertex] for vertex in nvertices]
        )
        z_allowed = tuple(defects[vertex] for vertex in zvertices)
        n_allowed = tuple(defects[vertex] for vertex in nvertices)

        for z_supports in propagation.labeled_support_families(z_allowed):
            counters["labeled_support_families"] += 1
            propagated = propagation.analyze_support_assignment(
                graph_n,
                z_supports,
                n_allowed,
                zero_forcing,
                clique_solver,
            )
            if propagated.failure is not None:
                propagation_failures[propagated.failure] += 1
                continue
            counters["propagation_passing_families"] += 1

            link_failure = link_caps.first_link_cap_failure(
                z_supports, propagated.propagated_masks
            )
            if link_failure is not None:
                counters["link_cap_failing_families"] += 1
                link_failure_sizes[str(link_failure.subset_size)] += 1
                continue
            counters["link_cap_passing_families"] += 1

            sparse = sparse_value.check_small_support_masks(
                graph_n, propagated.propagated_masks
            )
            counters["sparse_calls"] += 1
            counters["sparse_assignments_checked"] += sparse.assignments_checked
            sparse_failures.update(dict(sparse.failures))
            if not sparse.feasible:
                counters["sparse_failing_families"] += 1
                continue
            counters["sparse_passing_families"] += 1
            survivor_witness = {
                "zmask": zmask,
                "z_vertices": [outside[vertex] for vertex in zvertices],
                "z_supports": list(z_supports),
                "propagated_masks": list(propagated.propagated_masks),
                "small_support_witness": (
                    None if sparse.first_witness is None else list(sparse.first_witness)
                ),
            }
            break
        if survivor_witness is not None:
            break

    if survivor_witness is not None:
        status = STATUS_SURVIVOR
    elif not counters["baseline_passing_covers"]:
        status = STATUS_GRAPH
    elif not counters["propagation_passing_families"]:
        status = STATUS_SUPPORT
    elif not counters["link_cap_passing_families"]:
        status = STATUS_LINK
    else:
        status = STATUS_SPARSE

    return {
        "seed_mask": seed_mask,
        "apex": apex,
        "added_edges": added_edges,
        "status": status,
        "counters": dict(counters),
        "baseline_failure_counts": dict(baseline_failures),
        "propagation_failure_counts": dict(propagation_failures),
        "link_failure_subset_sizes": dict(link_failure_sizes),
        "sparse_failure_counts": dict(sparse_failures),
        "survivor_witness": survivor_witness,
    }


def compact_augmentation(result: dict) -> dict:
    counters = result["counters"]
    return {
        "seed_mask": result["seed_mask"],
        "apex": result["apex"],
        "added_edges": result["added_edges"],
        "status": result["status"],
        "eligible_covers": counters.get("eligible_covers", 0),
        "small_covers": counters.get("small_covers", 0),
        "baseline_passing_covers": counters.get("baseline_passing_covers", 0),
        "labeled_support_families": counters.get("labeled_support_families", 0),
        "propagation_passing_families": counters.get(
            "propagation_passing_families", 0
        ),
        "link_cap_passing_families": counters.get("link_cap_passing_families", 0),
        "sparse_calls": counters.get("sparse_calls", 0),
        "sparse_assignments_checked": counters.get(
            "sparse_assignments_checked", 0
        ),
    }


def aggregate_detail(target: Counter[str], values: dict[str, int]) -> None:
    target.update({str(name): int(count) for name, count in values.items()})


def evaluate_record(record: dict) -> dict:
    started = time.perf_counter()
    adjacency = tuple(map(int, record["adjacency"]))
    validate_k6_graph(adjacency)
    k7.validate_graph(adjacency)
    if next(k7.clique_masks(adjacency, 7), None) is not None:
        raise ValueError("virtual-K7 pilot input unexpectedly already contains K7")
    seeds = [tuple(vertices(mask)) for mask in k6_clique_masks(adjacency, 6)]
    if len(seeds) != int(record["k6_count"]):
        raise AssertionError("stored K6 count changed")

    support_solver = k7.SupportSolver()
    zero_forcing = k7.ZeroForcingSolver()
    clique_solver = k7.CliqueStructureSolver()
    status_counts: Counter[str] = Counter()
    effort: Counter[str] = Counter()
    baseline_failures: Counter[str] = Counter()
    propagation_failures: Counter[str] = Counter()
    link_failure_sizes: Counter[str] = Counter()
    sparse_failures: Counter[str] = Counter()
    representatives: dict[str, dict] = {}
    augmentations = []
    seed_summaries = []

    full = (1 << len(adjacency)) - 1
    for seed in seeds:
        seed_mask = sum(1 << vertex for vertex in seed)
        viable = []
        for apex in vertices(full & ~seed_mask):
            result = analyze_virtual_seed(
                adjacency,
                seed,
                apex,
                support_solver,
                zero_forcing,
                clique_solver,
            )
            status = result["status"]
            status_counts[status] += 1
            effort.update(result["counters"])
            aggregate_detail(baseline_failures, result["baseline_failure_counts"])
            aggregate_detail(
                propagation_failures, result["propagation_failure_counts"]
            )
            aggregate_detail(link_failure_sizes, result["link_failure_subset_sizes"])
            aggregate_detail(sparse_failures, result["sparse_failure_counts"])
            if status == STATUS_SURVIVOR:
                viable.append(apex)
            if status not in representatives:
                representatives[status] = {
                    "seed": list(seed),
                    "apex": apex,
                    "detail": result,
                }
            augmentations.append(compact_augmentation(result))

        possible_pairs = [
            [first, second]
            for first, second in combinations(viable, 2)
            if not (adjacency[first] & (1 << second))
        ]
        seed_summaries.append(
            {
                "seed_mask": seed_mask,
                "viable_apices": viable,
                "viable_apex_count": len(viable),
                "possible_opposite_empty_pairs": possible_pairs,
                "possible_opposite_empty_pair_count": len(possible_pairs),
            }
        )

    expected_augmentations = 13 * len(seeds)
    if len(augmentations) != expected_augmentations:
        raise AssertionError("not every K6-seed/outside-apex augmentation was tested")
    if sum(status_counts.values()) != expected_augmentations:
        raise AssertionError("augmentation status partition changed")
    return {
        "index": int(record["index"]),
        "stratum": int(record["stratum"]),
        "k6_seeds": len(seeds),
        "augmentations_checked": expected_augmentations,
        "augmentation_status_counts": dict(status_counts),
        "seeds_forcing_all_defects_nonempty": sum(
            not summary["viable_apices"] for summary in seed_summaries
        ),
        "seeds_allowing_at_most_one_empty": sum(
            not summary["possible_opposite_empty_pairs"] for summary in seed_summaries
        ),
        "seed_summaries": seed_summaries,
        "augmentations": augmentations,
        "effort": dict(effort),
        "baseline_failure_counts": dict(baseline_failures),
        "propagation_failure_counts": dict(propagation_failures),
        "link_failure_subset_sizes": dict(link_failure_sizes),
        "sparse_failure_counts": dict(sparse_failures),
        "representatives": representatives,
        "elapsed_seconds": time.perf_counter() - started,
    }


def aggregate_graph_results(results: Sequence[dict]) -> dict:
    totals: Counter[str] = Counter()
    statuses: Counter[str] = Counter()
    effort: Counter[str] = Counter()
    baseline: Counter[str] = Counter()
    propagated: Counter[str] = Counter()
    links: Counter[str] = Counter()
    sparse: Counter[str] = Counter()
    representatives: dict[str, dict] = {}
    for graph in results:
        totals.update(
            graphs=1,
            k6_seeds=graph["k6_seeds"],
            augmentations=graph["augmentations_checked"],
            seeds_forcing_all_defects_nonempty=graph[
                "seeds_forcing_all_defects_nonempty"
            ],
            seeds_allowing_at_most_one_empty=graph[
                "seeds_allowing_at_most_one_empty"
            ],
        )
        statuses.update(graph["augmentation_status_counts"])
        effort.update(graph["effort"])
        baseline.update(graph["baseline_failure_counts"])
        propagated.update(graph["propagation_failure_counts"])
        links.update(graph["link_failure_subset_sizes"])
        sparse.update(graph["sparse_failure_counts"])
        for status, witness in graph["representatives"].items():
            representatives.setdefault(status, {"index": graph["index"], **witness})
    eliminated = sum(statuses[status] for status in STATUSES[:-1])
    totals["augmentation_branches_eliminated"] = eliminated
    totals["augmentation_branches_surviving"] = statuses[STATUS_SURVIVOR]
    return {
        "totals": dict(totals),
        "augmentation_status_counts": dict(statuses),
        "effort": dict(effort),
        "baseline_failure_counts": dict(baseline),
        "propagation_failure_counts": dict(propagated),
        "link_failure_subset_sizes": dict(links),
        "sparse_failure_counts": dict(sparse),
        "representatives": representatives,
    }


def run(output: Path, workers: int) -> dict:
    records, residue_indices, parent_report = load_residue()
    selected = stratified_selection(records)
    selection = [int(record["index"]) for record in selected]
    started = time.perf_counter()
    if workers == 1:
        results = [evaluate_record(record) for record in selected]
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=worker_init) as pool:
            results = list(pool.map(evaluate_record, selected, chunksize=1))
    aggregate = aggregate_graph_results(results)

    augmentation_keys = [
        [graph["index"], row["seed_mask"], row["apex"]]
        for graph in results
        for row in graph["augmentations"]
    ]
    eliminated_keys = [
        [graph["index"], row["seed_mask"], row["apex"]]
        for graph in results
        for row in graph["augmentations"]
        if row["status"] != STATUS_SURVIVOR
    ]
    report = {
        "schema": "d6-k6-virtual-k7-pilot-v1",
        "status": "PILOT_COMPLETE",
        "scope": {
            "input": "ordered 822-graph exact K6 residue",
            "selection": (
                "four K6-count quartiles; five minimum SHA256(decimal index) "
                "priorities per quartile"
            ),
            "virtual_seed_only": True,
            "filters": [
                "generic K7 rank/support/Schur masks",
                "labeled zero-factor support propagation",
                "propagated clique-link caps",
                "one/two-defect sparse-value CSP",
            ],
            "archive_conditioned_filters_used": False,
        },
        "input_graphs": len(residue_indices),
        "input_indices_sha256": stable_hash(residue_indices),
        "selection": selection,
        "selection_sha256": stable_hash(selection),
        "augmentation_keys_sha256": stable_hash(augmentation_keys),
        "eliminated_augmentation_keys_sha256": stable_hash(eliminated_keys),
        "aggregate": aggregate,
        "graph_results": results,
        "runtime": {
            "command": (
                f"{sys.executable} {Path(__file__).name} --workers {workers} "
                f"--output {output}"
            ),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "sources": {
            Path(__file__).name: sha256(Path(__file__)),
            PARENT_REPORT.name: sha256(PARENT_REPORT),
            "d6_k6_bipartite_rank_input.json": sha256(
                ROOT / "d6_k6_bipartite_rank_input.json"
            ),
            "d6_k7_rank_reference.py": sha256(ROOT / "d6_k7_rank_reference.py"),
            "d6_k7_support_propagation.py": sha256(
                ROOT / "d6_k7_support_propagation.py"
            ),
            "d6_k7_propagated_link_caps.py": sha256(
                ROOT / "d6_k7_propagated_link_caps.py"
            ),
            "d6_k7_small_support_value.py": sha256(
                ROOT / "d6_k7_small_support_value.py"
            ),
        },
        "parent_rejected_indices": list(parent_report["rejected_indices"]),
        "semantics": {
            "augmentation": "adds only branch-forced unit edges from x to Q",
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defect_coordinates": "upper bounds and may be zero",
            "eliminated_augmentation": "proves only u_Q(x) is nonzero",
            "original_graph_rejections_claimed": 0,
            "survivor": "necessary filters pass; not a realization",
            "unresolved_is_rejected": False,
        },
    }
    atomic_json(output, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_virtual_k7_pilot_report.json",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    report = run(args.output, args.workers)
    print(
        json.dumps(
            {
                "selection": len(report["selection"]),
                "aggregate": report["aggregate"]["totals"],
                "status_counts": report["aggregate"]["augmentation_status_counts"],
                "wall_seconds": report["runtime"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
