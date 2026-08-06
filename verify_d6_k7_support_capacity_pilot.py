#!/usr/bin/env python3
"""Independent replay of the bounded K7 support-capacity pilot.

This checker does not import the new capacity kernel or pilot runner.  It
reconstructs the deterministic sample, independently replays every reported
joint pre-capacity rejection using the previously independent support and
sparse-value implementations, and directly checks every recorded capacity
witness from the defining finite constraints.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import d6_k7_special_h_reference as strict_h
import verify_d6_k7_sparse_value_full as sparse_independent
import verify_d6_k7_support_full as support_independent


ROOT = Path(__file__).resolve().parent
RANK_INPUT = ROOT / ".runs/d6_k7_rank_survivors.json"
UNION = ROOT / "d6_k7_rankone_pattern_union.json"
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"

EXPECTED_HASHES = {
    ".runs/d6_k7_rank_survivors.json": (
        "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
    ),
    "d6_k7_rankone_pattern_union.json": (
        "1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599"
    ),
    "d6_k7_rankone_tetrad_full_decisions.tsv.gz": (
        "2552ce91f52727d9beda60d99d534c1d5a0be3f93686ab321f223d8e0d11dfc1"
    ),
    "d6_k7_rankone_tetrad_full_certificates.jsonl.gz": (
        "f719dbb6492fc20fb3103cb79079567aa2cad163835f97e750412c3d27897c50"
    ),
    "d6_k7_positive_dual_full_certificates.jsonl.gz": (
        "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
    ),
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_special_h_reference.py": (
        "e649633d7262941c9c4608ea694ae0a5d6f29c146d47060c3dae71a036a9d190"
    ),
    "verify_d6_k7_support_full.py": (
        "c6b9dce47618bd9a4247eb5a9f093682265447729e1868c4dc0677a24f51c628"
    ),
    "verify_d6_k7_sparse_value_full.py": (
        "6c0012e66d0b89e280f3cebed7a61f24cc1c6e128e726663a7b246e7da0edbbc"
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
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def signature(profile: dict) -> tuple[bool, bool, bool]:
    return (
        bool(profile["has_saturating_cover"]),
        bool(profile["has_no_near_clique_cover"]),
        bool(profile["has_tetrad_resistant_near_clique_cover"]),
    )


def independent_sample(profiles: Sequence[dict], limit: int) -> list[int]:
    groups: dict[tuple[bool, bool, bool], list[int]] = defaultdict(list)
    for profile in profiles:
        groups[signature(profile)].append(int(profile["index"]))
    for key, values in groups.items():
        values.sort(key=lambda index: hashlib.sha256(
            f"d6-support-capacity-v1:{key}:{index}".encode("ascii")
        ).digest())
    keys = sorted(groups)
    output = []
    offset = 0
    while len(output) < min(limit, len(profiles)):
        changed = False
        for key in keys:
            if offset < len(groups[key]):
                output.append(groups[key][offset])
                changed = True
                if len(output) == min(limit, len(profiles)):
                    break
        if not changed:
            break
        offset += 1
    return output


def load_witness_keys(
    path: Path, selected: set[int], field: str
) -> dict[int, set[tuple[tuple[int, ...], int]]]:
    output: dict[int, set[tuple[tuple[int, ...], int]]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            item = json.loads(line)
            index = int(item["index"])
            if index not in selected:
                continue
            for witness in item[field]:
                key = (tuple(int(v) for v in witness["seed"]), int(witness["zmask"]))
                if key in output[index]:
                    raise ValueError("duplicate inherited cover witness")
                output[index].add(key)
    return output


def load_tetrad_rows(selected: set[int]) -> dict[int, dict]:
    output = {}
    with gzip.open(TETRAD_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            if index in selected:
                output[index] = row
    if set(output) != selected:
        raise ValueError("selected graph absent from tetrad decisions")
    return output


def direct_capacity_witness_check(
    graph_n: Sequence[int],
    propagated_masks: Sequence[int],
    fixed_supports: Sequence[int],
    actual_supports: Sequence[int],
) -> None:
    """Check the finite constraints without the production capacity module."""

    graph = tuple(int(row) for row in graph_n)
    masks = tuple(int(mask) for mask in propagated_masks)
    fixed = tuple(int(mask) for mask in fixed_supports)
    actual = tuple(int(mask) for mask in actual_supports)
    if len(graph) != len(masks) or len(actual) != len(masks):
        raise ValueError("capacity witness order mismatch")
    full = (1 << len(graph)) - 1
    for vertex, row in enumerate(graph):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError("bad induced graph")
        for neighbour in support_independent.bit_positions(row):
            if not graph[neighbour] & (1 << vertex):
                raise ValueError("asymmetric induced graph")
    if any(mask <= 0 or mask >= 128 for mask in (*masks, *fixed, *actual)):
        raise ValueError("bad seven-bit support")
    if any(support & ~mask for support, mask in zip(actual, masks)):
        raise ValueError("actual support escapes propagated mask")
    for first, second in combinations(range(len(graph)), 2):
        if graph[first] & (1 << second) and not actual[first] & actual[second]:
            raise ValueError("required edge has disjoint actual supports")
    counts = Counter((*fixed, *actual))
    for support, number in counts.items():
        if number > support.bit_count():
            raise ValueError("exact type multiplicity exceeds |S|")


def independently_current_cover_passes(
    adj: tuple[int, ...],
    outside: tuple[int, ...],
    defects: tuple[int, ...],
    seed: tuple[int, ...],
    zmask: int,
    baseline: reference.CoverAnalysis,
    prior_failures: set[tuple[tuple[int, ...], int]],
    tetrad_failures: set[tuple[tuple[int, ...], int]],
) -> bool:
    if baseline.enhanced_joint_failed:
        return False
    nvertices = tuple(
        outside[position] for position in range(len(outside))
        if not zmask & (1 << position)
    )
    graph_n = support_independent.independent_induced_graph(adj, nvertices)
    saturating = tuple(reference.clique_masks(
        graph_n, baseline.k_rank_upper
    )) if baseline.k_rank_upper else ()
    if any(strict_h.assess_saturating_clique(graph_n, clique).failed
           for clique in saturating):
        return False
    key = (seed, zmask)
    return key not in prior_failures and key not in tetrad_failures


def replay_failing_seed(
    graph: dict,
    seed_mask: int,
    prior_failures: set[tuple[tuple[int, ...], int]],
    tetrad_failures: set[tuple[tuple[int, ...], int]],
) -> dict:
    """Exhaust every current-passing cover/family with independent kernels."""

    index = int(graph["index"])
    adj = tuple(int(row) for row in graph["adjacency"])
    reference.validate_graph(adj)
    outside, defects, ladj, eligible = support_independent.independent_seed_instance(
        adj, seed_mask
    )
    seed = tuple(support_independent.bit_positions(seed_mask))
    covers = tuple(support_independent.independent_eligible_covers(ladj, eligible))
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    total_term_rank = reference.matching_size(defects)
    counts: Counter = Counter()
    for zmask in covers:
        baseline = reference.analyze_cover(
            adj, outside, defects, zmask,
            support_solver, zero_forcing, total_term_rank, clique_solver,
        )
        if not independently_current_cover_passes(
            adj, outside, defects, seed, zmask, baseline,
            prior_failures, tetrad_failures,
        ):
            continue
        counts["current_passing_covers"] += 1
        zvertices = tuple(support_independent.bit_positions(zmask))
        nlocal = tuple(
            vertex for vertex in range(len(outside))
            if not zmask & (1 << vertex)
        )
        z_allowed = tuple(defects[vertex] for vertex in zvertices)
        n_allowed = tuple(defects[vertex] for vertex in nlocal)
        graph_n = support_independent.independent_induced_graph(
            adj, tuple(outside[vertex] for vertex in nlocal)
        )
        for fixed in support_independent.independent_labeled_support_families(
            z_allowed
        ):
            counts["labeled_z_families"] += 1
            masks, _ = sparse_independent.propagated_masks(n_allowed, fixed)
            failure = sparse_independent.simple_propagation_failure(graph_n, masks)
            if failure is not None:
                counts[f"propagation_failure:{failure}"] += 1
                continue
            counts["propagation_passing_families"] += 1
            sparse = sparse_independent.independent_small_support_check(
                graph_n, masks
            )
            if sparse.feasible:
                raise AssertionError(
                    f"graph {index}, seed {seed_mask} independently survives "
                    f"at cover {zmask}, supports {fixed}"
                )
            counts["sparse_value_failing_families"] += 1
    if not counts["current_passing_covers"]:
        raise AssertionError("reported joint rejection has no current-passing cover")
    return {
        "index": index,
        "seed_mask": seed_mask,
        "seed": list(seed),
        "counts": dict(sorted(counts.items())),
    }


def verify_recorded_cover_witness(
    graph: dict, seed_record: dict, cover_record: dict
) -> None:
    adj = tuple(int(row) for row in graph["adjacency"])
    seed_mask = int(seed_record["seed_mask"])
    outside, defects, _, _ = support_independent.independent_seed_instance(
        adj, seed_mask
    )
    zmask = int(cover_record["zmask"])
    zvertices = tuple(support_independent.bit_positions(zmask))
    nlocal = tuple(
        vertex for vertex in range(len(outside))
        if not zmask & (1 << vertex)
    )
    z_allowed = tuple(defects[vertex] for vertex in zvertices)
    n_allowed = tuple(defects[vertex] for vertex in nlocal)
    graph_n = support_independent.independent_induced_graph(
        adj, tuple(outside[vertex] for vertex in nlocal)
    )
    witness = cover_record.get("first_capacity_witness")
    if not isinstance(witness, dict):
        raise ValueError("pre-capacity-passing cover lacks capacity witness")
    fixed = tuple(int(mask) for mask in witness["z_supports"])
    if len(fixed) != len(z_allowed) or any(
        support & ~allowed for support, allowed in zip(fixed, z_allowed)
    ) or not support_independent.independently_valid_support_family(fixed):
        raise ValueError("recorded fixed support family is invalid")
    masks, _ = sparse_independent.propagated_masks(n_allowed, fixed)
    if tuple(witness["propagated_masks"]) != masks:
        raise ValueError("recorded propagated masks do not reconstruct")
    if sparse_independent.simple_propagation_failure(graph_n, masks) is not None:
        raise ValueError("capacity witness family fails propagation")
    sparse = sparse_independent.independent_small_support_check(graph_n, masks)
    if not sparse.feasible:
        raise ValueError("capacity witness family fails sparse-value replay")
    direct_capacity_witness_check(
        graph_n, masks, fixed, witness["n_actual_supports"]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_support_capacity_pilot_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_support_capacity_pilot_verification.json",
    )
    args = parser.parse_args()
    if sha256(args.report) != args.report_sha256:
        raise SystemExit("pilot report hash differs from explicit pin")
    observed = {}
    for name, expected in EXPECTED_HASHES.items():
        value = sha256(ROOT / name)
        if value != expected:
            raise ValueError(f"hash mismatch for {name}: {value} != {expected}")
        observed[name] = value

    report = json.loads(args.report.read_text(encoding="utf-8"))
    if report.get("status") != "COMPLETE" or report.get("schema") != 1:
        raise ValueError("pilot report is not complete schema one")
    if report["source_sha256"] != sha256(
        ROOT / "run_d6_k7_support_capacity_pilot.py"
    ) or report["capacity_source_sha256"] != sha256(
        ROOT / "d6_k7_support_capacity.py"
    ):
        raise ValueError("pilot source bindings do not match")
    union = json.loads(UNION.read_text(encoding="utf-8"))
    profiles = list(union["cover_structure"]["residue_profiles"])
    residue = [int(v) for v in union["sets"]["exact_residue"]["indices"]]
    selection = report["selection"]
    selected_indices = [int(v) for v in selection["selected_indices"]]
    if selection["mode"] != "stratified_sample":
        raise ValueError("bounded verifier expects the stratified pilot")
    independently_selected = independent_sample(profiles, len(selected_indices))
    if selected_indices != independently_selected:
        raise ValueError("deterministic stratified selection mismatch")
    if selection["population_indices_sha256"] != stable_hash(residue):
        raise ValueError("residue population hash mismatch")
    if selection["selected_indices_sha256"] != stable_hash(selected_indices):
        raise ValueError("selected-index hash mismatch")
    selected = set(selected_indices)

    rank_payload = json.loads(RANK_INPUT.read_text(encoding="utf-8"))
    graph_by_index = {
        int(graph["index"]): graph for graph in rank_payload["graphs"]
    }
    prior_failures = load_witness_keys(
        PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
    )
    tetrad_failures = load_witness_keys(
        TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
    )
    tetrad_rows = load_tetrad_rows(selected)

    per_graph = report["per_graph"]
    if [int(item["index"]) for item in per_graph] != selected_indices:
        raise ValueError("per-graph report order differs from selection")
    replayed_rejections = []
    cover_witnesses = 0
    graph_decisions: Counter = Counter()
    aggregate: Counter = Counter()
    maximum_depth = 0
    for entry in per_graph:
        index = int(entry["index"])
        graph = graph_by_index[index]
        graph_decisions[entry["decision"]] += 1
        graph_counts = dict(entry["counts"])
        maximum_depth = max(
            maximum_depth,
            int(graph_counts.pop("capacity_maximum_depth_observed", 0)),
        )
        aggregate.update({key: int(value) for key, value in graph_counts.items()})
        aggregate[f"graph_{entry['decision'].lower()}"] += 1
        expected_passing = int(tetrad_rows[index]["tetrad_passing_covers"])
        reported_passing = sum(
            int(seed["counts"].get("cover_passing", 0))
            for seed in entry["seed_records"]
        )
        if reported_passing != expected_passing:
            raise ValueError("reported current cover population mismatch")
        for seed_record in entry["seed_records"]:
            seed_has_witness = False
            for cover_record in seed_record["current_cover_records"]:
                if cover_record["pre_capacity_pass"]:
                    verify_recorded_cover_witness(graph, seed_record, cover_record)
                    cover_witnesses += 1
                    seed_has_witness = True
            if seed_record["capacity_pass"] != seed_has_witness:
                raise ValueError("seed capacity decision/witness mismatch")
        if entry["decision"] == "JOINT_PRE_CAPACITY_REJECTED":
            seeds = [int(value) for value in entry["pre_capacity_rejected_seeds"]]
            if not seeds:
                raise ValueError("joint rejection has no failing seed")
            for seed_mask in seeds:
                replayed_rejections.append(replay_failing_seed(
                    graph, seed_mask,
                    prior_failures.get(index, set()),
                    tetrad_failures.get(index, set()),
                ))
        elif entry["decision"] != "SURVIVOR":
            raise ValueError("pilot contains an unexpected decision class")

    aggregate["capacity_maximum_depth_observed"] = maximum_depth
    if dict(sorted(aggregate.items())) != report["totals"]:
        raise ValueError("report aggregate totals do not reconstruct")
    if graph_decisions != Counter({
        "SURVIVOR": 25,
        "JOINT_PRE_CAPACITY_REJECTED": 7,
    }):
        raise ValueError(f"unexpected pilot decision counts: {graph_decisions}")
    if report["totals"].get("capacity_incremental_failing_covers", 0):
        raise ValueError("bounded pilot unexpectedly has a capacity-only hit")
    if report["totals"].get("capacity_unresolved_covers", 0):
        raise ValueError("bounded pilot contains unresolved capacity covers")

    output = {
        "schema": 1,
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": args.report_sha256},
        "dependency_sha256": observed,
        "verifier_source_sha256": sha256(Path(__file__)),
        "selected_graphs": len(selected_indices),
        "graph_decisions": dict(sorted(graph_decisions.items())),
        "joint_rejected_graphs_replayed": 7,
        "failing_seeds_replayed": len(replayed_rejections),
        "capacity_cover_witnesses_checked": cover_witnesses,
        "capacity_incremental_hits": 0,
        "unresolved": 0,
        "replayed_rejections": replayed_rejections,
    }
    atomic_json(args.output, output)
    print(
        f"PASS: {len(selected_indices)} graphs, "
        f"{len(replayed_rejections)} failing seeds, "
        f"{cover_witnesses} direct capacity witnesses; "
        f"SHA-256 {sha256(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
