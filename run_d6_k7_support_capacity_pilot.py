#!/usr/bin/env python3
"""One-core exact support-capacity pilot on the frozen K7 union residue.

The pilot combines only already certified cover eliminations (strict-H,
degree-one Schur, and degree-four tetrads) with the frozen labeled-support
propagation and sparse-value layers.  For every remaining labeled family it
then applies :mod:`d6_k7_support_capacity`.

The default is a deterministic 32-graph stratified sample of the exact
258-graph tetrad/pattern-954 complement.  ``--all-residue`` is available only
after the bounded profile justifies it.  Execution is deliberately serial.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_special_h_reference as strict_h
import d6_k7_support_capacity as capacity
import d6_k7_support_propagation as propagation
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
RANK_INPUT = ROOT / ".runs/d6_k7_rank_survivors.json"
UNION = ROOT / "d6_k7_rankone_pattern_union.json"
UNION_VERIFICATION = ROOT / "d6_k7_rankone_pattern_union_verification.json"
TETRAD_REPORT = ROOT / "d6_k7_rankone_tetrad_full_report.json"
TETRAD_VERIFICATION = ROOT / "d6_k7_rankone_tetrad_full_verification_report.json"
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
SPARSE_REPORT = ROOT / "d6_k7_small_support_value_full_report.json"
SPARSE_DECISIONS = ROOT / "d6_k7_small_support_value_full_decisions.tsv.gz"

EXPECTED_HASHES = {
    ".runs/d6_k7_rank_survivors.json": (
        "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
    ),
    "d6_k7_rankone_pattern_union.json": (
        "1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599"
    ),
    "d6_k7_rankone_pattern_union_verification.json": (
        "7b9d5e35bf787e5d05aa486f914416598c31836264a479fbbc1af117f1dfcc14"
    ),
    "d6_k7_rankone_tetrad_full_report.json": (
        "ae2075ac83abcdfc0b7d42f63b9515c4e48b40cf9c976178aa366d151a78996a"
    ),
    "d6_k7_rankone_tetrad_full_verification_report.json": (
        "1812524c835fd6635b9815c0f0c25e9a49d99ff5d3b042e4b194ac5d7898dc97"
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
    "d6_k7_small_support_value_full_report.json": (
        "cea64f3dde804e0766c5b77e373aa50338d5bee6e5e54f44749f4e387cf529aa"
    ),
    "d6_k7_small_support_value_full_decisions.tsv.gz": (
        "9f70e0e267fe97bc2f6a2890cae47d8b1c882974d5f054f16b4673bfb45ee003"
    ),
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_support_propagation.py": (
        "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
    ),
    "d6_k7_small_support_value.py": (
        "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f"
    ),
    "d6_k7_special_h_reference.py": (
        "e649633d7262941c9c4608ea694ae0a5d6f29c146d47060c3dae71a036a9d190"
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


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        ).stdout.strip()

    try:
        status = git("status", "--porcelain=v1", "--untracked-files=all")
        return {
            "available": True,
            "commit": git("rev-parse", "HEAD"),
            "branch": git("branch", "--show-current"),
            "dirty": bool(status),
            "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
            "porcelain_lines": status.splitlines(),
        }
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}


def verify_inputs() -> dict[str, str]:
    observed = {}
    for name, expected in EXPECTED_HASHES.items():
        path = ROOT / name
        value = sha256(path)
        if value != expected:
            raise ValueError(f"hash mismatch for {name}: {value} != {expected}")
        observed[name] = value
    union_verification = json.loads(UNION_VERIFICATION.read_text())
    tetrad_verification = json.loads(TETRAD_VERIFICATION.read_text())
    if union_verification.get("status") != "PASS":
        raise ValueError("rank-one/pattern union is not independently PASS")
    if tetrad_verification.get("errors") or (
        tetrad_verification.get("summary", {}).get("status_counts")
        != {"PASS": 12839}
    ):
        raise ValueError("tetrad campaign is not independently complete")
    sparse_report = json.loads(SPARSE_REPORT.read_text())
    if (
        sparse_report.get("decisions", {}).get("graphs") != 16_228
        or sparse_report.get("decisions", {}).get("survivors") != 13_033
    ):
        raise ValueError("frozen sparse-value campaign is not complete")
    return observed


def profile_signature(profile: dict) -> tuple[bool, bool, bool]:
    return (
        bool(profile["has_saturating_cover"]),
        bool(profile["has_no_near_clique_cover"]),
        bool(profile["has_tetrad_resistant_near_clique_cover"]),
    )


def stratified_indices(profiles: list[dict], limit: int) -> list[int]:
    """Round-robin deterministic hash order across cover-mechanism strata."""

    groups: dict[tuple[bool, bool, bool], list[int]] = defaultdict(list)
    for profile in profiles:
        groups[profile_signature(profile)].append(int(profile["index"]))
    for signature, indices in groups.items():
        indices.sort(key=lambda index: hashlib.sha256(
            f"d6-support-capacity-v1:{signature}:{index}".encode("ascii")
        ).digest())
    signatures = sorted(groups)
    selected = []
    offset = 0
    while len(selected) < min(limit, len(profiles)):
        progressed = False
        for signature in signatures:
            values = groups[signature]
            if offset < len(values):
                selected.append(values[offset])
                progressed = True
                if len(selected) == min(limit, len(profiles)):
                    break
        if not progressed:
            break
        offset += 1
    return selected


def load_jsonl_witness_keys(
    path: Path, selected: set[int], field: str
) -> dict[int, set[tuple[tuple[int, ...], int]]]:
    output: dict[int, set[tuple[tuple[int, ...], int]]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            entry = json.loads(line)
            index = int(entry["index"])
            if index not in selected:
                continue
            for witness in entry[field]:
                key = (tuple(int(v) for v in witness["seed"]), int(witness["zmask"]))
                if key in output[index]:
                    raise ValueError(f"duplicate {field} key at graph {index}")
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
        raise ValueError("selected graph missing from tetrad decisions")
    return output


def load_sparse_survivors(selected: set[int]) -> None:
    seen = set()
    with gzip.open(SPARSE_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            if index in selected:
                if row["decision"] != "SURVIVOR":
                    raise ValueError(f"selected graph {index} failed sparse value")
                seen.add(index)
    if seen != selected:
        raise ValueError("selected graph missing from sparse-value decisions")


def current_cover_status(
    adj: tuple[int, ...],
    outside: list[int],
    defects: tuple[int, ...],
    zmask: int,
    baseline: reference.CoverAnalysis,
    key: tuple[tuple[int, ...], int],
    prior_failures: set[tuple[tuple[int, ...], int]],
    tetrad_failures: set[tuple[tuple[int, ...], int]],
) -> str:
    if baseline.enhanced_joint_failed:
        return "baseline"
    nvertices = [
        outside[position] for position in range(len(outside))
        if not zmask & (1 << position)
    ]
    graph_n = reference.induced_graph(adj, nvertices)
    saturating = tuple(reference.clique_masks(
        graph_n, baseline.k_rank_upper
    )) if baseline.k_rank_upper else ()
    if any(strict_h.assess_saturating_clique(graph_n, clique).failed
           for clique in saturating):
        return "strict_h"
    if key in prior_failures:
        return "prior_dual"
    if key in tetrad_failures:
        return "tetrad"
    return "passing"


def add_search_counts(counts: Counter, result: capacity.CapacitySearchResult) -> None:
    counts["capacity_nodes"] += result.nodes
    counts["capacity_flow_checks"] += result.flow_checks
    counts["capacity_flow_prunes"] += result.flow_prunes
    counts["capacity_empty_domain_prunes"] += result.empty_domain_prunes
    counts["capacity_edge_prunes"] += result.edge_prunes
    counts["capacity_type_prunes"] += result.capacity_prunes
    counts["capacity_maximum_depth_observed"] = max(
        counts["capacity_maximum_depth_observed"], result.maximum_depth
    )


def analyze_graph(
    graph: dict,
    prior_failures: set[tuple[tuple[int, ...], int]],
    tetrad_failures: set[tuple[tuple[int, ...], int]],
    expected_tetrad_row: dict,
    node_limit: int,
) -> dict:
    started = time.perf_counter()
    index = int(graph["index"])
    adj = tuple(int(row) for row in graph["adjacency"])
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    counts: Counter = Counter()
    seed_records = []
    current_passing_cover_keys = set()

    for seed_mask in reference.clique_masks(adj, 7):
        counts["seeds"] += 1
        seed, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        seed_pre_pass = False
        seed_capacity_pass = False
        seed_unresolved = False
        seed_counts: Counter = Counter()
        cover_records = []
        for zmask in reference.eligible_covers(ladj, eligible):
            counts["covers"] += 1
            seed_counts["covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            key = (tuple(seed), zmask)
            status = current_cover_status(
                adj, outside, defects, zmask, baseline, key,
                prior_failures, tetrad_failures,
            )
            counts[f"cover_{status}"] += 1
            seed_counts[f"cover_{status}"] += 1
            if status != "passing":
                continue
            current_passing_cover_keys.add(key)
            zvertices = tuple(reference.bits(zmask))
            nlocal = tuple(
                vertex for vertex in range(len(outside))
                if not zmask & (1 << vertex)
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nlocal]
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nlocal)
            cover_pre_pass = False
            cover_capacity_pass = False
            cover_unresolved = False
            family_counts: Counter = Counter()
            first_capacity_witness = None
            first_pre_capacity_failure = None
            first_infeasible = None
            for z_supports in propagation.labeled_support_families(z_allowed):
                counts["labeled_z_families"] += 1
                family_counts["labeled_z_families"] += 1
                propagated = propagation.analyze_support_assignment(
                    graph_n, z_supports, n_allowed,
                    zero_forcing, clique_solver,
                )
                if propagated.failure is not None:
                    name = f"propagation_failure:{propagated.failure}"
                    counts[name] += 1
                    family_counts[name] += 1
                    if first_pre_capacity_failure is None:
                        first_pre_capacity_failure = {
                            "stage": "propagation",
                            "reason": propagated.failure,
                            "z_supports": list(z_supports),
                            "propagated_masks": list(
                                propagated.propagated_masks
                            ),
                        }
                    continue
                counts["propagation_passing_families"] += 1
                family_counts["propagation_passing_families"] += 1
                sparse = sparse_value.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                if not sparse.feasible:
                    counts["sparse_value_failing_families"] += 1
                    family_counts["sparse_value_failing_families"] += 1
                    if first_pre_capacity_failure is None:
                        first_pre_capacity_failure = {
                            "stage": "sparse_value",
                            "reason_counts": dict(sparse.failures),
                            "z_supports": list(z_supports),
                            "propagated_masks": list(
                                propagated.propagated_masks
                            ),
                        }
                    continue
                cover_pre_pass = seed_pre_pass = True
                counts["pre_capacity_passing_families"] += 1
                family_counts["pre_capacity_passing_families"] += 1
                result = capacity.solve_support_capacity(
                    graph_n,
                    propagated.propagated_masks,
                    z_supports,
                    node_limit=node_limit,
                )
                add_search_counts(counts, result)
                add_search_counts(family_counts, result)
                counts[f"capacity_{result.status.lower()}_families"] += 1
                family_counts[f"capacity_{result.status.lower()}_families"] += 1
                if result.status == capacity.STATUS_FEASIBLE:
                    cover_capacity_pass = seed_capacity_pass = True
                    if first_capacity_witness is None:
                        first_capacity_witness = {
                            "z_supports": list(z_supports),
                            "propagated_masks": list(propagated.propagated_masks),
                            "n_actual_supports": list(result.witness or ()),
                        }
                elif result.status == capacity.STATUS_UNRESOLVED:
                    cover_unresolved = seed_unresolved = True
                elif first_infeasible is None:
                    first_infeasible = {
                        "z_supports": list(z_supports),
                        "propagated_masks": list(propagated.propagated_masks),
                        "nodes": result.nodes,
                        "flow_checks": result.flow_checks,
                        "reason": result.reason,
                    }
            if cover_pre_pass:
                counts["pre_capacity_passing_covers"] += 1
                seed_counts["pre_capacity_passing_covers"] += 1
            else:
                counts["joint_pre_capacity_failing_covers"] += 1
                seed_counts["joint_pre_capacity_failing_covers"] += 1
            if cover_capacity_pass:
                counts["capacity_passing_covers"] += 1
                seed_counts["capacity_passing_covers"] += 1
                cover_decision = "FEASIBLE"
            elif not cover_pre_pass:
                cover_decision = "JOINT_PRE_CAPACITY_INFEASIBLE"
            elif cover_unresolved:
                counts["capacity_unresolved_covers"] += 1
                seed_counts["capacity_unresolved_covers"] += 1
                cover_decision = "UNRESOLVED"
            else:
                counts["capacity_incremental_failing_covers"] += 1
                seed_counts["capacity_incremental_failing_covers"] += 1
                cover_decision = "INFEASIBLE"
            cover_records.append({
                "seed": seed,
                "zmask": zmask,
                "decision": cover_decision,
                "pre_capacity_pass": cover_pre_pass,
                "counts": dict(sorted(family_counts.items())),
                "first_pre_capacity_failure": first_pre_capacity_failure,
                "first_capacity_witness": first_capacity_witness,
                "first_infeasible_family": first_infeasible,
            })
        seed_records.append({
            "seed": seed,
            "seed_mask": seed_mask,
            "pre_capacity_pass": seed_pre_pass,
            "capacity_pass": seed_capacity_pass,
            "unresolved": seed_unresolved,
            "counts": dict(sorted(seed_counts.items())),
            "current_cover_records": cover_records,
        })

    expected_passing = int(expected_tetrad_row["tetrad_passing_covers"])
    if len(current_passing_cover_keys) != expected_passing:
        raise ValueError(
            f"graph {index} reconstructed {len(current_passing_cover_keys)} "
            f"tetrad-passing covers, expected {expected_passing}"
        )
    pre_rejected_seeds = [
        entry["seed_mask"] for entry in seed_records
        if not entry["pre_capacity_pass"]
    ]
    capacity_rejected_seeds = [
        entry["seed_mask"] for entry in seed_records
        if entry["pre_capacity_pass"]
        and not entry["capacity_pass"]
        and not entry["unresolved"]
    ]
    unresolved_seeds = [
        entry["seed_mask"] for entry in seed_records
        if entry["pre_capacity_pass"]
        and not entry["capacity_pass"]
        and entry["unresolved"]
    ]
    if pre_rejected_seeds:
        decision = "JOINT_PRE_CAPACITY_REJECTED"
    elif capacity_rejected_seeds:
        decision = "CAPACITY_REJECTED"
    elif unresolved_seeds:
        decision = "UNRESOLVED"
    else:
        decision = "SURVIVOR"
    return {
        "index": index,
        "decision": decision,
        "pre_capacity_rejected_seeds": pre_rejected_seeds,
        "capacity_rejected_seeds": capacity_rejected_seeds,
        "unresolved_seeds": unresolved_seeds,
        "counts": dict(sorted(counts.items())),
        "seed_records": seed_records,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=32)
    parser.add_argument("--all-residue", action="store_true")
    parser.add_argument("--node-limit", type=int, default=200_000)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_support_capacity_pilot_report.json",
    )
    args = parser.parse_args()
    if args.limit <= 0:
        raise SystemExit("--limit must be positive")
    if args.node_limit < 0:
        raise SystemExit("--node-limit must be nonnegative")

    dependency_hashes = verify_inputs()
    union = json.loads(UNION.read_text(encoding="utf-8"))
    profiles = list(union["cover_structure"]["residue_profiles"])
    residue = [int(value) for value in union["sets"]["exact_residue"]["indices"]]
    if len(profiles) != 258 or len(residue) != 258:
        raise ValueError("unexpected exact union residue population")
    selected_indices = (
        residue if args.all_residue else stratified_indices(profiles, args.limit)
    )
    selected = set(selected_indices)
    profile_by_index = {int(item["index"]): item for item in profiles}
    selected_profiles = [profile_by_index[index] for index in selected_indices]
    signature_counts = Counter(
        str(profile_signature(profile)) for profile in selected_profiles
    )

    rank_payload = json.loads(RANK_INPUT.read_text(encoding="utf-8"))
    graph_by_index = {
        int(graph["index"]): graph for graph in rank_payload["graphs"]
    }
    if not selected <= set(graph_by_index):
        raise ValueError("selected union residue graph absent from rank input")
    prior_failures = load_jsonl_witness_keys(
        PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
    )
    tetrad_failures = load_jsonl_witness_keys(
        TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
    )
    tetrad_rows = load_tetrad_rows(selected)
    load_sparse_survivors(selected)

    # The standard realizable 18-point lower bound has no K7 seed, so this
    # K7-conditioned layer is correctly reported as not applicable.
    positive = tuple(lower_bound_18_graph())
    reference.validate_graph(positive)
    positive_k7 = tuple(reference.clique_masks(positive, 7))
    if positive_k7:
        raise AssertionError("18-point lower-bound control unexpectedly has K7")

    started_unix = time.time()
    started = time.monotonic()
    per_graph = []
    totals: Counter = Counter()
    for position, index in enumerate(selected_indices, 1):
        result = analyze_graph(
            graph_by_index[index],
            prior_failures.get(index, set()),
            tetrad_failures.get(index, set()),
            tetrad_rows[index],
            args.node_limit,
        )
        per_graph.append(result)
        graph_counts = dict(result["counts"])
        graph_maximum_depth = graph_counts.pop(
            "capacity_maximum_depth_observed", 0
        )
        totals.update(graph_counts)
        totals["capacity_maximum_depth_observed"] = max(
            totals["capacity_maximum_depth_observed"], graph_maximum_depth
        )
        totals[f"graph_{result['decision'].lower()}"] += 1
        print(
            f"{position}/{len(selected_indices)} graph {index}: "
            f"{result['decision']} ({result['elapsed_seconds']:.3f}s)",
            flush=True,
        )
    elapsed = time.monotonic() - started
    report = {
        "schema": 1,
        "kind": "d6_k7_support_capacity_pilot",
        "status": "COMPLETE",
        "claim_scope": (
            "Exact serial bounded pilot. INFEASIBLE uses only nonempty actual "
            "support containment, required-edge intersection, and exact-type "
            "capacity |S| after frozen exact filters. SURVIVOR is a nonclaim."
        ),
        "selection": {
            "population": "exact 258-graph tetrad/pattern-954 union residue",
            "population_count": len(residue),
            "population_indices_sha256": stable_hash(residue),
            "mode": "all_residue" if args.all_residue else "stratified_sample",
            "selected_count": len(selected_indices),
            "selected_indices": selected_indices,
            "selected_indices_sha256": stable_hash(selected_indices),
            "signature_counts": dict(sorted(signature_counts.items())),
        },
        "node_limit_per_family": args.node_limit,
        "dependency_sha256": dependency_hashes,
        "source_sha256": sha256(Path(__file__)),
        "capacity_source_sha256": sha256(ROOT / "d6_k7_support_capacity.py"),
        "capacity_proof_sha256": (
            sha256(ROOT / "d6_k7_exact_support_multiplicity.md")
        ),
        "test_source_sha256": sha256(
            ROOT / "test_d6_k7_support_capacity.py"
        ),
        "git": git_provenance(),
        "command": [sys.executable, *sys.argv],
        "environment": {
            "workers": 1,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
        "started_unix": started_unix,
        "elapsed_seconds": elapsed,
        "positive_18_control": {
            "vertices": len(positive),
            "K7_seeds": 0,
            "status": "NOT_APPLICABLE_NO_K7",
        },
        "totals": dict(sorted(totals.items())),
        "per_graph": per_graph,
    }
    atomic_json(args.output, report)
    print(
        f"wrote {args.output}; SHA-256 {sha256(args.output)}; "
        f"elapsed {elapsed:.3f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
