#!/usr/bin/env python3
"""Independent verifier for the K7 joint-support/pentad conjunction.

This checker does not import the conjunction builder or its production support
kernel.  It reads the 258 adjacencies embedded in the committed v2 residue,
independently enumerates every K7 seed and eligible cover, replays all 738
inherited-passing covers with the independent support and sparse-value
implementations, and reconstructs the seed-local disjunction.  It also
re-expands every one of the 36 pentad certificates in fresh subprocesses.

No ignored ``.runs`` artifact is read.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import d6_k7_special_h_reference as strict_h
import verify_d6_k7_ranktwo_pentad as pentad_checker
import verify_d6_k7_sparse_value_full as sparse_independent
import verify_d6_k7_support_full as support_independent
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
CURRENT_RESIDUE = ROOT / "d6_current_residue_manifest_v2.json"
JOINT_CHECKPOINT = ROOT / "d6_k7_joint_support_full_checkpoint.json"
JOINT_DECISIONS = ROOT / "d6_k7_joint_support_full_decisions.tsv.gz"
JOINT_CERTIFICATES = ROOT / "d6_k7_joint_support_full_certificates.jsonl.gz"
JOINT_REPORT = ROOT / "d6_k7_joint_support_full_report.json"
JOINT_VERIFICATION = ROOT / "d6_k7_joint_support_full_verification.json"
PENTAD_REPORT = ROOT / "d6_k7_ranktwo_pentad_full_report.json"
PENTAD_VERIFICATION = ROOT / "d6_k7_ranktwo_pentad_full_verification.json"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"

EXPECTED_HASHES = {
    "d6_current_residue_manifest_v2.json": (
        "961dac1f9b44bb541e2c5bd1626027826eeea7bc628bf5ff0a85ab26afa949f4"
    ),
    "d6_k7_joint_support_full_checkpoint.json": (
        "91638dbcde9e18feb598521c801578a5e95a6a3c6fe48a1568c079d1c1d8ef3a"
    ),
    "d6_k7_joint_support_full_decisions.tsv.gz": (
        "b6040e7796e3d7df511a4e460ad71082e63c074ac2cf2713c19d8f3b772bcc90"
    ),
    "d6_k7_joint_support_full_certificates.jsonl.gz": (
        "774d3d2119d5ce5eb965a689853d85704c44e081b12ff0947f360b4571d00522"
    ),
    "d6_k7_joint_support_full_report.json": (
        "b3fd2752a371f8909b79b50b8e2b6135f93423d0adf7286aee26d4a8c0643326"
    ),
    "d6_k7_joint_support_full_verification.json": (
        "41132fb0cb6d7fcb9c02d4e5171de9c0b8afd9c7eab4da2a2776f2518a1facd1"
    ),
    "d6_k7_ranktwo_pentad_full_report.json": (
        "b85067767a926861f6438d395b167538a02124a5c6e7f11471df2fb87c4839e4"
    ),
    "d6_k7_ranktwo_pentad_full_verification.json": (
        "8fdac4a9fa158e57ff8e9e79b19801baaf0f486525920e699f50d9c8c1340568"
    ),
    "d6_k7_positive_dual_full_certificates.jsonl.gz": (
        "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
    ),
    "d6_k7_rankone_tetrad_full_certificates.jsonl.gz": (
        "f719dbb6492fc20fb3103cb79079567aa2cad163835f97e750412c3d27897c50"
    ),
    "d6_k7_rankone_tetrad_full_decisions.tsv.gz": (
        "2552ce91f52727d9beda60d99d534c1d5a0be3f93686ab321f223d8e0d11dfc1"
    ),
}

EXPECTED_SOURCE_HASHES = {
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
    "verify_d6_k7_ranktwo_pentad.py": (
        "14cf8d3e12dd436bc164c7d17a7aa1e03aa2aa933d66ffe0a1335c605b45eeb3"
    ),
    "verify_profile_d6.py": (
        "ab74371f26a4ab3a8dd716ea5804bf86aa947a9f37f9f53e33fa5c4b481c549a"
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
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def read_witness_keys(
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
                key = (
                    tuple(int(vertex) for vertex in witness["seed"]),
                    int(witness["zmask"]),
                )
                if key in output[index]:
                    raise ValueError("duplicate inherited witness key")
                output[index].add(key)
    return output


def read_tetrad_rows(selected: set[int]) -> dict[int, dict]:
    output = {}
    with gzip.open(TETRAD_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            if index in selected:
                output[index] = row
    if set(output) != selected:
        raise ValueError("tetrad decision archive omits a K7 residue graph")
    return output


def read_joint_decisions(indices: Sequence[int]) -> dict[int, dict]:
    with gzip.open(JOINT_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != len(indices):
        raise ValueError("joint decision population mismatch")
    output = {}
    for ordinal, (index, row) in enumerate(zip(indices, rows)):
        if int(row["ordinal"]) != ordinal or int(row["index"]) != index:
            raise ValueError("joint decision order mismatch")
        output[index] = row
    return output


def independently_all_eligible_covers(
    ladj: Sequence[int], eligible: int
) -> tuple[int, ...]:
    """Enumerate the theorem's full size-at-most-seven cover universe."""

    full = (1 << len(ladj)) - 1
    output = []
    for zmask in range(1 << len(ladj)):
        if zmask & ~eligible or zmask.bit_count() > 7:
            continue
        remaining = full & ~zmask
        if all(
            not (ladj[vertex] & remaining)
            for vertex in support_independent.bit_positions(remaining)
        ):
            output.append(zmask)
    return tuple(output)


def inherited_cover_status(
    adj: tuple[int, ...],
    outside: tuple[int, ...],
    defects: tuple[int, ...],
    seed: tuple[int, ...],
    zmask: int,
    baseline: reference.CoverAnalysis,
    prior_failures: set[tuple[tuple[int, ...], int]],
    tetrad_failures: set[tuple[tuple[int, ...], int]],
) -> str:
    if baseline.enhanced_joint_failed:
        return "baseline"
    nvertices = tuple(
        outside[position] for position in range(len(outside))
        if not zmask & (1 << position)
    )
    graph_n = support_independent.independent_induced_graph(adj, nvertices)
    saturating = tuple(reference.clique_masks(
        graph_n, baseline.k_rank_upper
    )) if baseline.k_rank_upper else ()
    if any(
        strict_h.assess_saturating_clique(graph_n, clique).failed
        for clique in saturating
    ):
        return "strict_h"
    key = (seed, int(zmask))
    if key in prior_failures:
        return "prior_dual"
    if key in tetrad_failures:
        return "tetrad"
    return "passing"


def independently_classify_support(
    graph_n: tuple[int, ...],
    z_allowed: tuple[int, ...],
    n_allowed: tuple[int, ...],
) -> tuple[str, dict[str, int]]:
    counts: Counter = Counter()
    pre_capacity_pass = False
    for fixed in support_independent.independent_labeled_support_families(z_allowed):
        counts["labeled_z_families"] += 1
        masks, _ = sparse_independent.propagated_masks(n_allowed, fixed)
        failure = sparse_independent.simple_propagation_failure(graph_n, masks)
        if failure is not None:
            counts[f"propagation_failure:{failure}"] += 1
            continue
        counts["propagation_passing_families"] += 1
        sparse = sparse_independent.independent_small_support_check(graph_n, masks)
        if not sparse.feasible:
            counts["sparse_value_failing_families"] += 1
            continue
        counts["pre_capacity_passing_families"] += 1
        pre_capacity_pass = True
    return (
        "PASSING" if pre_capacity_pass else "INFEASIBLE",
        dict(sorted(counts.items())),
    )


def verify_graph(payload: tuple[dict, list, list]) -> dict:
    """Independent complete seed/cover/support quantifier for one graph."""

    graph, serialized_prior, serialized_tetrad = payload
    prior_failures = {
        (tuple(int(vertex) for vertex in seed), int(zmask))
        for seed, zmask in serialized_prior
    }
    tetrad_failures = {
        (tuple(int(vertex) for vertex in seed), int(zmask))
        for seed, zmask in serialized_tetrad
    }
    index = int(graph["index"])
    adj = tuple(int(row) for row in graph["adjacency"])
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    graph_counts: Counter = Counter()
    seed_quantifiers = []
    for seed_mask in reference.clique_masks(adj, 7):
        graph_counts["seeds"] += 1
        outside, defects, ladj, eligible = (
            support_independent.independent_seed_instance(adj, seed_mask)
        )
        outside = tuple(int(vertex) for vertex in outside)
        defects = tuple(int(mask) for mask in defects)
        seed = tuple(support_independent.bit_positions(seed_mask))
        covers = independently_all_eligible_covers(ladj, eligible)
        total_term_rank = reference.matching_size(defects)
        status_counts: Counter = Counter()
        current = []
        joint_failures = 0
        pre_capacity_passes = 0
        for zmask in covers:
            graph_counts["covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask,
                support_solver, zero_forcing, total_term_rank, clique_solver,
            )
            status = inherited_cover_status(
                adj, outside, defects, seed, zmask, baseline,
                prior_failures, tetrad_failures,
            )
            status_counts[status] += 1
            if status != "passing":
                continue
            graph_counts["current_passing_covers"] += 1
            zvertices = tuple(support_independent.bit_positions(zmask))
            nlocal = tuple(
                position for position in range(len(outside))
                if not zmask & (1 << position)
            )
            graph_n = support_independent.independent_induced_graph(
                adj, tuple(outside[position] for position in nlocal)
            )
            support_status, family_counts = independently_classify_support(
                graph_n,
                tuple(defects[position] for position in zvertices),
                tuple(defects[position] for position in nlocal),
            )
            if support_status == "INFEASIBLE":
                joint_failures += 1
                graph_counts["joint_pre_capacity_failing_covers"] += 1
            else:
                pre_capacity_passes += 1
                graph_counts["pre_capacity_passing_covers"] += 1
            current.append({
                "seed": list(seed),
                "zmask": int(zmask),
                "joint_pre_capacity_status": support_status,
                "family_counts": family_counts,
            })
        seed_quantifiers.append({
            "seed": list(seed),
            "seed_mask": int(seed_mask),
            "eligible_covers": len(covers),
            "cover_status_counts": dict(sorted(status_counts.items())),
            "joint_pre_capacity_failing_covers": joint_failures,
            "pre_capacity_passing_covers": pre_capacity_passes,
            "current_passing_covers": current,
        })
    return {
        "index": index,
        "counts": dict(sorted(graph_counts.items())),
        "seed_quantifiers": seed_quantifiers,
    }


def key_tuple(item: dict) -> tuple[int, tuple[int, ...], int]:
    return (
        int(item["graph_index"]),
        tuple(int(vertex) for vertex in item["seed"]),
        int(item["zmask"]),
    )


def key_payload(key: tuple[int, tuple[int, ...], int]) -> dict:
    return {"graph_index": key[0], "seed": list(key[1]), "zmask": key[2]}


def pentad_worker(payload: tuple[int, dict, dict]) -> dict:
    ordinal, target, record = payload
    return pentad_checker.verify_full_cover_record(target, record, ordinal)


def pentad_subprocess(payload: tuple[int, dict, dict]) -> dict:
    environment = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = "1"
    result = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--pentad-worker"],
        cwd=ROOT,
        env=environment,
        input=json.dumps(payload, sort_keys=True, separators=(",", ":")),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"pentad subprocess failed ({result.returncode}): "
            f"{result.stderr[-1000:]}"
        )
    return json.loads(result.stdout)


def validate_roots() -> tuple[list[int], dict[int, dict], dict[int, dict], dict, dict]:
    for name, expected in {**EXPECTED_HASHES, **EXPECTED_SOURCE_HASHES}.items():
        observed = sha256(ROOT / name)
        if observed != expected:
            raise ValueError(f"hash mismatch for {name}: {observed} != {expected}")
    current = json.loads(CURRENT_RESIDUE.read_text(encoding="utf-8"))
    if current.get("schema") != "d6-current-exact-residue-v2" or (
        current.get("status") != "COMPLETE"
    ):
        raise ValueError("bad v2 residue gate")
    k7 = current["classes"]["K7"]
    graphs = list(k7["graphs"])
    indices = [int(index) for index in k7["residue_indices"]]
    if (
        len(graphs) != 258
        or [int(graph["index"]) for graph in graphs] != indices
        or stable_hash(indices) != k7["residue_indices_sha256"]
        or stable_hash(graphs) != k7["graphs_sha256"]
    ):
        raise ValueError("embedded K7 corpus does not bind")
    graph_by_index = {int(graph["index"]): graph for graph in graphs}

    joint = json.loads(JOINT_REPORT.read_text(encoding="utf-8"))
    joint_verification = json.loads(JOINT_VERIFICATION.read_text(encoding="utf-8"))
    checkpoint = json.loads(JOINT_CHECKPOINT.read_text(encoding="utf-8"))
    if (
        joint.get("status") != "COMPLETE"
        or joint_verification.get("status") != "PASS"
        or joint_verification["report"]["sha256"] != EXPECTED_HASHES[
            "d6_k7_joint_support_full_report.json"
        ]
        or checkpoint.get("status") != "COMPLETE"
        or checkpoint.get("completed") != 258
        or checkpoint.get("config") != joint.get("configuration")
        or checkpoint.get("config_sha256") != joint.get("configuration_sha256")
    ):
        raise ValueError("joint-support root gate failed")
    for field, name in (
        ("checkpoint_copy_sha256", "d6_k7_joint_support_full_checkpoint.json"),
        ("decisions_sha256", "d6_k7_joint_support_full_decisions.tsv.gz"),
        ("certificates_sha256", "d6_k7_joint_support_full_certificates.jsonl.gz"),
    ):
        if joint["artifacts"][field] != EXPECTED_HASHES[name]:
            raise ValueError("joint artifact binding mismatch")
    rows = read_joint_decisions(indices)
    joint_rejected = [
        index for index in indices
        if rows[index]["status"] == "JOINT_PRE_CAPACITY_REJECTED"
    ]
    if joint_rejected != joint["summary"]["joint_rejected_indices"]:
        raise ValueError("joint rejection list mismatch")
    with gzip.open(JOINT_CERTIFICATES, "rt", encoding="utf-8") as stream:
        certificates = [json.loads(line) for line in stream if line.strip()]
    if [int(item["index"]) for item in certificates] != joint_rejected:
        raise ValueError("joint certificate list mismatch")

    pentad = json.loads(PENTAD_REPORT.read_text(encoding="utf-8"))
    pentad_verification = json.loads(
        PENTAD_VERIFICATION.read_text(encoding="utf-8")
    )
    if (
        pentad.get("status") != "COMPLETE"
        or pentad_verification.get("status") != "PASS"
        or pentad_verification["report"]["sha256"] != EXPECTED_HASHES[
            "d6_k7_ranktwo_pentad_full_report.json"
        ]
        or pentad_verification["checked"]["certificates"] != 36
    ):
        raise ValueError("pentad root gate failed")
    return indices, graph_by_index, rows, joint, pentad


def build_pentad_target(
    graph: dict,
    record: dict,
    prior_failures: set[tuple[tuple[int, ...], int]],
    tetrad_failures: set[tuple[tuple[int, ...], int]],
) -> dict:
    adj = tuple(int(row) for row in graph["adjacency"])
    seed = tuple(int(vertex) for vertex in record["seed"])
    seed_mask = sum(1 << vertex for vertex in seed)
    outside, defects, ladj, eligible = support_independent.independent_seed_instance(
        adj, seed_mask
    )
    outside = tuple(int(vertex) for vertex in outside)
    defects = tuple(int(mask) for mask in defects)
    zmask = int(record["zmask"])
    covers = independently_all_eligible_covers(ladj, eligible)
    if zmask not in covers:
        raise ValueError("pentad record is not an eligible cover")
    baseline = reference.analyze_cover(
        adj, outside, defects, zmask,
        reference.SupportSolver(), reference.ZeroForcingSolver(),
        reference.matching_size(defects), reference.CliqueStructureSolver(),
    )
    if inherited_cover_status(
        adj, outside, defects, seed, zmask, baseline,
        prior_failures, tetrad_failures,
    ) != "passing":
        raise ValueError("pentad certificate does not address a current cover")
    nvertices = tuple(
        outside[position] for position in range(len(outside))
        if not zmask & (1 << position)
    )
    graph_n = support_independent.independent_induced_graph(adj, nvertices)
    omega = reference.clique_number(graph_n)
    clique_mask = next(reference.clique_masks(graph_n, omega))
    return {
        "graph_index": int(graph["index"]),
        "seed": list(seed),
        "zmask": zmask,
        "nvertices": list(nvertices),
        "graph_n": list(graph_n),
        "rank_upper": int(baseline.k_rank_upper),
        "maximum_clique_size": omega,
        "fixed_maximum_clique": list(reference.bits(clique_mask)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_pentad_joint_conjunction_report.json",
    )
    parser.add_argument("--report-sha256")
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--pentad-workers", type=int, default=4)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_pentad_joint_conjunction_verification.json",
    )
    parser.add_argument("--pentad-worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.pentad_worker:
        payload = json.loads(sys.stdin.read())
        print(json.dumps(pentad_worker(tuple(payload)), sort_keys=True))
        return
    if not args.report_sha256:
        parser.error("--report-sha256 is required")
    if args.workers <= 0 or args.pentad_workers <= 0:
        parser.error("worker counts must be positive")

    started = time.monotonic()
    report_hash = sha256(args.report)
    if report_hash != args.report_sha256:
        raise ValueError("conjunction report differs from explicit hash pin")
    report = json.loads(args.report.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind") != "d6_k7_pentad_joint_seed_conjunction"
        or report.get("status") != "COMPLETE"
        or report.get("upstream_sha256") != dict(sorted(EXPECTED_HASHES.items()))
    ):
        raise ValueError("unexpected conjunction report schema/root hashes")
    builder_path = ROOT / "build_d6_k7_pentad_joint_conjunction.py"
    if report["source_sha256"][builder_path.name] != sha256(builder_path):
        raise ValueError("conjunction builder source binding mismatch")
    for name, expected in report["source_sha256"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"reported source hash mismatch: {name}")

    indices, graph_by_index, joint_rows, joint_report, pentad_report = (
        validate_roots()
    )
    selected = set(indices)
    prior = read_witness_keys(
        PRIOR_CERTIFICATES, selected, "dual_failure_witnesses"
    )
    tetrad = read_witness_keys(
        TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses"
    )
    tetrad_rows = read_tetrad_rows(selected)
    payloads = [
        (
            graph_by_index[index],
            [[list(seed), zmask] for seed, zmask in sorted(prior.get(index, set()))],
            [[list(seed), zmask] for seed, zmask in sorted(tetrad.get(index, set()))],
        )
        for index in indices
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        verified = list(executor.map(verify_graph, payloads, chunksize=1))
    if [int(item["index"]) for item in verified] != indices:
        raise ValueError("independent graph replay order mismatch")
    verified_by_index = {int(item["index"]): item for item in verified}

    graph_summaries = []
    totals: Counter = Counter()
    for item in verified:
        index = int(item["index"])
        counts = item["counts"]
        summary = {
            "index": index,
            "seeds": int(counts.get("seeds", 0)),
            "covers": int(counts.get("covers", 0)),
            "current_passing_covers": int(counts.get("current_passing_covers", 0)),
            "joint_pre_capacity_failing_covers": int(
                counts.get("joint_pre_capacity_failing_covers", 0)
            ),
            "pre_capacity_passing_covers": int(
                counts.get("pre_capacity_passing_covers", 0)
            ),
        }
        for name, value in summary.items():
            if name != "index":
                if int(joint_rows[index][name]) != value:
                    raise ValueError(f"graph {index} archive mismatch in {name}")
                totals[name] += value
        if int(tetrad_rows[index]["tetrad_passing_covers"]) != summary[
            "current_passing_covers"
        ]:
            raise ValueError("tetrad/current cover count mismatch")
        graph_summaries.append(summary)
    replay = report["joint_cover_quantifier_replay"]
    if (
        replay["ordered_graph_summaries"] != graph_summaries
        or replay["ordered_graph_summaries_sha256"] != stable_hash(graph_summaries)
        or replay["required_cover_partition"] != {
            "current_passing_covers": 738,
            "joint_pre_capacity_failing_covers": 342,
            "pre_capacity_passing_covers": 396,
        }
        or totals["current_passing_covers"] != 738
        or totals["joint_pre_capacity_failing_covers"] != 342
        or totals["pre_capacity_passing_covers"] != 396
    ):
        raise ValueError("global independent cover partition mismatch")

    pentad_by_key = {}
    rejected_order = []
    pentad_tasks = []
    for record in pentad_report["covers"]:
        if record["status"] != "REJECTED":
            continue
        key = key_tuple(record)
        if key in pentad_by_key:
            raise ValueError("duplicate pentad rejection key")
        pentad_by_key[key] = record
        rejected_order.append(key)
        target = build_pentad_target(
            graph_by_index[key[0]], record,
            prior.get(key[0], set()), tetrad.get(key[0], set()),
        )
        pentad_tasks.append((int(record["ordinal"]), target, record))
    if len(pentad_tasks) != 36:
        raise ValueError("pentad certificate population is not 36")
    with ThreadPoolExecutor(max_workers=args.pentad_workers) as executor:
        pentad_results = list(executor.map(pentad_subprocess, pentad_tasks))
    if any(
        result.get("status") != "REJECTED"
        or not result.get("certificate_checked")
        for result in pentad_results
    ):
        raise ValueError("an independently expanded pentad did not certify")

    expected_key_payloads = [key_payload(key) for key in rejected_order]
    pentad_section = report["pentad_cover_rejections"]
    if (
        pentad_section["covers"] != 36
        or pentad_section["graphs"] != 34
        or pentad_section["ordered_keys"] != expected_key_payloads
        or pentad_section["ordered_keys_sha256"]
        != stable_hash(expected_key_payloads)
    ):
        raise ValueError("reported pentad key list mismatch")

    new_section = report["new_seed_conjunction_rejections"]
    candidate_indices = [int(index) for index in new_section["ordered_indices"]]
    joint_rejected = [
        int(index) for index in joint_report["summary"]["joint_rejected_indices"]
    ]
    expected_candidates = [
        index for index in indices
        if any(key[0] == index for key in pentad_by_key)
        and index not in set(joint_rejected)
    ]
    joint_section = report["joint_support_rejections"]
    if (
        joint_section["count"] != 69
        or joint_section["ordered_indices"] != joint_rejected
        or joint_section["ordered_indices_sha256"] != stable_hash(joint_rejected)
        or candidate_indices != expected_candidates
        or new_section["count"] != 34
        or new_section["ordered_indices_sha256"] != stable_hash(candidate_indices)
    ):
        raise ValueError("reported prior/new rejection index lists do not bind")
    witnesses = []
    conservative_rejected = []
    full_rejected = []
    for index in candidate_indices:
        item = verified_by_index[index]
        conservative_seeds = []
        full_seeds = []
        for seed_record in item["seed_quantifiers"]:
            current_keys = [
                (index, tuple(cover["seed"]), int(cover["zmask"]))
                for cover in seed_record["current_passing_covers"]
            ]
            if not current_keys:
                continue
            covers = []
            for cover, key in zip(
                seed_record["current_passing_covers"], current_keys
            ):
                entry = {
                    **key_payload(key),
                    "joint_pre_capacity_status": cover[
                        "joint_pre_capacity_status"
                    ],
                    "family_counts": cover["family_counts"],
                }
                pentad_record = pentad_by_key.get(key)
                entry["pentad_certificate"] = None if pentad_record is None else {
                    "report_ordinal": int(pentad_record["ordinal"]),
                    "certificate_sha256": stable_hash(
                        pentad_record["certificate"]
                    ),
                    "pentad_index": int(
                        pentad_record["certificate"]["pentad_index"]
                    ),
                }
                covers.append(entry)
            seed_payload = {
                "seed": seed_record["seed"],
                "seed_mask": int(seed_record["seed_mask"]),
                "eligible_covers": int(seed_record["eligible_covers"]),
                "cover_status_counts": seed_record["cover_status_counts"],
                "joint_pre_capacity_failing_covers": int(
                    seed_record["joint_pre_capacity_failing_covers"]
                ),
                "pre_capacity_passing_covers": int(
                    seed_record["pre_capacity_passing_covers"]
                ),
                "current_passing_covers": covers,
            }
            if all(key in pentad_by_key for key in current_keys):
                conservative_seeds.append(seed_payload)
            if all(
                cover["joint_pre_capacity_status"] == "INFEASIBLE"
                or key in pentad_by_key
                for cover, key in zip(
                    seed_record["current_passing_covers"], current_keys
                )
            ):
                full_seeds.append(seed_payload)
        if conservative_seeds:
            conservative_rejected.append(index)
        if full_seeds:
            full_rejected.append(index)
        witnesses.append({
            "ordinal": indices.index(index),
            "index": index,
            "archived_graph_counts": {
                name: int(joint_rows[index][name])
                for name in (
                    "seeds", "covers", "current_passing_covers",
                    "joint_pre_capacity_failing_covers",
                    "pre_capacity_passing_covers",
                )
            },
            "conservative_pentad_only_seeds": conservative_seeds,
            "full_joint_or_pentad_seeds": full_seeds,
        })
    if (
        witnesses != new_section["witnesses"]
        or stable_hash(witnesses) != new_section["witnesses_sha256"]
        or conservative_rejected
        != new_section["conservative_pentad_only_ordered_indices"]
        or full_rejected != new_section["full_joint_or_pentad_ordered_indices"]
        or conservative_rejected != full_rejected
        or len(full_rejected) != 34
    ):
        raise ValueError("independent seed-local conjunction differs from report")

    if set(joint_rejected) & set(full_rejected):
        raise ValueError("new and prior joint rejection lists overlap")
    rejected_set = set(joint_rejected) | set(full_rejected)
    combined = [index for index in indices if index in rejected_set]
    residue = [index for index in indices if index not in rejected_set]
    if (
        len(combined) != 103 or len(residue) != 155
        or report["combined_rejections"]["count"] != 103
        or report["combined_rejections"]["ordered_indices"] != combined
        or report["combined_rejections"]["ordered_indices_sha256"]
        != stable_hash(combined)
        or report["exact_residue"]["count"] != 155
        or report["exact_residue"]["ordered_indices"] != residue
        or report["exact_residue"]["ordered_indices_sha256"] != stable_hash(residue)
    ):
        raise ValueError("independent 103/155 partition mismatch")

    positive = tuple(int(row) for row in lower_bound_18_graph())
    reference.validate_graph(positive)
    positive_k7 = sum(1 for _ in reference.clique_masks(positive, 7))
    if positive_k7 or report["positive_18_control"] != {
        "vertices": 18,
        "K7_seeds": 0,
        "status": "PASS_NOT_APPLICABLE_NO_K7",
    }:
        raise ValueError("known realizable 18-point control failed")

    output = {
        "schema": 1,
        "kind": "d6_k7_pentad_joint_seed_conjunction_verification",
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": report_hash},
        "root_artifact_sha256": dict(sorted(EXPECTED_HASHES.items())),
        "source_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            **EXPECTED_SOURCE_HASHES,
        },
        "checked": {
            "graphs": 258,
            "K7_seeds": totals["seeds"],
            "eligible_covers": totals["covers"],
            "inherited_passing_covers": totals["current_passing_covers"],
            "joint_pre_capacity_failing_covers": totals[
                "joint_pre_capacity_failing_covers"
            ],
            "pre_capacity_passing_covers": totals[
                "pre_capacity_passing_covers"
            ],
            "pentad_certificates_reexpanded": len(pentad_results),
            "prior_joint_rejections": len(joint_rejected),
            "new_seed_conjunction_rejections": len(full_rejected),
            "combined_rejections": len(combined),
            "exact_residue": len(residue),
        },
        "ordered_rejections_sha256": stable_hash(combined),
        "ordered_residue_sha256": stable_hash(residue),
        "positive_18_control": {
            "vertices": 18, "K7_seeds": 0,
            "status": "PASS_NOT_APPLICABLE_NO_K7",
        },
        "execution": {
            "workers": args.workers,
            "pentad_workers": args.pentad_workers,
            "logical_cpus": os.cpu_count(),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
        "claim": (
            "The complete 738-cover support partition, every seed-local OR "
            "quantifier, all 36 exact pentad identities, the disjoint 69+34 "
            "rejection union, and the ordered 155-graph residue were replayed "
            "without the conjunction builder or any ignored .runs input."
        ),
    }
    atomic_json(args.output, output)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "status": output["status"],
        "checked": output["checked"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
