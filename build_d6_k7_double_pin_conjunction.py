#!/usr/bin/env python3
"""Build the exact two-defect double-pin conjunction on K7 residue v5."""

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
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation
import d6_k7_two_defect_double_pin as double_pin
import run_d6_k7_support_capacity_pilot as inherited


ROOT = Path(__file__).resolve().parent
RESIDUE = ROOT / "d6_current_residue_manifest_v5.json"
RESIDUE_VERIFICATION = ROOT / "d6_current_residue_manifest_v5_verification.json"
PRIOR_CERTIFICATES = ROOT / "d6_k7_positive_dual_full_certificates.jsonl.gz"
TETRAD_CERTIFICATES = ROOT / "d6_k7_rankone_tetrad_full_certificates.jsonl.gz"
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"

EXPECTED_HASHES = {
    "d6_current_residue_manifest_v5.json": (
        "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5"
    ),
    "d6_current_residue_manifest_v5_verification.json": (
        "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46"
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
    "run_d6_k7_support_capacity_pilot.py": (
        "9af542217d178bec2a71cb4c30faabc3dea8c8279053d729d924d9ca35eb128a"
    ),
    "d6_k7_support_propagation.py": (
        "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
    ),
    "d6_k7_small_support_value.py": (
        "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


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
            "branch": git("branch", "--show-current"),
            "commit": git("rev-parse", "HEAD"),
            "dirty": bool(status),
            "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
        }
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}


def read_failure_keys(path: Path, selected: set[int], field: str):
    output: dict[int, set[tuple[tuple[int, ...], int]]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            item = json.loads(line)
            index = int(item["index"])
            if index not in selected:
                continue
            for witness in item[field]:
                key = (tuple(map(int, witness["seed"])), int(witness["zmask"]))
                if key in output[index]:
                    raise ValueError("duplicate inherited cover certificate")
                output[index].add(key)
    return output


def read_tetrad_rows(selected: set[int]) -> dict[int, dict]:
    rows = {}
    with gzip.open(TETRAD_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            index = int(row["index"])
            if index in selected:
                rows[index] = row
    if set(rows) != selected:
        raise ValueError("v5 K7 graph missing from tetrad decisions")
    return rows


def evaluate_graph(payload: tuple) -> dict:
    graph, serialized_prior, serialized_tetrad, tetrad_row = payload
    prior = {(tuple(seed), int(zmask)) for seed, zmask in serialized_prior}
    tetrad = {(tuple(seed), int(zmask)) for seed, zmask in serialized_tetrad}
    index = int(graph["index"])
    adj = tuple(map(int, graph["adjacency"]))
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    totals: Counter = Counter()
    seed_records = []
    observed_current_keys = set()

    for seed_mask in reference.clique_masks(adj, 7):
        totals["seeds"] += 1
        seed, outside, defects, ladj, eligible = reference.seed_instance(adj, seed_mask)
        total_term_rank = reference.matching_size(defects)
        cover_records = []
        for zmask in reference.eligible_covers(ladj, eligible):
            totals["eligible_covers"] += 1
            baseline = reference.analyze_cover(
                adj, outside, defects, zmask, support_solver, zero_forcing,
                total_term_rank, clique_solver,
            )
            key = (tuple(seed), zmask)
            status = inherited.current_cover_status(
                adj, outside, defects, zmask, baseline, key, prior, tetrad
            )
            totals[f"inherited_cover_{status}"] += 1
            if status != "passing":
                continue
            observed_current_keys.add(key)
            totals["current_covers"] += 1
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex for vertex in range(len(outside)) if not zmask & (1 << vertex)
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nvertices]
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            family_counts: Counter = Counter()
            certificates = []
            passing_witness = None
            for z_supports in propagation.labeled_support_families(z_allowed):
                totals["labeled_z_families"] += 1
                family_counts["labeled_z_families"] += 1
                propagated = propagation.analyze_support_assignment(
                    graph_n, z_supports, n_allowed, zero_forcing, clique_solver
                )
                if propagated.failure is not None:
                    name = f"propagation_failure:{propagated.failure}"
                    totals[name] += 1
                    family_counts[name] += 1
                    continue
                sparse = sparse_value.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                if not sparse.feasible:
                    totals["sparse_value_failure"] += 1
                    family_counts["sparse_value_failure"] += 1
                    continue
                totals["pre_double_pin_families"] += 1
                family_counts["pre_double_pin_families"] += 1
                pin = double_pin.find_double_pin(
                    graph_n, propagated.propagated_masks
                )
                if pin is not None:
                    double_pin.verify_certificate(
                        graph_n, propagated.propagated_masks, pin
                    )
                    totals["double_pin_infeasible_families"] += 1
                    family_counts["double_pin_infeasible_families"] += 1
                    certificates.append({
                        "z_supports": list(z_supports),
                        "propagated_masks": list(propagated.propagated_masks),
                        "certificate": double_pin.certificate_json(pin),
                    })
                else:
                    totals["double_pin_passing_families"] += 1
                    family_counts["double_pin_passing_families"] += 1
                    if passing_witness is None:
                        passing_witness = {
                            "z_supports": list(z_supports),
                            "propagated_masks": list(propagated.propagated_masks),
                        }
            cover_status = "PASSING" if passing_witness is not None else "INFEASIBLE"
            totals[f"double_pin_{cover_status.lower()}_covers"] += 1
            cover_records.append({
                "seed": seed,
                "zmask": zmask,
                "status": cover_status,
                "family_counts": dict(sorted(family_counts.items())),
                "double_pin_certificates": certificates,
                "first_passing_witness": passing_witness,
            })
        if not cover_records:
            raise ValueError("v5 survivor unexpectedly has no inherited cover")
        seed_status = (
            "PASSING" if any(cover["status"] == "PASSING" for cover in cover_records)
            else "INFEASIBLE"
        )
        seed_records.append({
            "seed": seed,
            "seed_mask": seed_mask,
            "status": seed_status,
            "current_covers": cover_records,
        })

    if len(observed_current_keys) != int(tetrad_row["tetrad_passing_covers"]):
        raise ValueError("inherited current-cover reconstruction mismatch")
    rejecting = next(
        (seed for seed in seed_records if seed["status"] == "INFEASIBLE"), None
    )
    return {
        "index": index,
        "decision": "REJECTED" if rejecting is not None else "SURVIVOR",
        "counts": dict(sorted(totals.items())),
        "first_rejecting_seed_mask": (
            int(rejecting["seed_mask"]) if rejecting is not None else 0
        ),
        "seeds": seed_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_double_pin_conjunction_report.json",
    )
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    for name, expected in {**EXPECTED_HASHES, **EXPECTED_SOURCE_HASHES}.items():
        observed = sha256(ROOT / name)
        if observed != expected:
            raise SystemExit(f"hash mismatch for {name}: {observed} != {expected}")
    verification = json.loads(RESIDUE_VERIFICATION.read_text(encoding="utf-8"))
    if verification.get("status") != "PASS":
        raise SystemExit("v5 residue verification is not PASS")
    residue = json.loads(RESIDUE.read_text(encoding="utf-8"))
    graphs = residue["classes"]["K7"]["graphs"]
    indices = [int(graph["index"]) for graph in graphs]
    k7 = residue["classes"]["K7"]
    if (
        len(graphs) != 155
        or indices != k7["indices"]
        or stable_hash(indices) != k7["indices_sha256"]
        or stable_hash(graphs) != k7["graphs_sha256"]
    ):
        raise SystemExit("v5 K7 corpus does not bind")
    selected = set(indices)
    prior = read_failure_keys(PRIOR_CERTIFICATES, selected, "dual_failure_witnesses")
    tetrad = read_failure_keys(TETRAD_CERTIFICATES, selected, "tetrad_failure_witnesses")
    tetrad_rows = read_tetrad_rows(selected)
    payloads = [(
        graph,
        [[list(seed), zmask] for seed, zmask in sorted(prior.get(index, set()))],
        [[list(seed), zmask] for seed, zmask in sorted(tetrad.get(index, set()))],
        tetrad_rows[index],
    ) for graph, index in zip(graphs, indices, strict=True)]

    started_at = datetime.now(UTC).isoformat()
    started = time.monotonic()
    if args.workers == 1:
        records = list(map(evaluate_graph, payloads))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            records = list(executor.map(evaluate_graph, payloads, chunksize=1))
    if [record["index"] for record in records] != indices:
        raise SystemExit("parallel evaluation changed graph order")
    rejected = [record["index"] for record in records if record["decision"] == "REJECTED"]
    survivors = [record["index"] for record in records if record["decision"] == "SURVIVOR"]
    totals: Counter = Counter()
    for record in records:
        totals.update(record["counts"])
    if len(rejected) != 130 or len(survivors) != 25:
        raise SystemExit(
            f"audited expected split changed: {len(rejected)} rejected, "
            f"{len(survivors)} survivors"
        )
    source_hashes = {
        Path(__file__).name: sha256(Path(__file__)),
        "d6_k7_two_defect_double_pin.py": sha256(
            ROOT / "d6_k7_two_defect_double_pin.py"
        ),
        **EXPECTED_SOURCE_HASHES,
    }
    report = {
        "schema": 1,
        "kind": "d6_k7_two_defect_double_pin_seed_conjunction",
        "status": "COMPLETE",
        "claim": (
            "The exact two-defect double-pin obstruction rejects 130 of the "
            "155 K7-containing v5 residue graphs and leaves 25 unresolved."
        ),
        "semantics": {
            "candidate_nonedges_optional": True,
            "all_vertices_distinct": True,
            "only_required_edges_enter_pin_triangles": True,
            "all_seed_cover_support_quantifiers_exhausted": True,
            "survivor_is_not_realizable_claim": True,
            "floating_point_enters_rejection": False,
        },
        "input": {
            "graphs": len(graphs),
            "indices_sha256": stable_hash(indices),
            "graphs_sha256": stable_hash(graphs),
            "residue_file_sha256": EXPECTED_HASHES[
                "d6_current_residue_manifest_v5.json"
            ],
        },
        "upstream_sha256": dict(sorted(EXPECTED_HASHES.items())),
        "source_sha256": dict(sorted(source_hashes.items())),
        "summary": {
            "rejected": len(rejected),
            "rejected_indices": rejected,
            "rejected_indices_sha256": stable_hash(rejected),
            "survivors": len(survivors),
            "survivor_indices": survivors,
            "survivor_indices_sha256": stable_hash(survivors),
            "totals": dict(sorted(totals.items())),
        },
        "records": records,
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "workers": args.workers,
            "logical_cpus": os.cpu_count(),
            "started_at": started_at,
            "finished_at": datetime.now(UTC).isoformat(),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
        "git": git_provenance(),
    }
    atomic_json(args.output, report)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "rejected": len(rejected),
        "survivors": len(survivors),
        "totals": dict(sorted(totals.items())),
        "elapsed_seconds": report["execution"]["elapsed_seconds"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
