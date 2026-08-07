#!/usr/bin/env python3
"""Exact fixed-maximum-clique basis-extension pilot on K7 residue covers.

For each selected no-near-clique cover, fix one deterministic maximum
required clique C.  If a feasible normalized Gram matrix K has rank r, then
K[C,C] is positive definite and C extends to an r-element principal basis.
Consequently it is enough to eliminate every C-subset basis of every size
|C| <= r <= U; arbitrary bases not containing C need not be enumerated.

HiGHS is only a locator for the restricted PSD-atom certificates.  Every
reported rejection has either a direct zero/duplicate-column witness or an
exact rational certificate checked by d6_k7_arbitrary_basis_psd_dual.
Unresolved bases and covers make no realizability claim.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import d6_k7_arbitrary_basis_psd_dual as dual
import d6_k7_rank_reference as prior


ROOT = Path(__file__).resolve().parent
EXPECTED_RANK_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_SELECTION_SHA256 = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
EXPECTED_PRIOR_CERTIFICATES_SHA256 = (
    "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
)
EXPECTED_UNION_SHA256 = (
    "1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599"
)


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


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def bits(mask: int) -> tuple[int, ...]:
    return tuple(prior.bits(mask))


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        ).stdout.strip()

    try:
        commit = git("rev-parse", "HEAD")
        branch = git("branch", "--show-current")
        status = git("status", "--porcelain=v1", "--untracked-files=all")
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}
    return {
        "available": True,
        "commit": commit,
        "branch": branch,
        "dirty": bool(status),
        "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
        "porcelain_lines": status.splitlines(),
    }


def load_campaign_inputs(
    rank_input_path: Path,
    selection_path: Path,
    union_path: Path,
) -> tuple[dict[int, dict], list[int], dict[int, list[dict]], dict]:
    """Load and hash-bind the frozen graph, selection, and residue inputs."""

    observed = {
        "rank_input": sha256(rank_input_path),
        "selection": sha256(selection_path),
        "union": sha256(union_path),
    }
    expected = {
        "rank_input": EXPECTED_RANK_INPUT_SHA256,
        "selection": EXPECTED_SELECTION_SHA256,
        "union": EXPECTED_UNION_SHA256,
    }
    if observed != expected:
        raise ValueError(f"frozen input hash mismatch: {observed} != {expected}")

    rank_payload = json.loads(rank_input_path.read_text(encoding="utf-8"))
    graphs = rank_payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != 17_764:
        raise ValueError("rank input graph count mismatch")
    graph_by_index: dict[int, dict] = {}
    for graph in graphs:
        index = int(graph["index"])
        if index in graph_by_index:
            raise ValueError("rank input repeats a graph index")
        adjacency = tuple(int(row) for row in graph["adjacency"])
        prior.validate_graph(adjacency)
        graph_by_index[index] = graph

    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if (
        selection.get("schema") != 1
        or selection.get("kind") != "d6_k7_rankone_tetrad_full_selection"
        or selection["rank_input"]["sha256"] != observed["rank_input"]
    ):
        raise ValueError("selection schema/input mismatch")
    selected_indices = [int(index) for index in selection["selected_indices"]]
    if (
        len(selected_indices) != 12_839
        or len(set(selected_indices)) != len(selected_indices)
        or selection["selected_indices_sha256"] != stable_hash(selected_indices)
    ):
        raise ValueError("selection index list mismatch")
    selected_set = set(selected_indices)

    certificate_path = selection_path.parent / selection[
        "prior_degree_one"
    ]["certificates"]["path"]
    certificate_hash = sha256(certificate_path)
    if certificate_hash != EXPECTED_PRIOR_CERTIFICATES_SHA256:
        raise ValueError("prior certificate archive hash mismatch")
    prior_witnesses: dict[int, list[dict]] = {}
    with gzip.open(certificate_path, "rt", encoding="utf-8") as stream:
        for line in stream:
            entry = json.loads(line)
            prior_witnesses[int(entry["index"])] = entry[
                "dual_failure_witnesses"
            ]

    union = json.loads(union_path.read_text(encoding="utf-8"))
    if (
        union.get("schema") != "d6-k7-rankone-pattern-union-v1"
        or union.get("status") != "COMPLETE"
    ):
        raise ValueError("union manifest schema/status mismatch")
    residue_set = union["sets"]["exact_residue"]
    residue = [int(index) for index in residue_set["indices"]]
    if residue_set["indices_sha256"] != stable_hash(residue):
        raise ValueError("union residue index hash mismatch")
    profiles = union["cover_structure"]["residue_profiles"]
    if [int(item["index"]) for item in profiles] != residue:
        raise ValueError("residue profiles are not in exact residue order")
    no_near_indices = [
        int(item["index"])
        for item in profiles if item["has_no_near_clique_cover"]
    ]
    if len(no_near_indices) != 50:
        raise ValueError("expected exactly 50 no-near residue graphs")
    if any(index not in selected_set for index in no_near_indices):
        raise ValueError("a no-near residue graph is absent from the selection")

    provenance = {
        "rank_input": {"path": str(rank_input_path), "sha256": observed["rank_input"]},
        "selection": {"path": str(selection_path), "sha256": observed["selection"]},
        "prior_certificates": {
            "path": str(certificate_path), "sha256": certificate_hash,
        },
        "union": {"path": str(union_path), "sha256": observed["union"]},
        "no_near_residue_graphs": 50,
        "no_near_residue_indices_sha256": stable_hash(no_near_indices),
    }
    return graph_by_index, no_near_indices, prior_witnesses, provenance


def extract_no_near_covers(
    graph_by_index: dict[int, dict],
    graph_indices: Sequence[int],
    prior_witnesses: dict[int, list[dict]],
) -> list[dict]:
    """Replay exact quantifiers and return all prior-passing no-near covers."""

    targets: list[dict] = []
    for graph_index in graph_indices:
        adjacency = tuple(int(row) for row in graph_by_index[graph_index]["adjacency"])
        witness_keys = {
            (tuple(int(vertex) for vertex in item["seed"]), int(item["zmask"]))
            for item in prior_witnesses.get(graph_index, [])
        }
        support_solver = prior.SupportSolver()
        zero_forcing = prior.ZeroForcingSolver()
        clique_solver = prior.CliqueStructureSolver()
        for seed_mask in prior.clique_masks(adjacency, 7):
            seed, outside, defects, ladj, eligible = prior.seed_instance(
                adjacency, seed_mask
            )
            total_term_rank = prior.matching_size(defects)
            for zmask in prior.eligible_covers(ladj, eligible):
                analysis = prior.analyze_cover(
                    adjacency, outside, defects, zmask,
                    support_solver, zero_forcing, total_term_rank, clique_solver,
                )
                if analysis.enhanced_joint_failed:
                    continue
                nvertices = [
                    outside[position] for position in range(len(outside))
                    if not (zmask & (1 << position))
                ]
                graph_n = prior.induced_graph(adjacency, nvertices)
                upper = int(analysis.k_rank_upper)
                if upper and next(prior.clique_masks(graph_n, upper), None) is not None:
                    continue
                if (tuple(seed), zmask) in witness_keys:
                    raise ValueError("no-saturating cover unexpectedly has a prior dual witness")
                near_size = upper - 1
                if near_size > 0 and next(
                    prior.clique_masks(graph_n, near_size), None
                ) is not None:
                    continue
                omega = prior.clique_number(graph_n)
                clique_mask = next(prior.clique_masks(graph_n, omega))
                clique = bits(clique_mask)
                targets.append({
                    "graph_index": int(graph_index),
                    "seed": [int(vertex) for vertex in seed],
                    "zmask": int(zmask),
                    "nvertices": [int(vertex) for vertex in nvertices],
                    "graph_n": [int(row) for row in graph_n],
                    "rank_upper": upper,
                    "maximum_clique_size": omega,
                    "fixed_maximum_clique": list(clique),
                })
    return targets


def extension_cores(target: dict) -> tuple[tuple[int, ...], ...]:
    """Enumerate exactly the fixed-clique principal-basis candidates."""

    clique = tuple(int(vertex) for vertex in target["fixed_maximum_clique"])
    graph_n = tuple(int(row) for row in target["graph_n"])
    remainder = tuple(vertex for vertex in range(len(graph_n)) if vertex not in clique)
    output = []
    for rank in range(len(clique), int(target["rank_upper"]) + 1):
        for extension in itertools.combinations(remainder, rank - len(clique)):
            output.append(tuple(sorted(clique + extension)))
    return tuple(output)


def easy_basis_certificate(system: dual.BasisAffineSystem) -> dict | None:
    """Return direct exact Schur contradictions visible in basis columns."""

    for index, column in enumerate(system.columns):
        if not any(column):
            return {
                "schema": 1,
                "kind": "zero_outside_column",
                "remainder_index": index,
                "vertex": system.remainder[index],
                "column": list(column),
            }
    first_by_column: dict[tuple[int, ...], int] = {}
    for index, column in enumerate(system.columns):
        if column in first_by_column:
            first = first_by_column[column]
            return {
                "schema": 1,
                "kind": "identical_outside_columns",
                "first_remainder_index": first,
                "second_remainder_index": index,
                "first_vertex": system.remainder[first],
                "second_vertex": system.remainder[index],
                "column": list(column),
            }
        first_by_column[column] = index
    return None


def verify_easy_basis_certificate(
    system: dual.BasisAffineSystem, certificate: dict
) -> None:
    """Check a direct basis-column contradiction without floating point."""

    kind = certificate.get("kind")
    if kind == "zero_outside_column":
        index = int(certificate["remainder_index"])
        if not 0 <= index < len(system.columns):
            raise ValueError("zero-column index out of range")
        if certificate.get("vertex") != system.remainder[index]:
            raise ValueError("zero-column vertex mismatch")
        if certificate.get("column") != list(system.columns[index]):
            raise ValueError("zero-column payload mismatch")
        if any(system.columns[index]):
            raise ValueError("claimed zero column is nonzero")
        return
    if kind == "identical_outside_columns":
        first = int(certificate["first_remainder_index"])
        second = int(certificate["second_remainder_index"])
        if not (0 <= first < second < len(system.columns)):
            raise ValueError("identical-column indices out of range/order")
        if certificate.get("first_vertex") != system.remainder[first]:
            raise ValueError("first identical-column vertex mismatch")
        if certificate.get("second_vertex") != system.remainder[second]:
            raise ValueError("second identical-column vertex mismatch")
        if certificate.get("column") != list(system.columns[first]):
            raise ValueError("identical-column payload mismatch")
        if system.columns[first] != system.columns[second]:
            raise ValueError("claimed identical columns differ")
        return
    raise ValueError(f"unknown easy certificate kind: {kind!r}")


def initialize_worker() -> None:
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def evaluate_basis(payload: tuple[int, list[int], list[int], int]) -> dict:
    cover_ordinal, graph_n_list, core_list, maximum_atom_support = payload
    started = time.perf_counter()
    graph_n = tuple(int(row) for row in graph_n_list)
    core = tuple(int(vertex) for vertex in core_list)
    system = dual.basis_affine_system(graph_n, core)
    easy = easy_basis_certificate(system)
    if easy is not None:
        verify_easy_basis_certificate(system, easy)
        outcome = "EASY_REJECTED"
        certificate = easy
    else:
        certificate = dual.find_certificate(system, maximum_atom_support)
        if certificate is None:
            outcome = "UNRESOLVED"
        else:
            dual.verify_certificate(system, certificate)
            outcome = "DUAL_REJECTED"
    return {
        "cover_ordinal": int(cover_ordinal),
        "rank": len(core),
        "core": list(core),
        "outcome": outcome,
        "certificate": certificate,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rank-input", type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    parser.add_argument(
        "--selection", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_selection.json"),
    )
    parser.add_argument(
        "--union", type=Path,
        default=Path("d6_k7_rankone_pattern_union.json"),
    )
    parser.add_argument("--covers", type=int, default=8)
    parser.add_argument(
        "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)
    )
    parser.add_argument("--maximum-atom-support", type=int, default=3)
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_arbitrary_basis_extension_pilot_report.json"),
    )
    args = parser.parse_args()
    if args.covers <= 0 or args.workers <= 0 or args.maximum_atom_support <= 0:
        parser.error("covers, workers, and maximum atom support must be positive")

    started_wall = time.perf_counter()
    started_at = utc_now()
    graph_by_index, graph_indices, witnesses, provenance = load_campaign_inputs(
        args.rank_input, args.selection, args.union
    )
    all_targets = extract_no_near_covers(
        graph_by_index, graph_indices, witnesses
    )
    if len(all_targets) != 53:
        raise ValueError(f"expected 53 exact no-near covers, found {len(all_targets)}")
    targets = all_targets[:args.covers]
    if len(targets) != args.covers:
        raise ValueError("requested more covers than the exact target population")

    jobs = []
    for ordinal, target in enumerate(targets):
        target["ordinal"] = ordinal
        cores = extension_cores(target)
        target["extension_core_count"] = len(cores)
        target["extension_cores_sha256"] = stable_hash([list(core) for core in cores])
        for core in cores:
            jobs.append((ordinal, target["graph_n"], list(core), args.maximum_atom_support))

    with ProcessPoolExecutor(
        max_workers=min(args.workers, len(jobs)), initializer=initialize_worker
    ) as pool:
        basis_records = list(pool.map(evaluate_basis, jobs, chunksize=1))
    basis_records.sort(key=lambda item: (
        int(item["cover_ordinal"]), int(item["rank"]), tuple(item["core"])
    ))

    by_cover: dict[int, list[dict]] = {ordinal: [] for ordinal in range(len(targets))}
    for record in basis_records:
        by_cover[int(record["cover_ordinal"])].append(record)
    for ordinal, target in enumerate(targets):
        records = by_cover[ordinal]
        counts = Counter(record["outcome"] for record in records)
        by_rank = {}
        for rank in range(
            target["maximum_clique_size"], target["rank_upper"] + 1
        ):
            rank_records = [record for record in records if record["rank"] == rank]
            by_rank[str(rank)] = dict(Counter(
                record["outcome"] for record in rank_records
            ))
        target["basis_outcome_counts"] = dict(counts)
        target["basis_outcome_counts_by_rank"] = by_rank
        target["cover_rejected"] = counts["UNRESOLVED"] == 0
        target["basis_records"] = records

    aggregate = Counter(record["outcome"] for record in basis_records)
    total_worker_seconds = sum(
        float(record["elapsed_seconds"]) for record in basis_records
    )
    source_paths = [
        Path(__file__).resolve(),
        ROOT / "d6_k7_arbitrary_basis_psd_dual.py",
        ROOT / "d6_k7_rank_reference.py",
    ]
    report = {
        "schema": 1,
        "kind": "d6_k7_arbitrary_basis_fixed_maximum_clique_extension_pilot",
        "status": "COMPLETE",
        "claim": (
            "Only bases carrying exact direct or rational PSD-atom certificates "
            "are rejected. UNRESOLVED is not a realizability claim. A cover is "
            "rejected only if the complete fixed-maximum-clique extension "
            "universe through its exact rank upper bound has no unresolved basis."
        ),
        "started_at": started_at,
        "finished_at": utc_now(),
        "configuration": {
            "selected_cover_prefix": args.covers,
            "workers": min(args.workers, len(jobs)),
            "logical_cpus": os.cpu_count(),
            "maximum_atom_support": args.maximum_atom_support,
            "command": " ".join(sys.argv),
        },
        "provenance": provenance,
        "git": git_provenance(),
        "machine": {
            "platform": platform.platform(),
            "python": sys.version,
            "implementation": platform.python_implementation(),
        },
        "source_hashes": {
            path.name: sha256(path) for path in source_paths
        },
        "population": {
            "no_near_residue_graphs": len(graph_indices),
            "no_near_covers": len(all_targets),
            "selected_covers": len(targets),
            "selected_cover_keys_sha256": stable_hash([
                [item["graph_index"], item["seed"], item["zmask"]]
                for item in targets
            ]),
        },
        "summary": {
            "basis_candidates": len(basis_records),
            "basis_outcome_counts": dict(aggregate),
            "covers_rejected": sum(bool(item["cover_rejected"]) for item in targets),
            "covers_unresolved": sum(not item["cover_rejected"] for item in targets),
            "total_worker_seconds": total_worker_seconds,
            "wall_seconds": time.perf_counter() - started_wall,
        },
        "covers": targets,
        "trust_scope": {
            "rejections": (
                "Exact integer graph parsing and Fraction certificate checks only; "
                "the floating-point LP is an untrusted certificate locator."
            ),
            "direct_zero_column": (
                "Schur equality would force K_yy=0, contradicting K_yy>1."
            ),
            "direct_identical_columns": (
                "Schur equality would force K_yy=K_yz in {0,1}, contradicting K_yy>1."
            ),
            "unresolved": "No mathematical conclusion.",
        },
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "summary": report["summary"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
