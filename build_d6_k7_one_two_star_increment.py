#!/usr/bin/env python3
"""Build the exact one-free-neighbour/two-free-center K7 increment."""

from __future__ import annotations

import argparse
import hashlib
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

import d6_k7_correlated_one_free_edge as correlated
import d6_k7_full_pin_odd_cycle as full_pin
import d6_k7_one_free_edge as one_free
import d6_k7_one_two_star as star
import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation
import d6_k7_two_defect_double_pin as double_pin
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "d6_current_residue_manifest_v6.json"
MANIFEST_VERIFICATION = ROOT / "d6_current_residue_manifest_v6_verification.json"
CURRENT = ROOT / "d6_k7_one_free_conjunction_report.json"
CURRENT_VERIFICATION = ROOT / "d6_k7_one_free_conjunction_verification.json"
EXPECTED_UPSTREAM = {
    MANIFEST.name: (
        "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
    ),
    MANIFEST_VERIFICATION.name: (
        "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
    ),
    CURRENT.name: (
        "181f63015d373fb69f37f4f390cdd8ad50bc32b5a2d97a36fd0bac876d0877a2"
    ),
    CURRENT_VERIFICATION.name: (
        "ebaacfe663728eff427c9ba2280fc4261d26dfe10fad21e1a00982b9d31ac062"
    ),
}
EXPECTED_INPUT_SHA256 = (
    "af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290"
)
EXPECTED_REJECTIONS = [2592657, 3785980, 3888410]
SOURCE_FILES = (
    "build_d6_k7_one_two_star_increment.py",
    "d6_k7_one_two_star.py",
    "d6_k7_one_free_edge.py",
    "d6_k7_correlated_one_free_edge.py",
    "d6_k7_two_defect_double_pin.py",
    "d6_k7_full_pin_odd_cycle.py",
    "d6_k7_support_propagation.py",
    "d6_k7_small_support_value.py",
    "d6_k7_rank_reference.py",
    "verify_profile_d6.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("ascii")
    ).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
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
        "porcelain_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "porcelain_lines": status.splitlines(),
    }


def validate_upstream() -> tuple[dict, dict]:
    paths = {
        MANIFEST.name: MANIFEST,
        MANIFEST_VERIFICATION.name: MANIFEST_VERIFICATION,
        CURRENT.name: CURRENT,
        CURRENT_VERIFICATION.name: CURRENT_VERIFICATION,
    }
    for name, expected in EXPECTED_UPSTREAM.items():
        if sha256(paths[name]) != expected:
            raise ValueError(f"upstream hash mismatch: {name}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest_verification = json.loads(
        MANIFEST_VERIFICATION.read_text(encoding="utf-8")
    )
    current = json.loads(CURRENT.read_text(encoding="utf-8"))
    current_verification = json.loads(
        CURRENT_VERIFICATION.read_text(encoding="utf-8")
    )
    if (
        manifest.get("schema") != "d6-current-certified-residue-v6"
        or manifest.get("status") != "COMPLETE_MIXED_CERTIFICATE_UNION"
        or manifest_verification.get("status") != "PASS"
        or manifest_verification.get("manifest", {}).get("sha256")
        != EXPECTED_UPSTREAM[MANIFEST.name]
        or current.get("kind") != "d6_k7_one_free_edge_seed_conjunction"
        or current.get("status") != "COMPLETE"
        or current_verification.get("status") != "PASS"
        or current_verification.get("report", {}).get("sha256")
        != EXPECTED_UPSTREAM[CURRENT.name]
        or current.get("semantics", {}).get("candidate_nonedges_optional")
        is not True
        or current.get("semantics", {}).get("floating_point_enters_rejection")
        is not False
    ):
        raise ValueError("upstream semantic or verification gate failed")
    return manifest, current


def evaluate_graph(payload: tuple[dict, dict]) -> dict:
    graph, current_record = payload
    index = int(graph["index"])
    adjacency = tuple(map(int, graph["adjacency"]))
    reference.validate_graph(adjacency)
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    counts: Counter[str] = Counter()
    certificates = []
    seed_records = []

    for archived_seed in current_record["seeds"]:
        seed_mask = int(archived_seed["seed_mask"])
        seed, outside, defects, _ladj, _eligible = reference.seed_instance(
            adjacency, seed_mask
        )
        if seed != archived_seed["seed"]:
            raise ValueError("seed reconstruction mismatch")
        cover_records = []
        for archived_cover in archived_seed["current_covers"]:
            zmask = int(archived_cover["zmask"])
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex
                for vertex in range(len(outside))
                if not zmask & (1 << vertex)
            )
            graph_n = reference.induced_graph(
                adjacency, [outside[vertex] for vertex in nvertices]
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            current_passing = 0
            star_passing = 0
            first_passing_witness = None
            cover_certificates = []
            for z_supports in propagation.labeled_support_families(z_allowed):
                counts["labeled_support_families"] += 1
                propagated = propagation.analyze_support_assignment(
                    graph_n,
                    z_supports,
                    n_allowed,
                    zero_forcing,
                    clique_solver,
                )
                if propagated.failure is not None:
                    counts[f"propagation_failure:{propagated.failure}"] += 1
                    continue
                masks = propagated.propagated_masks
                if not sparse_value.check_small_support_masks(graph_n, masks).feasible:
                    counts["sparse_value_failure"] += 1
                    continue
                if double_pin.find_double_pin(graph_n, masks) is not None:
                    counts["double_pin_failure"] += 1
                    continue
                if full_pin.find_full_pin(graph_n, masks) is not None:
                    counts["full_pin_failure"] += 1
                    continue
                if one_free.find_one_free_edge(graph_n, masks) is not None:
                    counts["singleton_one_free_failure"] += 1
                    continue
                if correlated.find_correlated_one_free_edge(graph_n, masks) is not None:
                    counts["correlated_one_free_failure"] += 1
                    continue
                current_passing += 1
                counts["current_passing_families"] += 1
                certificate = star.find_one_two_star(graph_n, masks)
                if certificate is None:
                    star_passing += 1
                    counts["star_passing_families"] += 1
                    if first_passing_witness is None:
                        first_passing_witness = {
                            "z_supports": list(z_supports),
                            "propagated_masks": list(masks),
                        }
                    continue
                star.verify_certificate(graph_n, masks, certificate)
                counts["star_failed_families"] += 1
                witness = {
                    "seed": seed,
                    "seed_mask": seed_mask,
                    "zmask": zmask,
                    "z_supports": list(z_supports),
                    "propagated_masks": list(masks),
                    "certificate": star.certificate_json(certificate),
                }
                certificates.append(witness)
                cover_certificates.append(witness)

            expected_current = int(
                archived_cover["family_counts"].get(
                    "combined_passing_families", 0
                )
            )
            if current_passing != expected_current:
                raise ValueError(
                    f"current family replay mismatch {index}/{seed_mask}/{zmask}: "
                    f"{current_passing} != {expected_current}"
                )
            cover_records.append(
                {
                    "zmask": zmask,
                    "current_passing_families": current_passing,
                    "star_failed_families": len(cover_certificates),
                    "star_passing_families": star_passing,
                    "status": "PASSING" if star_passing else "INFEASIBLE",
                    "first_passing_witness": first_passing_witness,
                }
            )
        if not any(cover["current_passing_families"] for cover in cover_records):
            raise ValueError("input graph was already rejected at a current seed")
        seed_records.append(
            {
                "seed": seed,
                "seed_mask": seed_mask,
                "status": (
                    "PASSING"
                    if any(cover["status"] == "PASSING" for cover in cover_records)
                    else "INFEASIBLE"
                ),
                "current_covers": cover_records,
            }
        )

    rejecting_seed = next(
        (record for record in seed_records if record["status"] == "INFEASIBLE"),
        None,
    )
    return {
        "index": index,
        "decision": "REJECTED" if rejecting_seed is not None else "SURVIVOR",
        "first_rejecting_seed_mask": (
            0 if rejecting_seed is None else int(rejecting_seed["seed_mask"])
        ),
        "counts": dict(sorted(counts.items())),
        "certificates": certificates,
        "seeds": seed_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_one_two_star_increment_report.json",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    manifest, current = validate_upstream()
    graphs = manifest["classes"]["K7"]["graphs"]
    indices = [int(graph["index"]) for graph in graphs]
    if len(indices) != 19 or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("v6 K7 input is not the frozen 19-list")
    current_by_index = {
        int(record["index"]): record for record in current["records"]
    }
    if any(index not in current_by_index for index in indices):
        raise ValueError("v6 K7 graph absent from current exact report")
    payloads = [
        (graph, current_by_index[int(graph["index"])]) for graph in graphs
    ]

    started_utc = datetime.now(UTC).isoformat()
    started = time.monotonic()
    if args.workers == 1:
        records = list(map(evaluate_graph, payloads))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            records = list(executor.map(evaluate_graph, payloads, chunksize=1))
    if [record["index"] for record in records] != indices:
        raise ValueError("worker output order changed")
    rejected = [
        record["index"] for record in records if record["decision"] == "REJECTED"
    ]
    if rejected != EXPECTED_REJECTIONS:
        raise ValueError(f"unexpected exact star marginal: {rejected}")
    rejected_set = set(rejected)
    survivors = [index for index in indices if index not in rejected_set]
    totals: Counter[str] = Counter()
    for record in records:
        totals.update(record["counts"])
    if (
        totals["current_passing_families"] != 88
        or totals["star_failed_families"] != 4
    ):
        raise ValueError("unexpected family profile")

    positive = tuple(map(int, lower_bound_18_graph()))
    reference.validate_graph(positive)
    positive_k7 = sum(1 for _ in reference.clique_masks(positive, 7))
    if positive_k7:
        raise ValueError("known realizable 18-point graph unexpectedly has K7")
    optional_zero = star.lines_hyperbola_decision(
        ((star.ONE, star.ZERO, star.ZERO),),
        star.ONE,
        star.q(-1),
    )
    if not optional_zero.feasible:
        raise ValueError("optional-zero line/hyperbola positive control failed")

    report = {
        "schema": 1,
        "kind": "d6_k7_one_free_neighbour_two_free_center_increment",
        "status": "COMPLETE",
        "claim": (
            "Every current propagated family on the frozen 19-graph K7 "
            "residue was replayed. Three graphs have a K7 seed for which "
            "every cover/family fails the exact correlated-sign affine-line "
            "and two-free-diagonal star system."
        ),
        "semantics": {
            "candidate_nonedges_optional": True,
            "allowed_unpinned_coordinates_may_be_zero": True,
            "only_required_edges_enter_star_equations": True,
            "propagated_masks_are_support_supersets": True,
            "normalization_factor_t_nonzero_on_cover_complement_N": True,
            "positive_sqrt7_embedding_checked_exactly": True,
            "floating_point_enters_rejection": False,
            "graph_rejected_if_any_k7_seed_is_infeasible": True,
            "cap500000_campaign_used": False,
        },
        "upstream_artifact_sha256": dict(sorted(EXPECTED_UPSTREAM.items())),
        "source_sha256": {
            name: sha256(ROOT / name) for name in SOURCE_FILES
        },
        "input": {
            "graphs": len(indices),
            "ordered_indices": indices,
            "ordered_indices_sha256": stable_hash(indices),
            "embedded_graphs_sha256": stable_hash(graphs),
        },
        "summary": {
            "totals": dict(sorted(totals.items())),
            "graphs_rejected": len(rejected),
            "rejected_indices": rejected,
            "rejected_indices_sha256": stable_hash(rejected),
            "graphs_surviving": len(survivors),
            "ordered_survivor_indices": survivors,
            "ordered_survivor_indices_sha256": stable_hash(survivors),
        },
        "records": records,
        "positive_controls": {
            "known_realizable_18": {
                "vertices": 18,
                "K7_seeds": 0,
                "status": "PASS_NOT_APPLICABLE_NO_K7",
            },
            "optional_zero_line_hyperbola": {
                "identity": "(x-1)(y-1)=-1 with required line x=0",
                "witness": {"x": 0, "y": 2},
                "classification": optional_zero.reason,
                "status": "PASS_FEASIBLE",
            },
        },
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "workers": args.workers,
            "logical_cpus": os.cpu_count(),
            "started_utc": started_utc,
            "finished_utc": datetime.now(UTC).isoformat(),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "git": git_provenance(),
        },
    }
    atomic_json(args.output.resolve(), report)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
                "graphs": len(indices),
                "current_passing_families": totals["current_passing_families"],
                "star_failed_families": totals["star_failed_families"],
                "rejected_indices": rejected,
                "survivors": len(survivors),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
