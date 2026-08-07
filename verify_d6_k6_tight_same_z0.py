#!/usr/bin/env python3
"""Independent verifier for the strongest K6 same-Z0 conjunction.

This checker imports neither the new producer nor its discovery probe.  It
combines the frozen independent repeated-arm/tight-Hall kernel with the older
independent Hall-subset/two-light-ray actual-support implementation, sharing
the identical independently enumerated Z0 in both systems.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import verify_d6_k6_repeated_two_support_arm as prior_verifier


strongest = prior_verifier.fan_verifier.parent
base = strongest.base
zf_reference = strongest.zf_reference
ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-tight-same-z0-v1"
CERTIFICATE_SCHEMA = "d6-k6-tight-same-z0-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-tight-same-z0-verification-v1"
EXPECTED_INPUT = 625
EXPECTED_INPUT_SHA256 = (
    "04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1"
)
EXPECTED_PARENT_RAW_CERTIFICATE_SHA256 = (
    "fbaebe286acce8b0c3603f89f05a5316010b51dceabfeface7e73ff55ce90ed2"
)
EXPECTED_PARENT_ARTIFACTS = {
    "d6_k6_repeated_two_support_arm.py": (
        "cf4a0a534504e312c47fcf1974d3362f19beeec270e4a49f04eba3080d06e1d4"
    ),
    "d6_k6_repeated_two_support_arm_report.json": (
        "ea8d9438d062b92857a1057b950e76b8ae08bbe1bf97b4327130fd4058171ebe"
    ),
    "d6_k6_repeated_two_support_arm_certificates.json.gz": (
        "604ec86b943ec909d177acc1fd36ef623c3e3bff67b4a8a5e5ffa601fc202af5"
    ),
    "d6_k6_repeated_two_support_arm_verification.json": (
        "3abf48bab092671a9dff03c2dc47c0063a75515cd219e1d5a4bfa9a5fa797b5b"
    ),
    "verify_d6_k6_repeated_two_support_arm.py": (
        "b23493bd62e79874117a2349462d66a12ef1577e720ff1aeb63224d14e2aecd7"
    ),
    "verify_d6_k6_same_z0.py": (
        "4db05ff06910a11395e4bbdd2fd13161d6c19f58f81cc99851b8c37558d7000f"
    ),
    "d6_k6_same_z0_verification.json": (
        "9d1f06c41827193f3a324bf2d08f450cf816567cbae727e21a6258afd68a9fbf"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gunzipped_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    prior_verifier.activate_kernel()


def assert_import_independence() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "d6_k6_tight_same_z0",
        "probe_d6_k6_tight_same_z0_conjunction",
        "d6_k6_support_reference",
        "d6_k6_bipartite_rank_reference",
    }
    if imports & forbidden:
        raise AssertionError(f"independent checker imports production code: {imports & forbidden}")


def load_input() -> tuple[list[dict], list[int]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_PARENT_ARTIFACTS
    }
    if observed != EXPECTED_PARENT_ARTIFACTS:
        raise ValueError(f"frozen parent/independent artifacts changed: {observed}")
    compressed = ROOT / "d6_k6_repeated_two_support_arm_certificates.json.gz"
    if gunzipped_sha256(compressed) != EXPECTED_PARENT_RAW_CERTIFICATE_SHA256:
        raise ValueError("repeated-arm compressed certificate payload changed")
    parent_report = json.loads(
        (ROOT / "d6_k6_repeated_two_support_arm_report.json").read_text(
            encoding="utf-8"
        )
    )
    parent_verification = json.loads(
        (ROOT / "d6_k6_repeated_two_support_arm_verification.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        parent_report.get("schema") != "d6-k6-repeated-two-support-arm-v1"
        or parent_report.get("status") != "COMPLETE"
        or parent_report.get("graphs_rejected") != 9
        or parent_report.get("graphs_surviving") != EXPECTED_INPUT
        or parent_report.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
        or parent_verification.get("schema")
        != "d6-k6-repeated-two-support-arm-verification-v1"
        or parent_verification.get("status") != "PASS"
        or parent_verification.get("graphs_rejected") != 9
        or parent_verification.get("graphs_surviving") != EXPECTED_INPUT
        or parent_verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen repeated-arm result boundary changed")
    parent_records, _ = prior_verifier.load_input()
    by_index = {int(record["index"]): record for record in parent_records}
    indices = [int(index) for index in parent_report["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("independent 625 input reconstruction changed")
    return [by_index[index] for index in indices], indices


def bipartite_states(adjacency, instance, z0, inertia, forcing, cache):
    states: set[tuple[int, ...]] = {()}
    local_components = base.components(instance, z0)
    for component in local_components:
        if not component.bipartite:
            continue
        options = strongest.component_options(
            adjacency, instance, z0, component,
            inertia, forcing, cache,
        )
        states = strongest.extend_states(states, options, adjacency)
        if not states:
            break
    return states, local_components


def nonbipartite_passes(instance, z0, local_components) -> bool:
    odd = [
        component.vertices
        for component in local_components
        if not component.bipartite
    ]
    return base.nonbipartite_system_passes(instance, z0, odd)


def cross_classification(adjacency, instance) -> dict:
    inertia = base.SympyInertiaCache()
    forcing = zf_reference.IndependentZeroForcing()
    cache = {}
    rows = []
    histogram: dict[str, int] = {}
    for z0 in base.z0_subsets(instance.eligible_z0):
        zvertices = [instance.outside[local] for local in base.bits(z0)]
        if not base.hall_matchable(z0, instance.defects):
            bipartite = False
            nonbipartite = False
        else:
            states, local_components = bipartite_states(
                adjacency, instance, z0, inertia, forcing, cache
            )
            bipartite = bool(states)
            nonbipartite = nonbipartite_passes(
                instance, z0, local_components
            )
        if bipartite and nonbipartite:
            category = "both"
        elif bipartite:
            category = "bipartite_only"
        elif nonbipartite:
            category = "nonbipartite_only"
        else:
            category = "neither"
        histogram[category] = histogram.get(category, 0) + 1
        rows.append({
            "Z0": zvertices,
            "tight_bipartite_passed": bipartite,
            "nonbipartite_actual_support_passed": nonbipartite,
            "category": category,
        })
    return {"histogram": dict(sorted(histogram.items())), "choices": rows}


def solve_seed(adjacency, instance) -> dict:
    inertia = base.SympyInertiaCache()
    forcing = zf_reference.IndependentZeroForcing()
    cache = {}
    failures = []
    for z0 in base.z0_subsets(instance.eligible_z0):
        zvertices = [instance.outside[local] for local in base.bits(z0)]
        if not base.hall_matchable(z0, instance.defects):
            failures.append({"Z0": zvertices, "reason": "Z0_matching"})
            continue
        states, local_components = bipartite_states(
            adjacency, instance, z0, inertia, forcing, cache
        )
        if not states:
            failures.append({
                "Z0": zvertices, "reason": "tight_bipartite_system"
            })
            continue
        if not nonbipartite_passes(instance, z0, local_components):
            failures.append({
                "Z0": zvertices,
                "reason": "nonbipartite_two_light_ray_actual_support",
            })
            continue
        return {"feasible": True, "failures": None}
    return {"feasible": False, "failures": failures}


def evaluate_record(record: dict) -> dict:
    prior_verifier.activate_kernel()
    adjacency = tuple(map(int, record["adjacency"]))
    base.validate_graph(adjacency)
    seeds_checked = 0
    for seed_mask in base.clique_masks(adjacency, base.COORDINATES):
        seeds_checked += 1
        instance = base.build_instance(adjacency, seed_mask)
        decision = solve_seed(adjacency, instance)
        if not decision["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "seed_mask": seed_mask,
                "seed": list(instance.seed),
                "failures": decision["failures"],
                "cross_classification": cross_classification(
                    adjacency, instance
                ),
            }
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "seed_mask": None,
        "seed": None,
        "failures": None,
        "cross_classification": None,
    }


def positive_control() -> dict:
    adjacency = base.lower_bound_18_graph()
    passed = 0
    for seed_mask in base.clique_masks(adjacency, base.COORDINATES):
        instance = base.build_instance(adjacency, seed_mask)
        if not solve_seed(adjacency, instance)["feasible"]:
            raise AssertionError("independent tight same-Z0 positive control failed")
        passed += 1
    if passed != 32:
        raise AssertionError("independent positive-control seed count changed")
    return {"passed": True, "K6_seeds": passed}


def verify(report_path: Path, certificates_path: Path, output: Path, workers: int) -> dict:
    assert_import_independence()
    records, indices = load_input()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    archive = json.loads(certificates_path.read_text(encoding="utf-8"))
    source_path = ROOT / "d6_k6_tight_same_z0.py"
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("production_source_sha256") != sha256(source_path)
        or report.get("dependencies")
        != {
            name: sha256(ROOT / name)
            for name in (
                "d6_k6_repeated_two_support_arm.py",
                "d6_k6_repeated_two_support_arm_report.json",
                "d6_k6_repeated_two_support_arm_certificates.json.gz",
                "d6_k6_repeated_two_support_arm_verification.json",
                "d6_k6_bipartite_rank_reference.py",
                "d6_k6_support_reference.py",
                "d6_k6_same_z0.py",
                "d6_k6_same_z0_report.json",
                "verify_d6_k6_same_z0.py",
                "d6_k6_same_z0_verification.json",
            )
        }
        or report.get("parent_raw_certificate_payload_sha256")
        != EXPECTED_PARENT_RAW_CERTIFICATE_SHA256
        or report.get("input_graphs") != EXPECTED_INPUT
        or report.get("ordered_input_indices") != indices
        or report.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(certificates_path)
        or report.get("positive_18_control")
        != {"passed": True, "K6_seeds": 32}
        or archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("production_source_sha256") != sha256(source_path)
        or archive.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("tight same-Z0 production boundary mismatch")
    started = time.monotonic()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [row["index"] for row in results if row["rejected"]]
    residue = [row["index"] for row in results if not row["rejected"]]
    if (
        rejected != report.get("rejected_indices")
        or stable_hash(rejected) != report.get("rejected_indices_sha256")
        or residue != report.get("ordered_residue_indices")
        or stable_hash(residue) != report.get("ordered_residue_indices_sha256")
        or len(rejected) != report.get("graphs_rejected")
        or len(residue) != report.get("graphs_surviving")
    ):
        raise AssertionError("independent tight same-Z0 decisions differ")
    saved_rows = {int(row["index"]): row for row in report["graph_results"]}
    certificates = {int(row["index"]): row for row in archive["rejected_graphs"]}
    if set(saved_rows) != set(indices) or len(saved_rows) != len(indices):
        raise AssertionError("production graph-result coverage differs")
    if (
        set(certificates) != set(rejected)
        or len(certificates) != len(archive["rejected_graphs"])
    ):
        raise AssertionError("certificate archive index coverage differs")
    for row in results:
        saved = saved_rows[row["index"]]
        if (
            row["rejected"] != saved["rejected"]
            or row["seeds_checked"] != saved["seeds_checked"]
            or row["seed_mask"] != saved["first_impossible_seed_mask"]
            or row["seed"] != saved["first_impossible_seed"]
        ):
            raise AssertionError(f"graph detail differs at {row['index']}")
        if row["rejected"]:
            certificate = certificates[row["index"]]
            if (
                certificate["seed_mask"] != row["seed_mask"]
                or certificate["seed"] != row["seed"]
                or certificate["Z0_failures"] != row["failures"]
                or certificate["same_Z0_cross_classification"]
                != row["cross_classification"]
            ):
                raise AssertionError(f"certificate differs at {row['index']}")
    control = positive_control()
    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": sha256(report_path)},
        "certificates": {
            "path": certificates_path.name,
            "sha256": sha256(certificates_path),
        },
        "sources": {
            Path(__file__).name: sha256(Path(__file__)),
            source_path.name: sha256(source_path),
            **{
                name: sha256(ROOT / name)
                for name in EXPECTED_PARENT_ARTIFACTS
            },
        },
        "checks": {
            "import_independence": True,
            "independent_repeated_arm_tight_Hall_kernel": True,
            "independent_Hall_subset_and_actual_support_DFS": True,
            "identical_Z0_quantifier": True,
            "all_graph_decisions": True,
            "all_rejected_Z0_partitions_and_cross_classifications": True,
            "positive_18_control": True,
        },
        "graphs_recomputed": len(results),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices_sha256": stable_hash(residue),
        "positive_18_control": control,
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report", type=Path, default=ROOT / "d6_k6_tight_same_z0_report.json"
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_tight_same_z0_certificates.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_tight_same_z0_verification.json",
    )
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    result = verify(args.report, args.certificates, args.output, args.workers)
    print(json.dumps({
        "status": result["status"],
        "rejected": result["graphs_rejected"],
        "surviving": result["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
