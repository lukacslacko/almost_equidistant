#!/usr/bin/env python3
"""Independent checker for ``d6_current_residue_manifest.json``.

The checker does not import the manifest builder.  It separately parses every
historical rejection list needed for overlap accounting and both new full
decision archives, then reconstructs all class-final and combined index lists.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import time
from itertools import combinations
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
EXPECTED = {
    ".runs/d6_k7_rank_survivors.json": "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9",
    "build_d6_k7_positive_dual_selection.py": "5717d24d448c1dc8b1da7dbb8dca53902d9b40546e9f4edc70c3183fde995dad",
    "d6_k7_positive_dual_selection.json": "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0",
    "d6_k7_support_rank_survivors_report.json": "5138d207543967b393baae9b14b47b49a276273cbc4468ee843856891917691e",
    "d6_k7_strict_h_full_report.json": "59be9c6a2d4cd2e8e8214d28d3198487846e58697046415ae1e15f4f69dc3add",
    "d6_k7_small_support_value_full_report.json": "cea64f3dde804e0766c5b77e373aa50338d5bee6e5e54f44749f4e387cf529aa",
    "d6_interval_residue_benchmark.json": "e246c3896e7ff2b9b40a598efdaffcc465ff005700a0d51fb17b66b2d6f68481",
    "d6_k7_positive_polynomial_dual_report.json": "51f29654f9514eec8f3911149a42407a7f197470720061528064ab77c9dec120",
    "d6_k7_positive_dual_full_report.json": "c2de7e06bbe4f3d163a4f747d42721f7f60da985b867ce74669664431738e345",
    "d6_k7_positive_dual_full_verification.json": "90df584f7a9b15438c3c0dd7ca04f49a586df6e73031da1cc55f8325b5ef950e",
    "d6_k7_positive_dual_full_decisions.tsv.gz": "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec",
    "d6_k7_positive_dual_full_certificates.jsonl.gz": "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f",
    "d6_k7_positive_dual_full_checkpoint.json": "971d1d38db145a49a9dc5cec61f9c2bf9ba17eccee9a26a90f59af49f6777ec4",
    "d6_k7_reflection_overlap_full_report.json": "d8294f1d1dd11206856ace4bfe210cbc79f4bb019969930645d95a7bad89cd07",
    "d6_k6_bipartite_rank_input.json": "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845",
    "d6_k6_bipartite_rank_report.json": "ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc",
    "d6_k6_normal_inertia_full_report.json": "8fd3d9d10038040f42407b25ea007f59f0e84e6eb0da8bec2694d1831e1991c5",
    "d6_k6_normal_inertia_full_verification_report.json": "61889185deb0e773bf899cdbb82c96467cc38a4688560a50f53c7761a7615e6b",
    "d6_k6_normal_inertia_full_decisions.jsonl.gz": "614bbb64989b5834ed0dbe1c12e7718583e53cfd68483668ffbbad44608b3828",
    "d6_k6_normal_inertia_full_checkpoint.json": "f4cb37faab90ed65b08b0e0ae2015a226de577caad5c74d23b530377d089962d",
    "build_d6_current_residue_manifest.py": "256d5341bdd1befa2007c66564f374c7785b554ee57b295b0855dacc5e8ef202",
    "d6_current_residue_manifest.json": "d06cdb23b03b1f8523fbe3c3ede5bc297ace21a844282830c5cf51a6cee3512d",
}
HISTORICAL = (
    ("support_propagation", "d6_k7_support_rank_survivors_report.json", 1536),
    ("strict_H", "d6_k7_strict_h_full_report.json", 603),
    ("sparse_value", "d6_k7_small_support_value_full_report.json", 3195),
    ("interval_cap20000", "d6_interval_residue_benchmark.json", 10),
    ("sample_positive_dual", "d6_k7_positive_polynomial_dual_report.json", 2),
)


def path(name: str) -> Path:
    return ROOT / name


def sha256(file_path: Path) -> str:
    digest = hashlib.sha256()
    with file_path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("ascii")).hexdigest()


def load(name: str) -> dict:
    value = json.loads(path(name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{name} is not an object")
    return value


def as_set(value: object, label: str) -> set[int]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    result = set(value)
    if len(result) != len(value):
        raise ValueError(f"{label} repeats an index")
    return result


def historical_sets() -> dict[str, set[int]]:
    support = load("d6_k7_support_rank_survivors_report.json")
    strict = load("d6_k7_strict_h_full_report.json")
    sparse = load("d6_k7_small_support_value_full_report.json")
    interval = load("d6_interval_residue_benchmark.json")
    sample = load("d6_k7_positive_polynomial_dual_report.json")
    caps = [row for row in interval["benchmarks"] if row["cap"] == 20_000]
    if len(caps) != 1:
        raise ValueError("interval report has no unique cap-20,000 row")
    answer = {
        "support_propagation": as_set(support["refined_rejected_indices"], "support"),
        "strict_H": as_set(strict["decisions"]["strict_H_rejected"], "strict H"),
        "sparse_value": as_set(sparse["decisions"]["rejected_indices"], "sparse"),
        "interval_cap20000": as_set(caps[0]["killed_indices"], "interval"),
        "sample_positive_dual": as_set(
            sample["decisions"]["marginal_dual_rejected"], "sample dual"
        ),
    }
    for name, _, expected in HISTORICAL:
        if len(answer[name]) != expected:
            raise ValueError(f"historical raw count mismatch at {name}")
    return answer


def parse_k7_dual(base: Sequence[int]) -> tuple[set[int], str, dict[str, int]]:
    with gzip.open(path("d6_k7_positive_dual_full_decisions.tsv.gz"), "rb") as stream:
        raw = stream.read()
    raw_hash = hashlib.sha256(raw).hexdigest()
    reader = csv.DictReader(io.StringIO(raw.decode("ascii")), delimiter="\t")
    rejected = set()
    statuses: dict[str, int] = {}
    rows = list(reader)
    if len(rows) != len(base):
        raise ValueError("K7 archive length mismatch")
    for ordinal, (row, index) in enumerate(zip(rows, base)):
        if int(row["ordinal"]) != ordinal or int(row["index"]) != index:
            raise ValueError(f"K7 archive order mismatch at {ordinal}")
        statuses[row["status"]] = statuses.get(row["status"], 0) + 1
        if row["error_type"]:
            raise ValueError("K7 archive has infrastructure errors")
        if int(row["marginal_dual_rejected"]):
            rejected.add(index)
    return rejected, raw_hash, statuses


def parse_k6_normal(base: Sequence[int]) -> tuple[set[int], str]:
    with gzip.open(path("d6_k6_normal_inertia_full_decisions.jsonl.gz"), "rb") as stream:
        raw = stream.read()
    rows = [json.loads(line) for line in raw.splitlines()]
    if len(rows) != len(base):
        raise ValueError("K6 archive length mismatch")
    rejected = set()
    for position, (row, index) in enumerate(zip(rows, base)):
        if row["index"] != index or row["schema"] != "d6-k6-normal-inertia-full-decision-v1":
            raise ValueError(f"K6 archive order mismatch at {position}")
        if row["decision"]["rejected"]:
            rejected.add(index)
    return rejected, hashlib.sha256(raw).hexdigest()


def layer(
    name: str,
    population: str,
    source: str,
    raw: set[int],
    universe: set[int],
    prior: set[int],
    reported: int,
    cumulative: set[int],
    include: bool,
) -> dict:
    inside = raw & universe
    overlap = inside & prior
    incremental = inside - prior
    record = {
        "layer": name,
        "population": population,
        "source": source,
        "reported_rejections": reported,
        "raw_rejections": len(raw),
        "raw_rejections_sha256": stable_hash(sorted(raw)),
        "rejections_in_universe": len(inside),
        "rejections_in_universe_sha256": stable_hash(sorted(inside)),
        "rejections_outside_universe": len(raw - universe),
        "overlap_with_prior_layers": len(overlap),
        "overlap_with_prior_layers_sha256": stable_hash(sorted(overlap)),
        "incremental_rejections": len(incremental),
        "incremental_rejections_sha256": stable_hash(sorted(incremental)),
        "cumulative_rejections": len(cumulative),
        "cumulative_survivors": len(universe - cumulative),
    }
    if include:
        record["incremental_rejected_indices"] = sorted(incremental)
    return record


def overlap_rows(layer_sets: dict[str, set[int]]) -> list[dict]:
    rows = []
    for left, right in combinations(layer_sets, 2):
        intersection = layer_sets[left] & layer_sets[right]
        rows.append({
            "left": left,
            "right": right,
            "count": len(intersection),
            "indices_sha256": stable_hash(sorted(intersection)),
        })
    return rows


def verify(manifest_path: Path) -> dict:
    observed = {name: sha256(path(name)) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError("one or more pinned source hashes changed")
    if manifest_path.resolve() != path("d6_current_residue_manifest.json").resolve():
        raise ValueError("checker accepts only the hash-pinned manifest path")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "d6-current-exact-residue-v1" or manifest.get("status") != "COMPLETE":
        raise ValueError("manifest schema/status mismatch")

    rank = load(".runs/d6_k7_rank_survivors.json")
    rank_indices = [int(row["index"]) for row in rank["graphs"]]
    rank_set = set(rank_indices)
    frozen = load("d6_k7_positive_dual_selection.json")
    k7_base = list(frozen["selected_indices"])
    if len(rank_set) != 17_764 or len(k7_base) != 12_941 or stable_hash(k7_base) != frozen["selected_indices_sha256"]:
        raise ValueError("frozen K7 base mismatch")
    historical = historical_sets()
    k7_seen = set()
    expected_layers = []
    layer_sets: dict[str, set[int]] = {}
    source_names = {name: source for name, source, _ in HISTORICAL}
    reported_counts = {name: count for name, _, count in HISTORICAL}
    for name in historical:
        prior = set(k7_seen)
        inside = historical[name] & rank_set
        k7_seen.update(inside)
        expected_layers.append(layer(
            name, "K7", source_names[name], historical[name], rank_set, prior,
            reported_counts[name], k7_seen, False,
        ))
        layer_sets[f"K7:{name}"] = inside
    if [index for index in rank_indices if index not in k7_seen] != k7_base:
        raise ValueError("historical sets fail to reconstruct K7 base")

    dual, dual_raw_hash, statuses = parse_k7_dual(k7_base)
    if len(dual) != 102 or statuses != {"REJECTED": 102, "SURVIVOR": 12_839}:
        raise ValueError("K7 dual archive counts mismatch")
    prior = set(k7_seen)
    k7_seen.update(dual)
    expected_layers.append(layer(
        "positive_dual_full", "K7", "d6_k7_positive_dual_full_decisions.tsv.gz",
        dual, rank_set, prior, 102, k7_seen, True,
    ))
    layer_sets["K7:positive_dual_full"] = dual
    reflection = as_set(load("d6_k7_reflection_overlap_full_report.json")["rejections"], "reflection")
    prior = set(k7_seen)
    k7_seen.update(reflection)
    expected_layers.append(layer(
        "reflection_overlap_full", "K7", "d6_k7_reflection_overlap_full_report.json",
        reflection, rank_set, prior, 0, k7_seen, True,
    ))
    layer_sets["K7:reflection_overlap_full"] = reflection
    k7_final = [index for index in k7_base if index not in dual and index not in reflection]

    k6_input = load("d6_k6_bipartite_rank_input.json")
    k6_indices = [int(row["index"]) for row in k6_input["graphs"]]
    k6_set = set(k6_indices)
    prior_report = load("d6_k6_bipartite_rank_report.json")
    if [int(row["index"]) for row in prior_report["graph_results"]] != k6_indices:
        raise ValueError("K6 prior row order mismatch")
    bip = {
        int(row["index"])
        for row in prior_report["graph_results"]
        if row["decision"]["rejected"]
    }
    k6_base = [index for index in k6_indices if index not in bip]
    normal, normal_raw_hash = parse_k6_normal(k6_base)
    if bip != {461_363} or len(normal) != 107:
        raise ValueError("K6 rejection counts mismatch")
    k6_seen = set(bip)
    expected_layers.append(layer(
        "bipartite_rank", "K6-only", "d6_k6_bipartite_rank_report.json",
        bip, k6_set, set(), 1, k6_seen, True,
    ))
    layer_sets["K6-only:bipartite_rank"] = bip
    prior = set(k6_seen)
    k6_seen.update(normal)
    expected_layers.append(layer(
        "normal_inertia_full", "K6-only", "d6_k6_normal_inertia_full_decisions.jsonl.gz",
        normal, k6_set, prior, 107, k6_seen, True,
    ))
    layer_sets["K6-only:normal_inertia_full"] = normal
    k6_final = [index for index in k6_base if index not in normal]
    combined = sorted(set(k7_final) | set(k6_final))
    if (len(k7_final), len(k6_final), len(combined)) != (12_839, 990, 13_829):
        raise ValueError("reconstructed final counts mismatch")

    if manifest["layers"] != expected_layers:
        raise ValueError("manifest per-layer accounting differs from independent reconstruction")
    if manifest["pairwise_layer_overlaps"] != overlap_rows(layer_sets):
        raise ValueError("manifest pairwise overlap table differs")
    k7_class = manifest["classes"]["K7"]
    k6_class = manifest["classes"]["K6_only"]
    if k7_class["final_indices"] != k7_final or k7_class["final_indices_sha256"] != stable_hash(k7_final):
        raise ValueError("manifest K7 final indices differ")
    if k6_class["final_indices"] != k6_final or k6_class["final_indices_sha256"] != stable_hash(k6_final):
        raise ValueError("manifest K6 final indices differ")
    if manifest["combined"]["indices"] != combined or manifest["combined"]["indices_sha256"] != stable_hash(combined):
        raise ValueError("manifest combined indices differ")
    if k7_class["archive_checks"]["decompressed_sha256"] != dual_raw_hash:
        raise ValueError("manifest K7 archive hash differs")
    if k6_class["archive_checks"]["decompressed_sha256"] != normal_raw_hash:
        raise ValueError("manifest K6 archive hash differs")
    source_table = {row["path"]: row["sha256"] for row in manifest["source_artifacts"]}
    expected_source_table = {
        name: digest
        for name, digest in EXPECTED.items()
        if name not in {"d6_current_residue_manifest.json"}
    }
    if source_table != expected_source_table:
        raise ValueError("manifest source-artifact table differs")
    return {
        "schema": "d6-current-exact-residue-verification-v1",
        "status": "PASS",
        "method": "independent source-list and decision-archive reconstruction",
        "K7": len(k7_final),
        "K6_only": len(k6_final),
        "combined": len(combined),
        "K7_final_indices_sha256": stable_hash(k7_final),
        "K6_only_final_indices_sha256": stable_hash(k6_final),
        "combined_indices_sha256": stable_hash(combined),
        "layers_checked": len(expected_layers),
        "pairwise_overlaps_checked": len(overlap_rows(layer_sets)),
        "inputs": observed,
        "verifier_sha256": sha256(Path(__file__).resolve()),
    }


def atomic_json(output: Path, value: object) -> None:
    temporary = output.with_name(f".{output.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=path("d6_current_residue_manifest.json"))
    parser.add_argument(
        "--output", type=Path,
        default=path("d6_current_residue_manifest_verification.json"),
    )
    args = parser.parse_args()
    result = verify(args.manifest.resolve())
    atomic_json(args.output.resolve(), result)
    print(json.dumps({
        "status": result["status"], "K7": result["K7"],
        "K6_only": result["K6_only"], "combined": result["combined"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
