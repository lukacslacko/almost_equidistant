#!/usr/bin/env python3
"""Build the hash-pinned current exact dimension-six residue manifest.

The manifest independently reconstructs two disjoint classes:

* the 12,941-graph frozen K7 residue, minus 102 verified full positive-dual
  rejections, leaving 12,839;
* the 1,097-graph frozen K6-only residue, minus 107 verified normal-inertia
  rejections, leaving 990.

Every final corpus index is written explicitly.  Production decision archives
are parsed rather than merely trusting report totals, and every source is
hash-pinned.  Historical K7 layers are replayed through the already frozen
selection builder and exposed individually in the new accounting.
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
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parent
EXPECTED = {
    ".runs/d6_k7_rank_survivors.json": (
        "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
    ),
    "build_d6_k7_positive_dual_selection.py": (
        "5717d24d448c1dc8b1da7dbb8dca53902d9b40546e9f4edc70c3183fde995dad"
    ),
    "d6_k7_positive_dual_selection.json": (
        "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0"
    ),
    "d6_k7_support_rank_survivors_report.json": (
        "5138d207543967b393baae9b14b47b49a276273cbc4468ee843856891917691e"
    ),
    "d6_k7_strict_h_full_report.json": (
        "59be9c6a2d4cd2e8e8214d28d3198487846e58697046415ae1e15f4f69dc3add"
    ),
    "d6_k7_small_support_value_full_report.json": (
        "cea64f3dde804e0766c5b77e373aa50338d5bee6e5e54f44749f4e387cf529aa"
    ),
    "d6_interval_residue_benchmark.json": (
        "e246c3896e7ff2b9b40a598efdaffcc465ff005700a0d51fb17b66b2d6f68481"
    ),
    "d6_k7_positive_polynomial_dual_report.json": (
        "51f29654f9514eec8f3911149a42407a7f197470720061528064ab77c9dec120"
    ),
    "d6_k7_positive_dual_full_report.json": (
        "c2de7e06bbe4f3d163a4f747d42721f7f60da985b867ce74669664431738e345"
    ),
    "d6_k7_positive_dual_full_verification.json": (
        "90df584f7a9b15438c3c0dd7ca04f49a586df6e73031da1cc55f8325b5ef950e"
    ),
    "d6_k7_positive_dual_full_decisions.tsv.gz": (
        "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
    ),
    "d6_k7_positive_dual_full_certificates.jsonl.gz": (
        "3adaa7e562dfb9241f205e3e7d2384805dbe7296f9677d76f2f23ce913a4ce7f"
    ),
    "d6_k7_positive_dual_full_checkpoint.json": (
        "971d1d38db145a49a9dc5cec61f9c2bf9ba17eccee9a26a90f59af49f6777ec4"
    ),
    "d6_k7_reflection_overlap_full_report.json": (
        "d8294f1d1dd11206856ace4bfe210cbc79f4bb019969930645d95a7bad89cd07"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "d6_k6_bipartite_rank_report.json": (
        "ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc"
    ),
    "d6_k6_normal_inertia_full_report.json": (
        "8fd3d9d10038040f42407b25ea007f59f0e84e6eb0da8bec2694d1831e1991c5"
    ),
    "d6_k6_normal_inertia_full_verification_report.json": (
        "61889185deb0e773bf899cdbb82c96467cc38a4688560a50f53c7761a7615e6b"
    ),
    "d6_k6_normal_inertia_full_decisions.jsonl.gz": (
        "614bbb64989b5834ed0dbe1c12e7718583e53cfd68483668ffbbad44608b3828"
    ),
    "d6_k6_normal_inertia_full_checkpoint.json": (
        "f4cb37faab90ed65b08b0e0ae2015a226de577caad5c74d23b530377d089962d"
    ),
}
EXPECTED_K7_RANK = 17_764
EXPECTED_K7_PRIOR = 12_941
EXPECTED_K7_NEW_REJECTIONS = 102
EXPECTED_K7_FINAL = 12_839
EXPECTED_K6_INPUT = 1_098
EXPECTED_K6_PRIOR = 1_097
EXPECTED_K6_NEW_REJECTIONS = 107
EXPECTED_K6_FINAL = 990
EXPECTED_COMBINED = 13_829


def source_path(name: str) -> Path:
    return ROOT / name


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


def verify_sources() -> dict[str, str]:
    observed = {name: sha256(source_path(name)) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"source hash mismatch: observed={observed}, expected={EXPECTED}")
    return observed


def read_json(name: str) -> dict:
    value = json.loads(source_path(name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{name} is not a JSON object")
    return value


def integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or not all(
        isinstance(item, int) and not isinstance(item, bool) for item in value
    ):
        raise ValueError(f"{label} is not an integer list")
    if len(set(value)) != len(value):
        raise ValueError(f"{label} contains duplicates")
    return value


def decompressed_bytes(name: str) -> tuple[bytes, str]:
    with gzip.open(source_path(name), "rb") as stream:
        raw = stream.read()
    return raw, hashlib.sha256(raw).hexdigest()


def k7_dual_decisions(base_indices: Sequence[int]) -> tuple[set[int], dict]:
    raw, raw_hash = decompressed_bytes("d6_k7_positive_dual_full_decisions.tsv.gz")
    report = read_json("d6_k7_positive_dual_full_report.json")
    if raw_hash != report.get("artifacts", {}).get("decisions_uncompressed_sha256"):
        raise ValueError("K7 decompressed decision hash does not match report")
    text = io.StringIO(raw.decode("ascii"), newline="")
    reader = csv.DictReader(text, delimiter="\t")
    expected_fields = [
        "ordinal", "index", "status", "applicable", "seeds", "covers",
        "enhanced_passing_covers", "strict_h_passing_covers",
        "dual_passing_covers", "strict_h_passing_cliques",
        "dual_failing_cliques", "strict_h_rejected", "dual_rejected",
        "marginal_dual_rejected", "error_type",
    ]
    if reader.fieldnames != expected_fields:
        raise ValueError(f"unexpected K7 decision schema: {reader.fieldnames}")
    rejected: set[int] = set()
    rows = 0
    status_counts: dict[str, int] = {}
    for ordinal, (row, expected_index) in enumerate(zip(reader, base_indices)):
        rows += 1
        if int(row["ordinal"]) != ordinal or int(row["index"]) != expected_index:
            raise ValueError(f"K7 decision prefix mismatch at ordinal {ordinal}")
        status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
        if row["error_type"] or row["status"] == "INFRA_ERROR":
            raise ValueError(f"K7 infrastructure error at index {expected_index}")
        marginal = int(row["marginal_dual_rejected"])
        if marginal not in (0, 1):
            raise ValueError("K7 marginal decision is not Boolean")
        if marginal:
            if row["status"] != "REJECTED" or int(row["dual_rejected"]) != 1:
                raise ValueError("K7 marginal rejection/status inconsistency")
            rejected.add(expected_index)
    if next(reader, None) is not None or rows != len(base_indices):
        raise ValueError("K7 decision archive row count mismatch")
    return rejected, {
        "rows": rows,
        "status_counts": status_counts,
        "decompressed_sha256": raw_hash,
    }


def k6_normal_decisions(base_indices: Sequence[int]) -> tuple[set[int], dict]:
    raw, raw_hash = decompressed_bytes("d6_k6_normal_inertia_full_decisions.jsonl.gz")
    report = read_json("d6_k6_normal_inertia_full_report.json")
    if raw_hash != report.get("artifacts", {}).get("decisions_sha256"):
        raise ValueError("K6 decompressed decision hash does not match report")
    rejected: set[int] = set()
    rows = []
    for line_number, line in enumerate(raw.splitlines(), 1):
        try:
            row = json.loads(line)
        except json.JSONDecodeError as error:
            raise ValueError(f"invalid K6 decision JSONL line {line_number}") from error
        rows.append(row)
    if len(rows) != len(base_indices):
        raise ValueError("K6 decision archive row count mismatch")
    for position, (row, expected_index) in enumerate(zip(rows, base_indices)):
        if (
            row.get("schema") != "d6-k6-normal-inertia-full-decision-v1"
            or row.get("index") != expected_index
        ):
            raise ValueError(f"K6 decision prefix mismatch at row {position}")
        decision = row.get("decision", {})
        if decision.get("applicable") is not True:
            raise ValueError(f"K6 inapplicable row at index {expected_index}")
        if decision.get("rejected"):
            rejected.add(expected_index)
    return rejected, {"rows": len(rows), "decompressed_sha256": raw_hash}


def layer_record(
    layer: str,
    population: str,
    source: str,
    raw: set[int],
    universe: set[int],
    prior: set[int],
    reported_rejections: int,
    include_indices: bool = False,
) -> dict:
    in_universe = raw & universe
    incremental = in_universe - prior
    overlap = in_universe & prior
    record = {
        "layer": layer,
        "population": population,
        "source": source,
        "reported_rejections": reported_rejections,
        "raw_rejections": len(raw),
        "raw_rejections_sha256": stable_hash(sorted(raw)),
        "rejections_in_universe": len(in_universe),
        "rejections_in_universe_sha256": stable_hash(sorted(in_universe)),
        "rejections_outside_universe": len(raw - universe),
        "overlap_with_prior_layers": len(overlap),
        "overlap_with_prior_layers_sha256": stable_hash(sorted(overlap)),
        "incremental_rejections": len(incremental),
        "incremental_rejections_sha256": stable_hash(sorted(incremental)),
    }
    if include_indices:
        record["incremental_rejected_indices"] = sorted(incremental)
    return record


def overlap_table(layer_sets: dict[str, set[int]]) -> list[dict]:
    answer = []
    for left, right in combinations(layer_sets, 2):
        overlap = layer_sets[left] & layer_sets[right]
        answer.append({
            "left": left,
            "right": right,
            "count": len(overlap),
            "indices_sha256": stable_hash(sorted(overlap)),
        })
    return answer


def build_manifest() -> dict:
    observed_sources = verify_sources()

    # Re-run the frozen K7 five-layer builder and demand byte-level semantic
    # equality with its committed explicit 12,941-index selection.
    import build_d6_k7_positive_dual_selection as frozen_k7

    rank_path = source_path(".runs/d6_k7_rank_survivors.json")
    rebuilt_selection = frozen_k7.build_selection(rank_path)
    frozen_selection = read_json("d6_k7_positive_dual_selection.json")
    if rebuilt_selection != frozen_selection:
        raise ValueError("frozen K7 selection does not match an independent rebuild")
    k7_payload = json.loads(rank_path.read_text(encoding="utf-8"))
    k7_universe_list = [int(record["index"]) for record in k7_payload["graphs"]]
    k7_universe = set(k7_universe_list)
    if len(k7_universe) != EXPECTED_K7_RANK:
        raise ValueError("K7 rank-survivor universe count mismatch")
    k7_base = integer_list(frozen_selection.get("selected_indices"), "K7 base")
    if len(k7_base) != EXPECTED_K7_PRIOR or not set(k7_base) <= k7_universe:
        raise ValueError("K7 frozen residue mismatch")

    historical_sources = frozen_k7.source_rejection_sets(
        frozen_selection["input_sha256"]
    )
    k7_layers = []
    k7_layer_sets: dict[str, set[int]] = {}
    prior_k7: set[int] = set()
    frozen_accounting = {
        item["layer"]: item for item in frozen_selection["layer_accounting"]
    }
    historical_report_names = {
        "support_propagation": "d6_k7_support_rank_survivors_report.json",
        "strict_H": "d6_k7_strict_h_full_report.json",
        "sparse_value": "d6_k7_small_support_value_full_report.json",
        "interval_cap20000": "d6_interval_residue_benchmark.json",
        "sample_positive_dual": "d6_k7_positive_polynomial_dual_report.json",
    }
    for name, raw, metadata in historical_sources:
        in_universe = raw & k7_universe
        record = layer_record(
            name,
            "K7",
            historical_report_names[name],
            raw,
            k7_universe,
            prior_k7,
            int(metadata["reported_rejections"]),
        )
        expected_accounting = frozen_accounting[name]
        for key in (
            "rejections_in_rank_survivors",
            "rejections_outside_rank_survivors",
            "overlap_with_prior_layers",
            "incremental_rejections",
            "rejections_in_rank_survivors_sha256",
            "incremental_rejections_sha256",
        ):
            translated = {
                "rejections_in_rank_survivors": "rejections_in_universe",
                "rejections_outside_rank_survivors": "rejections_outside_universe",
                "rejections_in_rank_survivors_sha256": "rejections_in_universe_sha256",
            }.get(key, key)
            if record[translated] != expected_accounting[key]:
                raise ValueError(f"historical K7 accounting mismatch at {name}:{key}")
        prior_k7.update(in_universe)
        record["cumulative_rejections"] = len(prior_k7)
        record["cumulative_survivors"] = len(k7_universe - prior_k7)
        k7_layers.append(record)
        k7_layer_sets[f"K7:{name}"] = in_universe
    if [index for index in k7_universe_list if index not in prior_k7] != k7_base:
        raise ValueError("historical K7 layer union does not reproduce frozen base")

    dual_report = read_json("d6_k7_positive_dual_full_report.json")
    dual_verification = read_json("d6_k7_positive_dual_full_verification.json")
    if (
        dual_report.get("summary", {}).get("graphs") != EXPECTED_K7_PRIOR
        or dual_report.get("summary", {}).get("marginal_dual_rejected")
        != EXPECTED_K7_NEW_REJECTIONS
        or dual_report.get("summary", {}).get("infra_errors") != 0
        or dual_verification.get("status") != "PASS"
        or dual_verification.get("result", {}).get("rejected_graph_quantifiers_rebuilt")
        != EXPECTED_K7_NEW_REJECTIONS
    ):
        raise ValueError("K7 full dual report/verification mismatch")
    dual_rejected, dual_archive = k7_dual_decisions(k7_base)
    report_dual_set = set(integer_list(
        dual_report["summary"]["marginal_dual_rejected_indices"],
        "K7 full dual report rejections",
    ))
    if dual_rejected != report_dual_set or len(dual_rejected) != EXPECTED_K7_NEW_REJECTIONS:
        raise ValueError("K7 dual archive/report rejection set mismatch")
    k7_seen = set(prior_k7)
    dual_layer = layer_record(
        "positive_dual_full",
        "K7",
        "d6_k7_positive_dual_full_decisions.tsv.gz",
        dual_rejected,
        k7_universe,
        k7_seen,
        EXPECTED_K7_NEW_REJECTIONS,
        include_indices=True,
    )
    k7_seen.update(dual_rejected)
    dual_layer["cumulative_rejections"] = len(k7_seen)
    dual_layer["cumulative_survivors"] = len(k7_universe - k7_seen)
    k7_layers.append(dual_layer)
    k7_layer_sets["K7:positive_dual_full"] = dual_rejected

    reflection = read_json("d6_k7_reflection_overlap_full_report.json")
    reflection_set = set(integer_list(reflection.get("rejections"), "reflection rejections"))
    if (
        reflection.get("counts", {}).get("graphs") != EXPECTED_K7_PRIOR
        or reflection.get("counts", {}).get("rejected") != 0
        or reflection.get("scope", {}).get("sample_indices_sha256")
        != frozen_selection["selected_indices_sha256"]
        or reflection_set
    ):
        raise ValueError("K7 reflection full report mismatch")
    reflection_layer = layer_record(
        "reflection_overlap_full",
        "K7",
        "d6_k7_reflection_overlap_full_report.json",
        reflection_set,
        k7_universe,
        k7_seen,
        0,
        include_indices=True,
    )
    reflection_layer["cumulative_rejections"] = len(k7_seen)
    reflection_layer["cumulative_survivors"] = len(k7_universe - k7_seen)
    k7_layers.append(reflection_layer)
    k7_layer_sets["K7:reflection_overlap_full"] = reflection_set

    k7_final = [index for index in k7_base if index not in dual_rejected]
    if len(k7_final) != EXPECTED_K7_FINAL:
        raise ValueError("K7 final residue count mismatch")

    # Reconstruct the frozen 1,097 K6-only base directly from its 1,098-row
    # input/report, then derive the new rejection set from the JSONL archive.
    k6_input = read_json("d6_k6_bipartite_rank_input.json")
    k6_prior_report = read_json("d6_k6_bipartite_rank_report.json")
    k6_universe_list = [int(record["index"]) for record in k6_input["graphs"]]
    k6_universe = set(k6_universe_list)
    if len(k6_universe) != EXPECTED_K6_INPUT:
        raise ValueError("K6 input universe count mismatch")
    prior_rows = k6_prior_report.get("graph_results", [])
    if [int(row["index"]) for row in prior_rows] != k6_universe_list:
        raise ValueError("K6 prior report row order mismatch")
    bip_rejected = {
        int(row["index"]) for row in prior_rows if row["decision"]["rejected"]
    }
    if bip_rejected != {461_363} or set(k6_prior_report["rejected_indices"]) != bip_rejected:
        raise ValueError("K6 prior rejection set mismatch")
    k6_base = [index for index in k6_universe_list if index not in bip_rejected]
    if len(k6_base) != EXPECTED_K6_PRIOR:
        raise ValueError("K6 frozen base count mismatch")

    normal_report = read_json("d6_k6_normal_inertia_full_report.json")
    normal_verification = read_json("d6_k6_normal_inertia_full_verification_report.json")
    if (
        normal_report.get("status") != "COMPLETE"
        or normal_report.get("input_graphs") != EXPECTED_K6_PRIOR
        or normal_report.get("graphs_rejected") != EXPECTED_K6_NEW_REJECTIONS
        or normal_verification.get("status") != "PASS"
        or normal_verification.get("graphs_recomputed") != EXPECTED_K6_PRIOR
        or normal_verification.get("graphs_rejected") != EXPECTED_K6_NEW_REJECTIONS
    ):
        raise ValueError("K6 normal-inertia report/verification mismatch")
    normal_rejected, normal_archive = k6_normal_decisions(k6_base)
    if (
        normal_rejected != set(integer_list(normal_report["rejected_indices"], "K6 report rejections"))
        or normal_rejected
        != set(integer_list(normal_verification["rejected_indices"], "K6 verifier rejections"))
        or len(normal_rejected) != EXPECTED_K6_NEW_REJECTIONS
    ):
        raise ValueError("K6 normal archive/report rejection set mismatch")

    k6_layers = []
    k6_layer_sets = {}
    bip_layer = layer_record(
        "bipartite_rank",
        "K6-only",
        "d6_k6_bipartite_rank_report.json",
        bip_rejected,
        k6_universe,
        set(),
        1,
        include_indices=True,
    )
    bip_layer["cumulative_rejections"] = 1
    bip_layer["cumulative_survivors"] = EXPECTED_K6_PRIOR
    k6_layers.append(bip_layer)
    k6_layer_sets["K6-only:bipartite_rank"] = bip_rejected
    normal_layer = layer_record(
        "normal_inertia_full",
        "K6-only",
        "d6_k6_normal_inertia_full_decisions.jsonl.gz",
        normal_rejected,
        k6_universe,
        bip_rejected,
        EXPECTED_K6_NEW_REJECTIONS,
        include_indices=True,
    )
    normal_layer["cumulative_rejections"] = 1 + len(normal_rejected)
    normal_layer["cumulative_survivors"] = EXPECTED_K6_FINAL
    k6_layers.append(normal_layer)
    k6_layer_sets["K6-only:normal_inertia_full"] = normal_rejected
    k6_final = [index for index in k6_base if index not in normal_rejected]
    if len(k6_final) != EXPECTED_K6_FINAL:
        raise ValueError("K6 final residue count mismatch")

    if set(k7_universe) & set(k6_universe):
        raise ValueError("K7 and K6-only source populations overlap")
    combined = sorted(set(k7_final) | set(k6_final))
    if len(combined) != EXPECTED_COMBINED:
        raise ValueError("combined residue count mismatch")
    all_layer_sets = {**k7_layer_sets, **k6_layer_sets}

    sources = [
        {"path": name, "sha256": digest}
        for name, digest in observed_sources.items()
    ]
    sources.append({
        "path": Path(__file__).name,
        "sha256": sha256(Path(__file__).resolve()),
    })
    return {
        "schema": "d6-current-exact-residue-v1",
        "status": "COMPLETE",
        "description": (
            "Exact unresolved dimension-six level-19 residue after the frozen "
            "K7 layers including full positive duals and the frozen K6-only "
            "layers including full normal inertia. Survivors are nonclaims."
        ),
        "classes": {
            "K7": {
                "source_universe": EXPECTED_K7_RANK,
                "source_indices_sha256": stable_hash(k7_universe_list),
                "prior_exact_residue": EXPECTED_K7_PRIOR,
                "prior_exact_residue_indices_sha256": stable_hash(k7_base),
                "new_full_dual_rejections": len(dual_rejected),
                "new_full_dual_rejected_indices_sha256": stable_hash(sorted(dual_rejected)),
                "final_count": len(k7_final),
                "final_indices": k7_final,
                "final_indices_sha256": stable_hash(k7_final),
                "final_sorted_indices_sha256": stable_hash(sorted(k7_final)),
                "adjacency_source": ".runs/d6_k7_rank_survivors.json",
                "adjacency_source_sha256": EXPECTED[".runs/d6_k7_rank_survivors.json"],
                "archive_checks": dual_archive,
            },
            "K6_only": {
                "source_universe": EXPECTED_K6_INPUT,
                "source_indices_sha256": stable_hash(k6_universe_list),
                "prior_exact_residue": EXPECTED_K6_PRIOR,
                "prior_exact_residue_indices_sha256": stable_hash(k6_base),
                "new_normal_inertia_rejections": len(normal_rejected),
                "new_normal_inertia_rejected_indices_sha256": stable_hash(sorted(normal_rejected)),
                "final_count": len(k6_final),
                "final_indices": k6_final,
                "final_indices_sha256": stable_hash(k6_final),
                "final_sorted_indices_sha256": stable_hash(sorted(k6_final)),
                "adjacency_source": "d6_k6_bipartite_rank_input.json",
                "adjacency_source_sha256": EXPECTED["d6_k6_bipartite_rank_input.json"],
                "archive_checks": normal_archive,
            },
        },
        "layers": k7_layers + k6_layers,
        "pairwise_layer_overlaps": overlap_table(all_layer_sets),
        "cross_class": {
            "source_universe_overlap": 0,
            "source_universe_overlap_indices_sha256": stable_hash([]),
            "final_overlap": 0,
            "final_overlap_indices_sha256": stable_hash([]),
        },
        "combined": {
            "count": len(combined),
            "indices": combined,
            "indices_sha256": stable_hash(combined),
            "K7_count": len(k7_final),
            "K6_only_count": len(k6_final),
        },
        "source_artifacts": sources,
        "semantics": {
            "survivors": "unresolved; no realizability claim",
            "candidate_nonedges": "unconstrained and may also have unit distance",
            "allowed_defect_coordinates": "may be zero",
        },
    }


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest.json",
    )
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output.resolve(), manifest)
    print(json.dumps({
        "status": manifest["status"],
        "K7": manifest["classes"]["K7"]["final_count"],
        "K6_only": manifest["classes"]["K6_only"]["final_count"],
        "combined": manifest["combined"]["count"],
        "output": str(args.output.resolve()),
        "sha256": sha256(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
