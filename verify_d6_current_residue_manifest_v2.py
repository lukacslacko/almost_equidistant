#!/usr/bin/env python3
"""Independent checker for ``d6_current_residue_manifest_v2.json``.

This checker imports neither the v2 builder nor any K7/K6 production engine.
It separately reconstructs both ordered complements from the frozen upstream
JSON artifacts, checks both independent-verification gates, validates every
embedded adjacency bitmask, and compares all manifest accounting fields.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
BOUNDARY_COMMIT = "fc51874458196c0391530d0e98968f530ab738c9"
BUILDER = "build_d6_current_residue_manifest_v2.py"
EXPECTED_BUILDER_SHA256 = (
    "f3db5ce7bbca61f432a73ca56184528504184ddce2d2c2cacf9e70db4d7ad5f1"
)
EXPECTED_SOURCES = {
    ".runs/d6_k7_rank_survivors.json": (
        "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
    ),
    "d6_current_residue_manifest.json": (
        "d06cdb23b03b1f8523fbe3c3ede5bc297ace21a844282830c5cf51a6cee3512d"
    ),
    "d6_current_residue_manifest_verification.json": (
        "15f4797bf7a3095dbf0a90ec234b07bdaa5a9ec35599a0913084dd5d17d11ae7"
    ),
    "d6_k7_rankone_tetrad_full_selection.json": (
        "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
    ),
    "d6_k7_rankone_pattern_union.json": (
        "1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599"
    ),
    "d6_k7_rankone_pattern_union_verification.json": (
        "7b9d5e35bf787e5d05aa486f914416598c31836264a479fbbc1af117f1dfcc14"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "d6_k6_fused_side_rank_report.json": (
        "60911844102dfeb494170f6b7aa5d5142b8f545fcc968e43078eb07c7ed102ec"
    ),
    "d6_k6_fused_side_rank_verification.json": (
        "0572e8088d75a4e4fe575fda185f8be0bbb46290c22570f760b92e5d954cdd92"
    ),
}


def path(name: str) -> Path:
    return ROOT / name


def sha256(file_path: Path) -> str:
    digest = hashlib.sha256()
    with file_path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load(file_path: Path) -> dict:
    value = json.loads(file_path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{file_path.name} is not a JSON object")
    return value


def load_named(name: str) -> dict:
    return load(path(name))


def as_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    if len(value) != len(set(value)):
        raise ValueError(f"{label} has duplicates")
    return list(value)


def set_payload(union: dict, name: str) -> tuple[list[int], set[int]]:
    payload = union.get("sets", {}).get(name)
    if not isinstance(payload, dict):
        raise ValueError(f"upstream K7 union has no {name}")
    sequence = as_list(payload.get("indices"), f"upstream K7 {name}")
    if payload.get("count") != len(sequence) or payload.get("indices_sha256") != stable_hash(sequence):
        raise ValueError(f"upstream K7 {name} count/hash mismatch")
    return sequence, set(sequence)


def graph_lookup(records: object, label: str) -> dict[int, list[int]]:
    if not isinstance(records, list):
        raise ValueError(f"{label} graph collection is not a list")
    result: dict[int, list[int]] = {}
    for position, record in enumerate(records):
        if not isinstance(record, dict) or type(record.get("index")) is not int:
            raise ValueError(f"{label} record {position} lacks an integer index")
        index = record["index"]
        adjacency = record.get("adjacency")
        if (
            index in result
            or not isinstance(adjacency, list)
            or len(adjacency) != 19
            or any(type(mask) is not int for mask in adjacency)
        ):
            raise ValueError(f"invalid {label} adjacency record {position}")
        for vertex, mask in enumerate(adjacency):
            if mask < 0 or mask >= 1 << 19 or mask & (1 << vertex):
                raise ValueError(f"invalid {label} row {vertex} at index {index}")
        for left in range(19):
            for right in range(left + 1, 19):
                if ((adjacency[left] >> right) & 1) != ((adjacency[right] >> left) & 1):
                    raise ValueError(f"asymmetric {label} adjacency at index {index}")
        result[index] = list(adjacency)
    return result


def expected_graphs(indices: Sequence[int], lookup: dict[int, list[int]], label: str) -> list[dict]:
    try:
        return [{"index": index, "adjacency": lookup[index]} for index in indices]
    except KeyError as error:
        raise ValueError(f"{label} source misses index {error.args[0]}") from error


def require_equal(observed: object, expected: object, label: str) -> None:
    if observed != expected:
        raise ValueError(f"{label} mismatch")


def verify(manifest_path: Path) -> dict:
    observed_sources = {
        name: sha256(path(name)) for name in EXPECTED_SOURCES
    }
    if observed_sources != EXPECTED_SOURCES:
        raise ValueError("one or more independently pinned source hashes changed")
    builder_hash = sha256(path(BUILDER))
    if builder_hash != EXPECTED_BUILDER_SHA256:
        raise ValueError("v2 builder source hash changed")
    manifest = load(manifest_path)
    require_equal(manifest.get("schema"), "d6-current-exact-residue-v2", "schema")
    require_equal(manifest.get("status"), "COMPLETE", "status")
    require_equal(manifest.get("boundary_commit"), BOUNDARY_COMMIT, "boundary commit")
    require_equal(
        manifest.get("description"),
        (
            "Self-contained exact unresolved dimension-six level-19 boundary "
            "after the independently verified K7 tetrad/pattern union and "
            "K6-only fused side-rank layers."
        ),
        "description",
    )
    require_equal(manifest.get("class_order"), ["K7", "K6_only"], "class order")

    # Rebuild the K7 union partition using only frozen JSON values.
    selection = load_named("d6_k7_rankone_tetrad_full_selection.json")
    selected = as_list(selection.get("selected_indices"), "K7 selection")
    if len(selected) != 12_839 or selection.get("selected_indices_sha256") != stable_hash(selected):
        raise ValueError("K7 selection gate failed")
    union = load_named("d6_k7_rankone_pattern_union.json")
    union_check = load_named("d6_k7_rankone_pattern_union_verification.json")
    if (
        union.get("schema") != "d6-k7-rankone-pattern-union-v1"
        or union.get("status") != "COMPLETE"
        or union_check.get("schema")
        != "d6-k7-rankone-pattern-union-verification-v1"
        or union_check.get("status") != "PASS"
        or union_check.get("manifest", {}).get("sha256")
        != observed_sources["d6_k7_rankone_pattern_union.json"]
    ):
        raise ValueError("K7 independent checker gate failed")
    union_indices, union_set = set_payload(union, "exact_union")
    k7_indices, k7_set = set_payload(union, "exact_residue")
    tetrad_indices, tetrad_set = set_payload(union, "tetrad_rejected")
    pattern_indices, pattern_set = set_payload(union, "pattern_954_rejected")
    both_indices, both_set = set_payload(union, "both")
    tetrad_only_indices, tetrad_only_set = set_payload(union, "tetrad_only")
    pattern_only_indices, pattern_only_set = set_payload(union, "pattern_954_only")
    if (
        len(union_indices) != 12_581
        or len(k7_indices) != 258
        or union_set != tetrad_set | pattern_set
        or both_set != tetrad_set & pattern_set
        or tetrad_only_set != tetrad_set - pattern_set
        or pattern_only_set != pattern_set - tetrad_set
        or union_set & k7_set
        or union_set | k7_set != set(selected)
        or k7_indices != [index for index in selected if index not in union_set]
        or union_check.get("graphs") != 12_839
        or union_check.get("exact_union") != 12_581
        or union_check.get("exact_residue") != 258
        or union_check.get("residue_indices_sha256") != stable_hash(k7_indices)
    ):
        raise ValueError("K7 partition reconstruction failed")
    k7_source = load_named(".runs/d6_k7_rank_survivors.json")
    k7_lookup = graph_lookup(k7_source.get("graphs"), "K7 source")
    k7_graphs = expected_graphs(k7_indices, k7_lookup, "K7")

    # Rebuild the K6-only complement and independently check its positive control.
    old = load_named("d6_current_residue_manifest.json")
    old_check = load_named("d6_current_residue_manifest_verification.json")
    if (
        old.get("schema") != "d6-current-exact-residue-v1"
        or old.get("status") != "COMPLETE"
        or old_check.get("schema") != "d6-current-exact-residue-verification-v1"
        or old_check.get("status") != "PASS"
        or old_check.get("combined") != 13_829
    ):
        raise ValueError("prior-boundary independent checker gate failed")
    k6_base = as_list(
        old.get("classes", {}).get("K6_only", {}).get("final_indices"),
        "K6 prior boundary",
    )
    if len(k6_base) != 990 or stable_hash(k6_base) != (
        "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
    ):
        raise ValueError("K6 prior boundary failed")
    k6_report = load_named("d6_k6_fused_side_rank_report.json")
    k6_check = load_named("d6_k6_fused_side_rank_verification.json")
    k6_rejected = as_list(k6_report.get("rejected_indices"), "K6 report rejections")
    checked_rejected = as_list(
        k6_check.get("rejected_indices"), "K6 verifier rejections"
    )
    report_order = [row.get("index") for row in k6_report.get("graph_results", [])]
    if (
        k6_report.get("schema") != "d6-k6-fused-side-rank-v1"
        or k6_report.get("status") != "COMPLETE"
        or k6_check.get("schema") != "d6-k6-fused-side-rank-verification-v1"
        or k6_check.get("status") != "PASS"
        or k6_check.get("inputs", {}).get("d6_k6_fused_side_rank_report.json")
        != observed_sources["d6_k6_fused_side_rank_report.json"]
        or report_order != k6_base
        or k6_report.get("input_graphs") != 990
        or k6_report.get("input_indices_sha256") != stable_hash(k6_base)
        or k6_report.get("graphs_rejected") != 13
        or k6_report.get("graphs_surviving") != 977
        or k6_check.get("graphs_recomputed") != 990
        or k6_check.get("graphs_rejected") != 13
        or k6_check.get("graphs_surviving") != 977
        or k6_rejected != checked_rejected
        or k6_report.get("positive_control", {}).get("passed") is not True
        or k6_check.get("positive_control", {}).get("passed") is not True
    ):
        raise ValueError("K6 independent checker/positive-control gate failed")
    k6_rejected_set = set(k6_rejected)
    if len(k6_rejected_set) != 13 or not k6_rejected_set <= set(k6_base):
        raise ValueError("K6 rejection set is not a 13-element subset")
    k6_indices = [index for index in k6_base if index not in k6_rejected_set]
    if len(k6_indices) != 977:
        raise ValueError("K6 ordered complement count mismatch")
    k6_source = load_named("d6_k6_bipartite_rank_input.json")
    k6_lookup = graph_lookup(k6_source.get("graphs"), "K6 source")
    k6_graphs = expected_graphs(k6_indices, k6_lookup, "K6")

    if set(selected) & set(k6_lookup) or k7_set & set(k6_indices):
        raise ValueError("K7 and K6-only classes are not disjoint")

    expected_k7 = {
        "input_boundary_count": 12_839,
        "input_boundary_indices_sha256": stable_hash(selected),
        "exact_rejections": 12_581,
        "exact_rejected_indices_sha256": stable_hash(union_indices),
        "residue_count": 258,
        "residue_indices": k7_indices,
        "residue_indices_sha256": stable_hash(k7_indices),
        "graphs": k7_graphs,
        "graphs_sha256": stable_hash(k7_graphs),
        "adjacency_origin": {
            "path": ".runs/d6_k7_rank_survivors.json",
            "sha256": observed_sources[".runs/d6_k7_rank_survivors.json"],
            "source_graphs": len(k7_lookup),
            "embedded": True,
        },
    }
    expected_k6 = {
        "input_boundary_count": 990,
        "input_boundary_indices_sha256": stable_hash(k6_base),
        "exact_rejections": 13,
        "exact_rejected_indices": k6_rejected,
        "exact_rejected_indices_sha256": stable_hash(k6_rejected),
        "residue_count": 977,
        "residue_indices": k6_indices,
        "residue_indices_sha256": stable_hash(k6_indices),
        "graphs": k6_graphs,
        "graphs_sha256": stable_hash(k6_graphs),
        "adjacency_origin": {
            "path": "d6_k6_bipartite_rank_input.json",
            "sha256": observed_sources["d6_k6_bipartite_rank_input.json"],
            "source_graphs": len(k6_lookup),
            "embedded": True,
        },
    }
    require_equal(manifest.get("classes", {}).get("K7"), expected_k7, "K7 class")
    require_equal(
        manifest.get("classes", {}).get("K6_only"), expected_k6, "K6-only class"
    )
    require_equal(set(manifest.get("classes", {})), {"K7", "K6_only"}, "class keys")

    concatenated = k7_indices + k6_indices
    tagged_indices = (
        [{"class": "K7", "index": index} for index in k7_indices]
        + [{"class": "K6_only", "index": index} for index in k6_indices]
    )
    tagged_graphs = (
        [{"class": "K7", **record} for record in k7_graphs]
        + [{"class": "K6_only", **record} for record in k6_graphs]
    )
    expected_combined = {
        "count": 1_235,
        "class_counts": {"K7": 258, "K6_only": 977},
        "ordered_indices_sha256": stable_hash(concatenated),
        "sorted_indices_sha256": stable_hash(sorted(concatenated)),
        "ordered_class_index_records_sha256": stable_hash(tagged_indices),
        "ordered_class_graph_records_sha256": stable_hash(tagged_graphs),
        "cross_class_overlap": 0,
        "cross_class_overlap_indices_sha256": stable_hash([]),
    }
    require_equal(manifest.get("combined"), expected_combined, "combined accounting")

    expected_gates = {
        "prior_boundary": {
            "manifest": "d6_current_residue_manifest.json",
            "manifest_sha256": observed_sources["d6_current_residue_manifest.json"],
            "verification": "d6_current_residue_manifest_verification.json",
            "verification_sha256": observed_sources[
                "d6_current_residue_manifest_verification.json"
            ],
            "verification_status": "PASS",
        },
        "K7": {
            "manifest": "d6_k7_rankone_pattern_union.json",
            "manifest_sha256": observed_sources["d6_k7_rankone_pattern_union.json"],
            "verification": "d6_k7_rankone_pattern_union_verification.json",
            "verification_sha256": observed_sources[
                "d6_k7_rankone_pattern_union_verification.json"
            ],
            "verification_status": "PASS",
            "exact_union": 12_581,
            "exact_residue": 258,
        },
        "K6_only": {
            "report": "d6_k6_fused_side_rank_report.json",
            "report_sha256": observed_sources["d6_k6_fused_side_rank_report.json"],
            "verification": "d6_k6_fused_side_rank_verification.json",
            "verification_sha256": observed_sources[
                "d6_k6_fused_side_rank_verification.json"
            ],
            "verification_status": "PASS",
            "positive_control_passed": True,
            "exact_rejections": 13,
            "exact_residue": 977,
        },
    }
    require_equal(manifest.get("verification_gates"), expected_gates, "verification gates")
    expected_hashes = dict(observed_sources)
    expected_hashes[BUILDER] = builder_hash
    require_equal(manifest.get("source_hashes"), expected_hashes, "source hashes")

    expected_semantics = {
        "unit_graph_edges": "required to have Euclidean distance exactly 1",
        "candidate_nonedges": "unconstrained and may also have distance 1",
        "points": "must be distinct",
        "allowed_defect_coordinates": "upper bounds only and may be zero",
        "membership_arithmetic": (
            "exact set filtering and integer bitmask copying; no floating-point decision"
        ),
        "residue": "filter non-rejection only; not a realizability claim",
    }
    expected_nonclaims = [
        "No graph in the residue is claimed realizable in R6.",
        "No candidate nonedge is constrained to be non-unit.",
        "The combined count 1235 is not yet a proof that f(6)=18.",
        "The manifest does not replace the upstream exact rejection certificates.",
    ]
    require_equal(manifest.get("semantics"), expected_semantics, "semantics")
    require_equal(manifest.get("nonclaims"), expected_nonclaims, "nonclaims")
    require_equal(
        set(manifest),
        {
            "schema", "status", "boundary_commit", "description", "class_order",
            "classes", "combined", "verification_gates", "source_hashes",
            "semantics", "nonclaims",
        },
        "top-level keys",
    )

    return {
        "schema": "d6-current-exact-residue-v2-verification-v1",
        "status": "PASS",
        "boundary_commit": BOUNDARY_COMMIT,
        "manifest": {
            "path": manifest_path.name,
            "sha256": sha256(manifest_path),
        },
        "classes": {
            "K7": {
                "input": 12_839,
                "exact_rejections": 12_581,
                "residue": 258,
                "residue_indices_sha256": stable_hash(k7_indices),
                "adjacency_records_checked": len(k7_graphs),
                "adjacency_records_sha256": stable_hash(k7_graphs),
            },
            "K6_only": {
                "input": 990,
                "exact_rejections": 13,
                "residue": 977,
                "residue_indices_sha256": stable_hash(k6_indices),
                "adjacency_records_checked": len(k6_graphs),
                "adjacency_records_sha256": stable_hash(k6_graphs),
            },
        },
        "combined": expected_combined,
        "checks": {
            "K7_union_partition": True,
            "K7_independent_verification_gate": True,
            "K6_ordered_complement": True,
            "K6_independent_verification_gate": True,
            "K6_positive_control": True,
            "class_disjointness": True,
            "embedded_adjacencies_match_sources": True,
            "all_adjacencies_symmetric_loopless_19_vertex": True,
            "source_hashes": True,
        },
        "inputs": expected_hashes,
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
    }


def atomic_json(output: Path, value: object) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v2.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v2_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.manifest.resolve())
    output = args.output.resolve()
    atomic_json(output, result)
    print(json.dumps({
        "status": result["status"],
        "K7": result["classes"]["K7"]["residue"],
        "K6_only": result["classes"]["K6_only"]["residue"],
        "combined": result["combined"]["count"],
        "output": str(output),
        "sha256": sha256(output),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
