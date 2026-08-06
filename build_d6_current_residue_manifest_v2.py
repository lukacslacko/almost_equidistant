#!/usr/bin/env python3
"""Build the self-contained, hash-pinned dimension-six residue manifest v2.

This is a focused boundary update, not a replay of the historical campaign.
It intersects two independently verified exact boundaries:

* the 258-graph complement of the exact K7 tetrad/pattern-954 union;
* the 977-graph complement of the exact fused K6-only rejection set.

The old ``d6_current_residue_manifest.json`` remains frozen because several
production tools pin its bytes.  This version embeds the surviving adjacency
records so downstream work does not depend on the ignored ``.runs`` source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parent
BOUNDARY_COMMIT = "fc51874458196c0391530d0e98968f530ab738c9"
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
EXPECTED_K7_INPUT = 12_839
EXPECTED_K7_REJECTED = 12_581
EXPECTED_K7_RESIDUE = 258
EXPECTED_K6_INPUT = 990
EXPECTED_K6_REJECTED = 13
EXPECTED_K6_RESIDUE = 977
EXPECTED_COMBINED = 1_235


def source_path(name: str) -> Path:
    return ROOT / name


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


def read_json(name: str) -> dict:
    value = json.loads(source_path(name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{name} is not a JSON object")
    return value


def integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    if len(value) != len(set(value)):
        raise ValueError(f"{label} contains duplicate indices")
    return list(value)


def verify_source_hashes() -> dict[str, str]:
    observed = {name: sha256(source_path(name)) for name in EXPECTED_SOURCES}
    if observed != EXPECTED_SOURCES:
        raise ValueError(
            f"source hash mismatch: observed={observed}, expected={EXPECTED_SOURCES}"
        )
    return observed


def validate_adjacency(adjacency: object, label: str) -> list[int]:
    if (
        not isinstance(adjacency, list)
        or len(adjacency) != 19
        or any(type(mask) is not int for mask in adjacency)
    ):
        raise ValueError(f"{label} is not a 19-row integer adjacency")
    bound = 1 << 19
    for vertex, mask in enumerate(adjacency):
        if mask < 0 or mask >= bound or mask & (1 << vertex):
            raise ValueError(f"{label} has an invalid row at vertex {vertex}")
        for other in range(19):
            if ((mask >> other) & 1) != ((adjacency[other] >> vertex) & 1):
                raise ValueError(f"{label} is asymmetric at ({vertex}, {other})")
    return list(adjacency)


def adjacency_map(records: object, label: str) -> tuple[list[int], dict[int, list[int]]]:
    if not isinstance(records, list):
        raise ValueError(f"{label} graphs are not a list")
    order: list[int] = []
    result: dict[int, list[int]] = {}
    for position, record in enumerate(records):
        if not isinstance(record, dict) or type(record.get("index")) is not int:
            raise ValueError(f"{label} record {position} has no integer index")
        index = record["index"]
        if index in result:
            raise ValueError(f"{label} repeats index {index}")
        order.append(index)
        result[index] = validate_adjacency(
            record.get("adjacency"), f"{label} index {index}"
        )
    return order, result


def records_for(indices: Sequence[int], lookup: dict[int, list[int]], label: str) -> list[dict]:
    missing = [index for index in indices if index not in lookup]
    if missing:
        raise ValueError(f"{label} adjacency lookup misses {missing[:5]}")
    return [{"index": index, "adjacency": lookup[index]} for index in indices]


def checked_set(union: dict, name: str) -> tuple[list[int], set[int]]:
    payload = union.get("sets", {}).get(name)
    if not isinstance(payload, dict):
        raise ValueError(f"K7 union has no {name} set")
    indices = integer_list(payload.get("indices"), f"K7 {name}")
    if payload.get("count") != len(indices) or payload.get("indices_sha256") != stable_hash(indices):
        raise ValueError(f"K7 {name} count/hash mismatch")
    return indices, set(indices)


def build_manifest() -> dict:
    observed = verify_source_hashes()

    old = read_json("d6_current_residue_manifest.json")
    old_verification = read_json("d6_current_residue_manifest_verification.json")
    if (
        old.get("schema") != "d6-current-exact-residue-v1"
        or old.get("status") != "COMPLETE"
        or old_verification.get("schema")
        != "d6-current-exact-residue-verification-v1"
        or old_verification.get("status") != "PASS"
        or old_verification.get("combined") != 13_829
    ):
        raise ValueError("frozen current-residue gate failed")

    # K7: independently verified exact union and its ordered complement.
    selection = read_json("d6_k7_rankone_tetrad_full_selection.json")
    selected = integer_list(selection.get("selected_indices"), "K7 selection")
    if (
        len(selected) != EXPECTED_K7_INPUT
        or selection.get("selected_indices_sha256") != stable_hash(selected)
    ):
        raise ValueError("K7 selection count/hash mismatch")
    union = read_json("d6_k7_rankone_pattern_union.json")
    union_verification = read_json("d6_k7_rankone_pattern_union_verification.json")
    if (
        union.get("schema") != "d6-k7-rankone-pattern-union-v1"
        or union.get("status") != "COMPLETE"
        or union_verification.get("schema")
        != "d6-k7-rankone-pattern-union-verification-v1"
        or union_verification.get("status") != "PASS"
        or union_verification.get("manifest", {}).get("sha256")
        != observed["d6_k7_rankone_pattern_union.json"]
        or union_verification.get("graphs") != EXPECTED_K7_INPUT
        or union_verification.get("exact_union") != EXPECTED_K7_REJECTED
        or union_verification.get("exact_residue") != EXPECTED_K7_RESIDUE
    ):
        raise ValueError("K7 independent-verification gate failed")
    exact_union, exact_union_set = checked_set(union, "exact_union")
    k7_residue, k7_residue_set = checked_set(union, "exact_residue")
    tetrad, tetrad_set = checked_set(union, "tetrad_rejected")
    pattern, pattern_set = checked_set(union, "pattern_954_rejected")
    both, both_set = checked_set(union, "both")
    tetrad_only, tetrad_only_set = checked_set(union, "tetrad_only")
    pattern_only, pattern_only_set = checked_set(union, "pattern_954_only")
    selected_set = set(selected)
    if (
        len(exact_union) != EXPECTED_K7_REJECTED
        or len(k7_residue) != EXPECTED_K7_RESIDUE
        or exact_union_set != tetrad_set | pattern_set
        or both_set != tetrad_set & pattern_set
        or tetrad_only_set != tetrad_set - pattern_set
        or pattern_only_set != pattern_set - tetrad_set
        or exact_union_set & k7_residue_set
        or exact_union_set | k7_residue_set != selected_set
        or k7_residue != [index for index in selected if index not in exact_union_set]
        or union_verification.get("residue_indices_sha256") != stable_hash(k7_residue)
    ):
        raise ValueError("K7 exact-union partition mismatch")
    k7_origin = read_json(".runs/d6_k7_rank_survivors.json")
    _, k7_lookup = adjacency_map(k7_origin.get("graphs"), "K7 rank source")
    if not selected_set <= set(k7_lookup):
        raise ValueError("K7 selection is not contained in its adjacency source")
    k7_graphs = records_for(k7_residue, k7_lookup, "K7 residue")

    # K6-only: filter the frozen 990-order boundary by the independently
    # verified fused exact rejection set.
    k6_base = integer_list(
        old.get("classes", {}).get("K6_only", {}).get("final_indices"),
        "K6-only frozen boundary",
    )
    if len(k6_base) != EXPECTED_K6_INPUT or stable_hash(k6_base) != (
        "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
    ):
        raise ValueError("K6-only frozen boundary mismatch")
    k6_report = read_json("d6_k6_fused_side_rank_report.json")
    k6_verification = read_json("d6_k6_fused_side_rank_verification.json")
    rejected = integer_list(k6_report.get("rejected_indices"), "K6 fused rejections")
    verified_rejected = integer_list(
        k6_verification.get("rejected_indices"), "K6 verified fused rejections"
    )
    report_order = [
        row.get("index") for row in k6_report.get("graph_results", [])
        if isinstance(row, dict)
    ]
    if (
        k6_report.get("schema") != "d6-k6-fused-side-rank-v1"
        or k6_report.get("status") != "COMPLETE"
        or k6_verification.get("schema")
        != "d6-k6-fused-side-rank-verification-v1"
        or k6_verification.get("status") != "PASS"
        or k6_verification.get("inputs", {}).get("d6_k6_fused_side_rank_report.json")
        != observed["d6_k6_fused_side_rank_report.json"]
        or k6_report.get("input_graphs") != EXPECTED_K6_INPUT
        or k6_report.get("input_indices_sha256") != stable_hash(k6_base)
        or report_order != k6_base
        or k6_report.get("graphs_rejected") != EXPECTED_K6_REJECTED
        or k6_report.get("graphs_surviving") != EXPECTED_K6_RESIDUE
        or k6_verification.get("graphs_recomputed") != EXPECTED_K6_INPUT
        or k6_verification.get("graphs_rejected") != EXPECTED_K6_REJECTED
        or k6_verification.get("graphs_surviving") != EXPECTED_K6_RESIDUE
        or rejected != verified_rejected
        or k6_report.get("positive_control", {}).get("passed") is not True
        or k6_verification.get("positive_control", {}).get("passed") is not True
    ):
        raise ValueError("K6 independent-verification gate failed")
    rejected_set = set(rejected)
    if len(rejected_set) != EXPECTED_K6_REJECTED or not rejected_set <= set(k6_base):
        raise ValueError("K6 fused rejection set mismatch")
    k6_residue = [index for index in k6_base if index not in rejected_set]
    if len(k6_residue) != EXPECTED_K6_RESIDUE:
        raise ValueError("K6 residue count mismatch")
    k6_input = read_json("d6_k6_bipartite_rank_input.json")
    _, k6_lookup = adjacency_map(k6_input.get("graphs"), "K6 adjacency source")
    if not set(k6_base) <= set(k6_lookup):
        raise ValueError("K6 boundary is not contained in its adjacency source")
    k6_graphs = records_for(k6_residue, k6_lookup, "K6 residue")

    if selected_set & set(k6_lookup):
        # The sources cover distinct graph classes, not merely distinct final sets.
        raise ValueError("K7 and K6-only adjacency populations overlap")
    if k7_residue_set & set(k6_residue):
        raise ValueError("K7 and K6-only residues overlap")

    class_index_records = (
        [{"class": "K7", "index": index} for index in k7_residue]
        + [{"class": "K6_only", "index": index} for index in k6_residue]
    )
    class_graph_records = (
        [{"class": "K7", **record} for record in k7_graphs]
        + [{"class": "K6_only", **record} for record in k6_graphs]
    )
    concatenated_indices = k7_residue + k6_residue
    sources = dict(observed)
    sources[Path(__file__).name] = sha256(Path(__file__).resolve())

    return {
        "schema": "d6-current-exact-residue-v2",
        "status": "COMPLETE",
        "boundary_commit": BOUNDARY_COMMIT,
        "description": (
            "Self-contained exact unresolved dimension-six level-19 boundary "
            "after the independently verified K7 tetrad/pattern union and "
            "K6-only fused side-rank layers."
        ),
        "class_order": ["K7", "K6_only"],
        "classes": {
            "K7": {
                "input_boundary_count": EXPECTED_K7_INPUT,
                "input_boundary_indices_sha256": stable_hash(selected),
                "exact_rejections": EXPECTED_K7_REJECTED,
                "exact_rejected_indices_sha256": stable_hash(exact_union),
                "residue_count": len(k7_residue),
                "residue_indices": k7_residue,
                "residue_indices_sha256": stable_hash(k7_residue),
                "graphs": k7_graphs,
                "graphs_sha256": stable_hash(k7_graphs),
                "adjacency_origin": {
                    "path": ".runs/d6_k7_rank_survivors.json",
                    "sha256": observed[".runs/d6_k7_rank_survivors.json"],
                    "source_graphs": len(k7_lookup),
                    "embedded": True,
                },
            },
            "K6_only": {
                "input_boundary_count": EXPECTED_K6_INPUT,
                "input_boundary_indices_sha256": stable_hash(k6_base),
                "exact_rejections": EXPECTED_K6_REJECTED,
                "exact_rejected_indices": rejected,
                "exact_rejected_indices_sha256": stable_hash(rejected),
                "residue_count": len(k6_residue),
                "residue_indices": k6_residue,
                "residue_indices_sha256": stable_hash(k6_residue),
                "graphs": k6_graphs,
                "graphs_sha256": stable_hash(k6_graphs),
                "adjacency_origin": {
                    "path": "d6_k6_bipartite_rank_input.json",
                    "sha256": observed["d6_k6_bipartite_rank_input.json"],
                    "source_graphs": len(k6_lookup),
                    "embedded": True,
                },
            },
        },
        "combined": {
            "count": len(class_index_records),
            "class_counts": {"K7": len(k7_residue), "K6_only": len(k6_residue)},
            "ordered_indices_sha256": stable_hash(concatenated_indices),
            "sorted_indices_sha256": stable_hash(sorted(concatenated_indices)),
            "ordered_class_index_records_sha256": stable_hash(class_index_records),
            "ordered_class_graph_records_sha256": stable_hash(class_graph_records),
            "cross_class_overlap": 0,
            "cross_class_overlap_indices_sha256": stable_hash([]),
        },
        "verification_gates": {
            "prior_boundary": {
                "manifest": "d6_current_residue_manifest.json",
                "manifest_sha256": observed["d6_current_residue_manifest.json"],
                "verification": "d6_current_residue_manifest_verification.json",
                "verification_sha256": observed[
                    "d6_current_residue_manifest_verification.json"
                ],
                "verification_status": "PASS",
            },
            "K7": {
                "manifest": "d6_k7_rankone_pattern_union.json",
                "manifest_sha256": observed["d6_k7_rankone_pattern_union.json"],
                "verification": "d6_k7_rankone_pattern_union_verification.json",
                "verification_sha256": observed[
                    "d6_k7_rankone_pattern_union_verification.json"
                ],
                "verification_status": "PASS",
                "exact_union": EXPECTED_K7_REJECTED,
                "exact_residue": EXPECTED_K7_RESIDUE,
            },
            "K6_only": {
                "report": "d6_k6_fused_side_rank_report.json",
                "report_sha256": observed["d6_k6_fused_side_rank_report.json"],
                "verification": "d6_k6_fused_side_rank_verification.json",
                "verification_sha256": observed[
                    "d6_k6_fused_side_rank_verification.json"
                ],
                "verification_status": "PASS",
                "positive_control_passed": True,
                "exact_rejections": EXPECTED_K6_REJECTED,
                "exact_residue": EXPECTED_K6_RESIDUE,
            },
        },
        "source_hashes": sources,
        "semantics": {
            "unit_graph_edges": "required to have Euclidean distance exactly 1",
            "candidate_nonedges": "unconstrained and may also have distance 1",
            "points": "must be distinct",
            "allowed_defect_coordinates": "upper bounds only and may be zero",
            "membership_arithmetic": (
                "exact set filtering and integer bitmask copying; no floating-point decision"
            ),
            "residue": "filter non-rejection only; not a realizability claim",
        },
        "nonclaims": [
            "No graph in the residue is claimed realizable in R6.",
            "No candidate nonedge is constrained to be non-unit.",
            "The combined count 1235 is not yet a proof that f(6)=18.",
            "The manifest does not replace the upstream exact rejection certificates.",
        ],
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
        default=ROOT / "d6_current_residue_manifest_v2.json",
    )
    args = parser.parse_args()
    manifest = build_manifest()
    output = args.output.resolve()
    atomic_json(output, manifest)
    print(json.dumps({
        "status": manifest["status"],
        "K7": manifest["classes"]["K7"]["residue_count"],
        "K6_only": manifest["classes"]["K6_only"]["residue_count"],
        "combined": manifest["combined"]["count"],
        "output": str(output),
        "sha256": sha256(output),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
