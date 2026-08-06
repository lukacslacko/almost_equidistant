#!/usr/bin/env python3
"""Build the hash-pinned current K7 residue for the full dual campaign.

The selection is reconstructed from the complete 17,764-graph rank-survivor
corpus.  It removes the union of every *already certified* exact layer that is
available at this commit: support propagation, strict-H, sparse values, the
cap-20,000 interval benchmark, and the frozen positive-polynomial sample.

This program is intentionally a small, independently runnable set-accounting
checker.  A changed source artifact is rejected until its hash and expected
accounting have been deliberately audited and updated here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_INPUT_GRAPHS = 17_764
EXPECTED_SOURCE_SHA256 = {
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
}
EXPECTED_INCREMENTAL_COUNTS = {
    "support_propagation": 1_536,
    "strict_H": 132,
    "sparse_value": 3_148,
    "interval_cap20000": 6,
    "sample_positive_dual": 1,
}
EXPECTED_SELECTED = 12_941
EXPECTED_SELECTED_INDICES_SHA256 = (
    "a07bafe688c45ff132921519320a19687d2917a79a6c9009bda51a19a6b7f290"
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


def read_pinned(name: str) -> dict:
    path = ROOT / name
    actual = sha256(path)
    if actual != EXPECTED_SOURCE_SHA256[name]:
        raise ValueError(f"source hash mismatch for {name}: {actual}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"source {name} is not a JSON object")
    return value


def integer_set(value: object, label: str) -> set[int]:
    if not isinstance(value, list) or not all(
        isinstance(item, int) and not isinstance(item, bool) for item in value
    ):
        raise ValueError(f"{label} is not an integer list")
    output = set(value)
    if len(output) != len(value):
        raise ValueError(f"{label} contains duplicates")
    return output


def source_rejection_sets(input_hash: str) -> list[tuple[str, set[int], dict]]:
    support = read_pinned("d6_k7_support_rank_survivors_report.json")
    if (
        support.get("input_sha256") != input_hash
        or support.get("graphs") != EXPECTED_INPUT_GRAPHS
        or support.get("refined_rejected") != 1_536
    ):
        raise ValueError("support report binding/count mismatch")
    support_set = integer_set(
        support.get("refined_rejected_indices"), "support rejections"
    )
    if len(support_set) != support["refined_rejected"]:
        raise ValueError("support rejection list/count mismatch")

    strict_h = read_pinned("d6_k7_strict_h_full_report.json")
    if (
        strict_h.get("graphs") != EXPECTED_INPUT_GRAPHS
        or strict_h.get("configuration", {}).get("graphs")
        != EXPECTED_INPUT_GRAPHS
        or strict_h.get("counts", {}).get("strict_H_rejected") != 603
    ):
        raise ValueError("strict-H report binding/count mismatch")
    strict_set = integer_set(
        strict_h.get("decisions", {}).get("strict_H_rejected"),
        "strict-H rejections",
    )
    if len(strict_set) != strict_h["counts"]["strict_H_rejected"]:
        raise ValueError("strict-H rejection list/count mismatch")

    sparse = read_pinned("d6_k7_small_support_value_full_report.json")
    if (
        sparse.get("selection", {}).get("rank_survivors")
        != EXPECTED_INPUT_GRAPHS
        or sparse.get("selection", {}).get("support_rejected") != 1_536
        or sparse.get("decisions", {}).get("graphs") != 16_228
        or sparse.get("decisions", {}).get("rejected") != 3_195
    ):
        raise ValueError("sparse-value report binding/count mismatch")
    sparse_set = integer_set(
        sparse.get("decisions", {}).get("rejected_indices"),
        "sparse-value rejections",
    )
    if len(sparse_set) != sparse["decisions"]["rejected"]:
        raise ValueError("sparse-value rejection list/count mismatch")

    interval = read_pinned("d6_interval_residue_benchmark.json")
    cap_records = [
        record for record in interval.get("benchmarks", [])
        if record.get("cap") == 20_000
    ]
    if (
        len(cap_records) != 1
        or cap_records[0].get("certified_killed") != 10
        or cap_records[0].get("status_counts", {}).get("INFRA_ERROR") != 0
        or interval.get("checks", {}).get("corpus_and_selection_reconstruction")
        != "PASS"
    ):
        raise ValueError("cap-20,000 interval benchmark binding/count mismatch")
    interval_set = integer_set(
        cap_records[0].get("killed_indices"), "cap-20,000 interval rejections"
    )
    if len(interval_set) != cap_records[0]["certified_killed"]:
        raise ValueError("interval rejection list/count mismatch")

    sample_dual = read_pinned("d6_k7_positive_polynomial_dual_report.json")
    if (
        sample_dual.get("graphs") != 84
        or sample_dual.get("counts", {}).get("marginal_dual_rejected") != 2
    ):
        raise ValueError("sample positive-dual report binding/count mismatch")
    sample_set = integer_set(
        sample_dual.get("decisions", {}).get("marginal_dual_rejected"),
        "sample positive-dual rejections",
    )
    if len(sample_set) != sample_dual["counts"]["marginal_dual_rejected"]:
        raise ValueError("sample positive-dual rejection list/count mismatch")

    return [
        ("support_propagation", support_set, {
            "reported_rejections": support["refined_rejected"]
        }),
        ("strict_H", strict_set, {
            "reported_rejections": strict_h["counts"]["strict_H_rejected"]
        }),
        ("sparse_value", sparse_set, {
            "reported_rejections": sparse["decisions"]["rejected"]
        }),
        ("interval_cap20000", interval_set, {
            "reported_rejections": cap_records[0]["certified_killed"],
            "benchmark_scope_graphs": interval.get("scope", {}).get(
                "sample_graphs"
            ),
        }),
        ("sample_positive_dual", sample_set, {
            "reported_rejections": sample_dual["counts"][
                "marginal_dual_rejected"
            ]
        }),
    ]


def build_selection(input_path: Path) -> dict:
    input_hash = sha256(input_path)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"rank-survivor input hash mismatch: {input_hash}")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list) or len(graphs) != EXPECTED_INPUT_GRAPHS:
        raise ValueError("rank-survivor corpus graph count mismatch")
    indices = [int(graph["index"]) for graph in graphs]
    universe = set(indices)
    if len(universe) != len(indices):
        raise ValueError("rank-survivor corpus repeats an index")

    excluded: set[int] = set()
    accounting: list[dict] = []
    for name, raw_rejections, metadata in source_rejection_sets(input_hash):
        in_universe = raw_rejections & universe
        incremental = in_universe - excluded
        if len(incremental) != EXPECTED_INCREMENTAL_COUNTS[name]:
            raise ValueError(
                f"{name} incremental count {len(incremental)} != "
                f"{EXPECTED_INCREMENTAL_COUNTS[name]}"
            )
        excluded.update(in_universe)
        record = {
            "layer": name,
            **metadata,
            "rejections_in_rank_survivors": len(in_universe),
            "rejections_outside_rank_survivors": len(raw_rejections - universe),
            "overlap_with_prior_layers": len(in_universe - incremental),
            "incremental_rejections": len(incremental),
            "rejections_in_rank_survivors_sha256": stable_hash(
                sorted(in_universe)
            ),
            "incremental_rejections_sha256": stable_hash(sorted(incremental)),
            "cumulative_rejections": len(excluded),
            "cumulative_survivors": len(universe - excluded),
        }
        if len(incremental) <= 64:
            record["incremental_rejected_indices"] = sorted(incremental)
        accounting.append(record)

    selected = [index for index in indices if index not in excluded]
    selected_hash = stable_hash(selected)
    if len(selected) != EXPECTED_SELECTED:
        raise ValueError(f"selected count {len(selected)} != {EXPECTED_SELECTED}")
    if selected_hash != EXPECTED_SELECTED_INDICES_SHA256:
        raise ValueError(f"selected-index hash mismatch: {selected_hash}")
    source_artifacts = [
        {"path": name, "sha256": digest}
        for name, digest in EXPECTED_SOURCE_SHA256.items()
    ]
    source_artifacts.append({
        "path": Path(__file__).name,
        "sha256": sha256(Path(__file__).resolve()),
    })
    try:
        input_display = str(input_path.resolve().relative_to(ROOT))
    except ValueError:
        input_display = str(input_path.resolve())
    return {
        "schema": 1,
        "kind": "d6_k7_positive_dual_selection",
        "description": (
            "Current K7 exact residue after support, strict-H, sparse-value, "
            "cap-20,000 interval-benchmark, and frozen sample-dual rejections."
        ),
        "input": input_display,
        "input_sha256": input_hash,
        "rank_survivors": len(indices),
        "layer_accounting": accounting,
        "excluded": len(excluded),
        "excluded_indices_sha256": stable_hash(sorted(excluded)),
        "selected": len(selected),
        "selected_indices": selected,
        "selected_indices_sha256": selected_hash,
        "source_artifacts": source_artifacts,
    }


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}"
    )
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(".runs/d6_k7_rank_survivors.json")
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_positive_dual_selection.json",
    )
    args = parser.parse_args()
    report = build_selection(args.input)
    atomic_json(args.output, report)
    print(
        f"selected {report['selected']} K7 graphs; "
        f"manifest sha256 {sha256(args.output)}"
    )


if __name__ == "__main__":
    main()
