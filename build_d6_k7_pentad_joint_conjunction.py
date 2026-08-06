#!/usr/bin/env python3
"""Build the exact seed-local conjunction of K7 cover and pentad filters.

The committed full joint-support campaign rejects 69 of the frozen 258 K7
residue graphs.  The rank-two pentad campaign certifies 36 individual covers.
This builder re-enumerates the exact current-cover quantifier and the full
propagation/sparse-value support layer from the 258 adjacencies embedded in
``d6_current_residue_manifest_v2.json``.  It records a new graph rejection
only when one K7 seed has a nonempty current-cover set and *every* member of
that set either has no support family or has a pentad certificate.  The more
conservative all-covers-have-pentads subset is reported separately.

No ignored ``.runs`` file is read.  In particular, the old aggregate
``has_no_near_clique_cover`` profile is not used for the theorem claim.
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
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import run_d6_k7_support_capacity_pilot as joint_kernel
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
    "run_d6_k7_support_capacity_pilot.py": (
        "9af542217d178bec2a71cb4c30faabc3dea8c8279053d729d924d9ca35eb128a"
    ),
    "d6_k7_support_propagation.py": (
        "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
    ),
    "d6_k7_small_support_value.py": (
        "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f"
    ),
    "d6_k7_support_capacity.py": (
        "df8c010faac37f9de5481316cc8ccf58e8df46230bd9d1e8542f7a4322b64996"
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
                    raise ValueError(f"duplicate {field} key at graph {index}")
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
        raise ValueError("candidate graph missing from tetrad decision archive")
    return output


def read_joint_decisions(indices: Sequence[int]) -> dict[int, dict]:
    with gzip.open(JOINT_DECISIONS, "rt", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != len(indices):
        raise ValueError("joint decision archive population mismatch")
    output = {}
    for ordinal, (index, row) in enumerate(zip(indices, rows)):
        if int(row["ordinal"]) != ordinal or int(row["index"]) != index:
            raise ValueError("joint decision archive order mismatch")
        if index in output:
            raise ValueError("duplicate joint decision")
        output[index] = row
    return output


def quantify_graph(payload: tuple[dict, list, list, dict]) -> dict:
    """Replay inherited covers plus exact propagation/sparse-value support."""

    graph, serialized_prior, serialized_tetrad, tetrad_row = payload
    prior_failures = {
        (tuple(int(vertex) for vertex in seed), int(zmask))
        for seed, zmask in serialized_prior
    }
    tetrad_failures = {
        (tuple(int(vertex) for vertex in seed), int(zmask))
        for seed, zmask in serialized_tetrad
    }
    result = joint_kernel.analyze_graph(
        graph, prior_failures, tetrad_failures, tetrad_row, node_limit=0
    )
    index = int(result["index"])
    seed_quantifiers = []
    for seed_record in result["seed_records"]:
        current = []
        for cover in seed_record["current_cover_records"]:
            pre_capacity_pass = bool(cover["pre_capacity_pass"])
            family_counts = {
                name: int(value)
                for name, value in cover["counts"].items()
                if name == "labeled_z_families"
                or name == "propagation_passing_families"
                or name == "sparse_value_failing_families"
                or name == "pre_capacity_passing_families"
                or name.startswith("propagation_failure:")
            }
            current.append({
                "seed": list(seed_record["seed"]),
                "zmask": int(cover["zmask"]),
                "joint_pre_capacity_status": (
                    "PASSING" if pre_capacity_pass else "INFEASIBLE"
                ),
                "family_counts": family_counts,
            })
        seed_quantifiers.append({
            "seed": list(seed_record["seed"]),
            "seed_mask": int(seed_record["seed_mask"]),
            "eligible_covers": int(seed_record["counts"].get("covers", 0)),
            "cover_status_counts": {
                name.removeprefix("cover_"): int(value)
                for name, value in seed_record["counts"].items()
                if name.startswith("cover_")
            },
            "joint_pre_capacity_failing_covers": int(
                seed_record["counts"].get(
                    "joint_pre_capacity_failing_covers", 0
                )
            ),
            "pre_capacity_passing_covers": int(
                seed_record["counts"].get("pre_capacity_passing_covers", 0)
            ),
            "current_passing_covers": current,
        })
    return {
        "index": index,
        "counts": result["counts"],
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


def validate_upstreams() -> tuple[
    list[int], dict[int, dict], dict[int, dict], dict, dict
]:
    for name, expected in {**EXPECTED_HASHES, **EXPECTED_SOURCE_HASHES}.items():
        observed = sha256(ROOT / name)
        if observed != expected:
            raise ValueError(f"hash mismatch for {name}: {observed} != {expected}")

    current = json.loads(CURRENT_RESIDUE.read_text(encoding="utf-8"))
    if current.get("schema") != "d6-current-exact-residue-v2" or (
        current.get("status") != "COMPLETE"
    ):
        raise ValueError("current residue manifest is not complete v2")
    k7 = current["classes"]["K7"]
    graphs = list(k7["graphs"])
    indices = [int(index) for index in k7["residue_indices"]]
    if (
        len(graphs) != 258
        or [int(graph["index"]) for graph in graphs] != indices
        or stable_hash(indices) != k7["residue_indices_sha256"]
        or stable_hash(graphs) != k7["graphs_sha256"]
    ):
        raise ValueError("embedded K7 adjacency corpus does not bind")
    graph_by_index = {int(graph["index"]): graph for graph in graphs}
    if len(graph_by_index) != 258:
        raise ValueError("embedded K7 adjacency corpus repeats an index")

    joint = json.loads(JOINT_REPORT.read_text(encoding="utf-8"))
    joint_verification = json.loads(JOINT_VERIFICATION.read_text(encoding="utf-8"))
    if (
        joint.get("schema") != 1 or joint.get("status") != "COMPLETE"
        or joint_verification.get("status") != "PASS"
        or joint_verification["report"]["sha256"] != EXPECTED_HASHES[
            "d6_k7_joint_support_full_report.json"
        ]
        or joint["selection"]["indices_sha256"] != stable_hash(indices)
    ):
        raise ValueError("full joint-support verification gate failed")
    for field, name in (
        ("checkpoint_copy_sha256", "d6_k7_joint_support_full_checkpoint.json"),
        ("decisions_sha256", "d6_k7_joint_support_full_decisions.tsv.gz"),
        ("certificates_sha256", "d6_k7_joint_support_full_certificates.jsonl.gz"),
    ):
        if joint["artifacts"][field] != EXPECTED_HASHES[name]:
            raise ValueError(f"joint report artifact binding mismatch: {name}")
    checkpoint = json.loads(JOINT_CHECKPOINT.read_text(encoding="utf-8"))
    if (
        checkpoint.get("status") != "COMPLETE"
        or int(checkpoint.get("completed", -1)) != 258
        or checkpoint.get("config") != joint.get("configuration")
        or checkpoint.get("config_sha256") != joint.get("configuration_sha256")
    ):
        raise ValueError("joint checkpoint-copy binding failed")
    joint_decisions = read_joint_decisions(indices)
    joint_rejected = [
        index for index in indices
        if joint_decisions[index]["status"] == "JOINT_PRE_CAPACITY_REJECTED"
    ]
    if (
        joint_rejected != joint["summary"]["joint_rejected_indices"]
        or len(joint_rejected) != 69
        or joint_verification.get("joint_rejections_verified") != 69
    ):
        raise ValueError("joint rejection population mismatch")
    with gzip.open(JOINT_CERTIFICATES, "rt", encoding="utf-8") as stream:
        joint_certificates = [json.loads(line) for line in stream if line.strip()]
    if (
        [int(item["index"]) for item in joint_certificates] != joint_rejected
        or any(
            item.get("kind") != "joint_k7_cover_support_rejection"
            for item in joint_certificates
        )
    ):
        raise ValueError("joint certificate archive does not bind rejection list")

    pentad = json.loads(PENTAD_REPORT.read_text(encoding="utf-8"))
    pentad_verification = json.loads(
        PENTAD_VERIFICATION.read_text(encoding="utf-8")
    )
    if (
        pentad.get("schema") != 1 or pentad.get("status") != "COMPLETE"
        or pentad_verification.get("status") != "PASS"
        or pentad_verification["report"]["sha256"] != EXPECTED_HASHES[
            "d6_k7_ranktwo_pentad_full_report.json"
        ]
        or pentad_verification["checked"]["certificates"] != 36
    ):
        raise ValueError("full pentad verification gate failed")
    return indices, graph_by_index, joint_decisions, joint, pentad


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_pentad_joint_conjunction_report.json",
    )
    args = parser.parse_args()
    if args.workers <= 0:
        raise SystemExit("workers must be positive")

    started_at = datetime.now(UTC).isoformat()
    started = time.monotonic()
    indices, graph_by_index, joint_rows, joint_report, pentad_report = (
        validate_upstreams()
    )
    selected = set(indices)
    joint_rejected = [
        int(index) for index in joint_report["summary"]["joint_rejected_indices"]
    ]
    joint_rejected_set = set(joint_rejected)

    pentad_by_key = {}
    rejected_key_order = []
    for record in pentad_report["covers"]:
        if record["status"] != "REJECTED":
            continue
        key = key_tuple(record)
        if key in pentad_by_key or key[0] not in selected:
            raise ValueError("duplicate or out-of-corpus pentad rejection")
        if not record.get("certificate"):
            raise ValueError("pentad rejection lacks its exact certificate")
        pentad_by_key[key] = record
        rejected_key_order.append(key)
    if len(pentad_by_key) != 36:
        raise ValueError("expected 36 pentad-rejected covers")
    candidate_indices = [
        index for index in indices
        if any(key[0] == index for key in pentad_by_key)
        and index not in joint_rejected_set
    ]
    if len(candidate_indices) != 34:
        raise ValueError("expected 34 pentad candidate graphs outside joint set")

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
            tetrad_rows[index],
        )
        for index in indices
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        quantified = list(executor.map(quantify_graph, payloads, chunksize=1))
    if [int(item["index"]) for item in quantified] != indices:
        raise ValueError("parallel quantifier replay changed graph order")

    quantified_by_index = {int(item["index"]): item for item in quantified}
    graph_summaries = []
    global_counts: Counter = Counter()
    all_current_keys = set()
    for item in quantified:
        index = int(item["index"])
        row = joint_rows[index]
        counts = item["counts"]
        observed_counts = {
            "seeds": int(counts.get("seeds", 0)),
            "covers": int(counts.get("covers", 0)),
            "current_passing_covers": int(counts.get("cover_passing", 0)),
            "joint_pre_capacity_failing_covers": int(
                counts.get("joint_pre_capacity_failing_covers", 0)
            ),
            "pre_capacity_passing_covers": int(
                counts.get("pre_capacity_passing_covers", 0)
            ),
        }
        for name, observed in observed_counts.items():
            if int(row[name]) != observed:
                raise ValueError(
                    f"graph {index} {name}: {observed} != archive {row[name]}"
                )
            global_counts[name] += observed
        if int(tetrad_rows[index]["tetrad_passing_covers"]) != (
            observed_counts["current_passing_covers"]
        ):
            raise ValueError(f"graph {index} tetrad-passing count mismatch")
        graph_summaries.append({"index": index, **observed_counts})
        for seed_record in item["seed_quantifiers"]:
            for cover in seed_record["current_passing_covers"]:
                all_current_keys.add((
                    index, tuple(cover["seed"]), int(cover["zmask"])
                ))
    expected_cover_totals = {
        "current_passing_covers": 738,
        "joint_pre_capacity_failing_covers": 342,
        "pre_capacity_passing_covers": 396,
    }
    joint_total_fields = {
        "current_passing_covers": "cover_passing",
        "joint_pre_capacity_failing_covers": (
            "joint_pre_capacity_failing_covers"
        ),
        "pre_capacity_passing_covers": "pre_capacity_passing_covers",
    }
    for name, expected in expected_cover_totals.items():
        if global_counts[name] != expected or (
            int(joint_report["summary"]["totals"][joint_total_fields[name]])
            != expected
        ):
            raise ValueError(f"global joint-cover total mismatch in {name}")

    witnesses = []
    conservative_rejected = []
    full_rejected = []
    for index in candidate_indices:
        item = quantified_by_index[index]
        row = joint_rows[index]
        if row["status"] != "SURVIVOR":
            raise ValueError(f"new graph {index} was not a joint survivor")

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
                if pentad_record is not None:
                    entry["pentad_certificate"] = {
                        "report_ordinal": int(pentad_record["ordinal"]),
                        "certificate_sha256": stable_hash(
                            pentad_record["certificate"]
                        ),
                        "pentad_index": int(
                            pentad_record["certificate"]["pentad_index"]
                        ),
                    }
                else:
                    entry["pentad_certificate"] = None
                covers.append(entry)
            record = {
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
                conservative_seeds.append(record)
            if all(
                cover["joint_pre_capacity_status"] == "INFEASIBLE"
                or key in pentad_by_key
                for cover, key in zip(
                    seed_record["current_passing_covers"], current_keys
                )
            ):
                full_seeds.append(record)
        if conservative_seeds:
            conservative_rejected.append(index)
        if full_seeds:
            full_rejected.append(index)
        if not full_seeds:
            raise ValueError(f"graph {index} has no full seed-local conjunction")
        witnesses.append({
            "ordinal": indices.index(index),
            "index": index,
            "archived_graph_counts": {
                "seeds": int(row["seeds"]),
                "covers": int(row["covers"]),
                "current_passing_covers": int(row["current_passing_covers"]),
                "joint_pre_capacity_failing_covers": int(
                    row["joint_pre_capacity_failing_covers"]
                ),
                "pre_capacity_passing_covers": int(
                    row["pre_capacity_passing_covers"]
                ),
            },
            "conservative_pentad_only_seeds": conservative_seeds,
            "full_joint_or_pentad_seeds": full_seeds,
        })

    missing_rejected_keys = set(pentad_by_key) - all_current_keys
    if missing_rejected_keys:
        raise ValueError(
            f"pentad-rejected covers absent from exact current quantifier: "
            f"{sorted(missing_rejected_keys)}"
        )
    if conservative_rejected != full_rejected:
        raise ValueError(
            "this fixed pentad library unexpectedly distinguishes full and "
            "conservative graph sets; inspect before freezing"
        )
    new_rejected = list(full_rejected)
    if len(new_rejected) != 34 or set(new_rejected) & joint_rejected_set:
        raise ValueError("new rejection count/disjointness failed")
    combined_set = joint_rejected_set | set(new_rejected)
    combined = [index for index in indices if index in combined_set]
    residue = [index for index in indices if index not in combined_set]
    if len(combined) != 103 or len(residue) != 155:
        raise ValueError("expected exact 103/155 partition")

    positive = tuple(int(row) for row in lower_bound_18_graph())
    reference.validate_graph(positive)
    positive_k7 = sum(1 for _ in reference.clique_masks(positive, 7))
    if positive_k7:
        raise ValueError("known realizable 18-point control unexpectedly has K7")

    source_hashes = {
        Path(__file__).name: sha256(Path(__file__)),
        **EXPECTED_SOURCE_HASHES,
    }
    report = {
        "schema": 1,
        "kind": "d6_k7_pentad_joint_seed_conjunction",
        "status": "COMPLETE",
        "claim": (
            "The 69 independently verified joint-support rejections and 34 "
            "disjoint seed-local pentad-cover rejections eliminate exactly "
            "103 of the frozen 258 K7 residue graphs. Exactly 155 remain."
        ),
        "semantics": {
            "seed_quantifier": (
                "A new graph is rejected only when one required K7 seed has "
                "a nonempty exact inherited-passing cover set and every such "
                "cover either has no propagation/sparse-value support family "
                "or has an exact coefficientwise-positive rank-two pentad "
                "certificate."
            ),
            "nonedges_optional": True,
            "aggregate_profile_used_for_claim": False,
            "joint_support_capacity_increment_credited": False,
            "unresolved_is_not_rejected": True,
        },
        "upstream_sha256": dict(sorted(EXPECTED_HASHES.items())),
        "source_sha256": source_hashes,
        "selection": {
            "graphs": len(indices),
            "ordered_indices_sha256": stable_hash(indices),
            "embedded_graphs_sha256": stable_hash([
                graph_by_index[index] for index in indices
            ]),
        },
        "joint_cover_quantifier_replay": {
            "graphs": len(graph_summaries),
            "ordered_graph_summaries": graph_summaries,
            "ordered_graph_summaries_sha256": stable_hash(graph_summaries),
            "totals": dict(sorted(global_counts.items())),
            "required_cover_partition": expected_cover_totals,
        },
        "pentad_cover_rejections": {
            "covers": len(rejected_key_order),
            "graphs": len({key[0] for key in rejected_key_order}),
            "ordered_keys": [key_payload(key) for key in rejected_key_order],
            "ordered_keys_sha256": stable_hash([
                key_payload(key) for key in rejected_key_order
            ]),
        },
        "joint_support_rejections": {
            "count": len(joint_rejected),
            "ordered_indices": joint_rejected,
            "ordered_indices_sha256": stable_hash(joint_rejected),
        },
        "new_seed_conjunction_rejections": {
            "count": len(new_rejected),
            "ordered_indices": new_rejected,
            "ordered_indices_sha256": stable_hash(new_rejected),
            "conservative_pentad_only_count": len(conservative_rejected),
            "conservative_pentad_only_ordered_indices": conservative_rejected,
            "conservative_pentad_only_indices_sha256": stable_hash(
                conservative_rejected
            ),
            "full_joint_or_pentad_count": len(full_rejected),
            "full_joint_or_pentad_ordered_indices": full_rejected,
            "full_joint_or_pentad_indices_sha256": stable_hash(full_rejected),
            "witnesses": witnesses,
            "witnesses_sha256": stable_hash(witnesses),
        },
        "combined_rejections": {
            "count": len(combined),
            "ordered_indices": combined,
            "ordered_indices_sha256": stable_hash(combined),
        },
        "exact_residue": {
            "count": len(residue),
            "ordered_indices": residue,
            "ordered_indices_sha256": stable_hash(residue),
        },
        "positive_18_control": {
            "vertices": 18,
            "K7_seeds": positive_k7,
            "status": "PASS_NOT_APPLICABLE_NO_K7",
        },
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
        "joint_rejections": len(joint_rejected),
        "new_rejections": len(new_rejected),
        "combined_rejections": len(combined),
        "residue": len(residue),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
