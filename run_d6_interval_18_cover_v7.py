#!/usr/bin/env python3
"""Run the certified d=6 interval kernel on a dense 18-deletion cover.

The input parents are exactly the independently verified v7 residue.  For
each parent this wrapper keeps all vertex deletions having the maximum number
of required unit edges, then greedily covers the parent set by canonical
deletion classes.  At every greedy tie the lexicographically largest class
identifier wins.  This fully specified rule selects 179 deletion classes for
the 263 v7 parents.  One selected class is a known realizable positive control
and is the sole selected coverer of two parents.  Production retains that
control and deterministically adds one non-positive dense backup for each of
those parents, for 181 deletion classes in total.

A ``KILLED`` deletion class is a sound obstruction for every parent listed in
its coverage record: only required unit edges are passed to the kernel and a
candidate nonedge remains unconstrained (it may also have distance one).
``ABORT``, ``UNRESOLVED``, and ``INFRA_ERROR`` make no claim.  Search and
atomic checkpoint mechanics are delegated unchanged to
``run_d6_interval_residue`` after temporarily setting its graph order to 18.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import run_d6_interval_residue as engine


ROOT = Path(__file__).resolve().parent
N = 18

# Independently checked immutable v7 artifact roots.
EXPECTED_V7_MANIFEST_SHA256 = (
    "1ab0948780d73cdd6aca2925107fc240dcae5df51d7cf5393accd8de3ae79f95"
)
EXPECTED_V7_VERIFICATION_SHA256 = (
    "ff02a9ca5bfe769ce77bbf67dedbc8e928561fb4f291a78d2ea577e06e784255"
)
EXPECTED_V7_COUNTS = {"K7": 12, "K6_only": 251}
EXPECTED_V7_CLASS_HASHES = {
    "K7": "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c",
    "K6_only": "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907",
}
EXPECTED_V7_COMBINED_HASH = (
    "a50baf54bfff642e7c1711fe7041b5184df12f90de2cf52557cd086f10875127"
)

EXPECTED_INPUT_HASHES = {
    "d6_residue_18_deletions_v2.json": (
        "7a13d7a204f866c2d7c421db8a46c5b042bbad8bad6a2cade4a94eb4d340521d"
    ),
    "d6_residue_18_deletions_v2_verification.json": (
        "99927f8b41b48fb2c0151f3f6d04eca37701861cfef9fce4d3bba3f5dab91090"
    ),
    "d6_standard18_geometry_report.json": (
        "97904946b382eb1e6cc325083e1c1ab20e5e8c426fcbaaed1b79357646c846ef"
    ),
    "d6_standard18_geometry_verification.json": (
        "12e3e1df3919dd465a54b1c470f524a3907db2a19110a34876137eeedeb926d6"
    ),
    "d6_18_exact_reconstruction.json": (
        "69fa1c0be5f46e285773f13cc05270219d606673703fd16b3e36a9551cb9a990"
    ),
    "d6_18_exact_reconstruction_verification.json": (
        "e86f8b538f4f78c8a0ace71a13d2e029be054e2135d19ed710faf4fc21aca887"
    ),
}

EXPECTED_PARENT_COUNT = 263
EXPECTED_SELECTED_COUNT = 179
EXPECTED_PRODUCTION_COUNT = 181
EXPECTED_GREEDY_GAIN_HISTOGRAM = {"1": 120, "2": 42, "3": 10, "4": 6, "5": 1}
EXPECTED_SELECTED_POPULATIONS = {"K7": 12, "K6": 167}
EXPECTED_PRODUCTION_POPULATIONS = {"K7": 12, "K6": 169}
EXPECTED_KNOWN_POSITIVE_CLASS_INDEX = 441
EXPECTED_KNOWN_POSITIVE_CLASS_ID = (
    "u18-81842b7afe827273043aa730eaec7a1fda421f033ad25e8b4a126f4fddd456f6"
)
EXPECTED_NO_ORDER_CLASS_INDEX = 7259
EXPECTED_NO_ORDER_CLASS_ID = (
    "u18-d47fb2b0a2f8d5f03f8bc4d37ca035690fd7a33a2c1d60bbc6be8fd912a83416"
)
EXPECTED_ACTIONABLE_BACKUPS = [
    {
        "blocked_parent": "K6_only:364827",
        "class_index": 126,
        "class_id": (
            "u18-f6d7f4c977df6677ea2ceb5682a0686e7f6255ebad113646aaca3eb8ce46cb43"
        ),
    },
    {
        "blocked_parent": "K6_only:3335955",
        "class_index": 5673,
        "class_id": (
            "u18-d3ce53599fbf37b60b6b5bd7f30468280eca97105df1044f097518392fd936d0"
        ),
    },
]
EXPECTED_MINIMUM_CIRCLE_HISTOGRAM = {"1": 123, "2": 57, "unreachable": 1}

# Filled after the deterministic selector below is frozen by its focused
# tests.  Unlike the two v7 artifact sentinels, these are source-level
# regression roots and must never be updated merely to make a test pass.
EXPECTED_SELECTED_CLASS_INDICES_SHA256 = (
    "ef4ac7e73ea924d3d9c6382f590dc73cdf2e65b357e86ddc69719372f717d6fe"
)
EXPECTED_SELECTED_CLASS_IDS_SHA256 = (
    "c3b74892ca75a407fa3b4455546acab5899c851493bdf9cfcd34cb3dcc4827ee"
)
EXPECTED_GREEDY_STEPS_SHA256 = (
    "e91fd4100d2d50bcabb023232b3dcffce52c5bf920fc34d51ee21c3e69b58257"
)
EXPECTED_PRODUCTION_CLASS_INDICES_SHA256 = (
    "8c624c969b517bb5663dc32bdf577ab3addbe8a108458b2803443c17f294d8e6"
)

TRUST_ASSUMPTIONS = {
    "basic_arithmetic": (
        "IEEE-754 binary64 basic operations and sqrt are correctly rounded; "
        "each interval endpoint is expanded with nextafter."
    ),
    "transcendentals": (
        "macOS libm cos endpoint values are assumed within 8 ulps; the "
        "kernel pads both directions by 8 nextafter steps and includes all "
        "interior extrema."
    ),
    "candidate_nonedges": "unconstrained and allowed to be unit distance",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def known_realizable_deletion(record: dict) -> bool:
    """Return the frozen exact-standard positive-control annotation."""

    return bool(record.get("standard18_compatible")) or bool(
        record.get("class_id")
        in {
            "u18-761add09c775ad5ba966a0eef97848e076faacc647744c4079b7bee48a65d92f",
            "u18-84a0a88e385b1aa4673cc149e0c2458f37b7cc6aef3def3cf2589f2075be4841",
        }
    )


def tagged_parent(class_name: str, index: int) -> str:
    return f"{class_name}:{index}"


def parent_sort_key(parent: tuple[str, int]) -> tuple[int, int]:
    return ({"K7": 0, "K6_only": 1}[parent[0]], parent[1])


def validate_adjacency(adjacency: Sequence[int], order: int) -> None:
    require(len(adjacency) == order, f"adjacency order is not {order}")
    full = (1 << order) - 1
    for vertex, row in enumerate(adjacency):
        require(type(row) is int, "adjacency entry is not an integer")
        require(not row & ~full and not row & (1 << vertex), "bad adjacency mask")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(adjacency[other] & (1 << vertex)),
                "asymmetric adjacency",
            )


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path.name} is not a JSON object")
    return value


def load_v7_parents(
    manifest_path: Path,
    verification_path: Path,
    *,
    expected_manifest_sha256: str = EXPECTED_V7_MANIFEST_SHA256,
    expected_verification_sha256: str = EXPECTED_V7_VERIFICATION_SHA256,
) -> tuple[list[dict], dict]:
    """Load exactly the pinned, independently checked v7 parent boundary."""

    require(
        len(expected_manifest_sha256) == 64
        and len(expected_verification_sha256) == 64,
        "v7 artifact hashes are not production-pinned yet",
    )
    manifest_path = manifest_path.resolve()
    verification_path = verification_path.resolve()
    manifest_hash = engine.sha256(manifest_path)
    verification_hash = engine.sha256(verification_path)
    require(manifest_hash == expected_manifest_sha256, "unexpected v7 manifest hash")
    require(
        verification_hash == expected_verification_sha256,
        "unexpected v7 verification hash",
    )
    manifest = load_json(manifest_path)
    verification = load_json(verification_path)
    require(
        manifest.get("schema") == "d6-current-certified-residue-v7"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and manifest.get("class_order") == ["K7", "K6_only"],
        "v7 manifest boundary",
    )
    require(
        verification.get("schema")
        == "d6-current-certified-residue-v7-verification-v1"
        and verification.get("status") == "PASS"
        and verification.get("manifest", {}).get("sha256") == manifest_hash
        and all(verification.get("checks", {}).values()),
        "v7 independent verification boundary",
    )

    parents: list[dict] = []
    ordered_indices: list[int] = []
    seen_indices: set[int] = set()
    for class_name in ("K7", "K6_only"):
        block = manifest.get("classes", {}).get(class_name, {})
        records = block.get("graphs")
        indices = block.get("indices")
        require(isinstance(records, list) and isinstance(indices, list), "v7 class data")
        require(
            len(records) == len(indices) == block.get("count")
            == EXPECTED_V7_COUNTS[class_name],
            f"v7 {class_name} count",
        )
        require(indices == [int(record["index"]) for record in records], "v7 graph order")
        require(
            stable_hash(indices) == block.get("indices_sha256")
            == EXPECTED_V7_CLASS_HASHES[class_name],
            f"v7 {class_name} index root",
        )
        require(stable_hash(records) == block.get("graphs_sha256"), "v7 graph root")
        for record in records:
            index = int(record["index"])
            require(index not in seen_indices, "cross-class v7 index overlap")
            seen_indices.add(index)
            adjacency = tuple(map(int, record["adjacency"]))
            validate_adjacency(adjacency, 19)
            parents.append(
                {"class": class_name, "index": index, "adjacency": adjacency}
            )
            ordered_indices.append(index)

    combined = manifest.get("combined", {})
    require(
        len(parents) == combined.get("count") == EXPECTED_PARENT_COUNT
        and combined.get("class_counts") == EXPECTED_V7_COUNTS
        and combined.get("cross_class_overlap") == 0
        and stable_hash(ordered_indices) == combined.get("ordered_indices_sha256")
        == EXPECTED_V7_COMBINED_HASH,
        "v7 combined boundary",
    )
    require(
        manifest.get("semantics", {}).get("nonedges")
        == "unconstrained and may also have distance one",
        "v7 optional-nonedge semantics",
    )
    return parents, {
        "kind": "independently_verified_d6_current_residue_v7",
        "manifest": {"path": str(manifest_path), "sha256": manifest_hash},
        "verification": {
            "path": str(verification_path),
            "sha256": verification_hash,
            "status": "PASS",
        },
        "counts": EXPECTED_V7_COUNTS,
        "combined_count": EXPECTED_PARENT_COUNT,
        "ordered_indices_sha256": stable_hash(ordered_indices),
    }


def load_deletion_boundary(
    manifest_path: Path, verification_path: Path
) -> tuple[dict, dict]:
    manifest_path = manifest_path.resolve()
    verification_path = verification_path.resolve()
    manifest_hash = engine.sha256(manifest_path)
    verification_hash = engine.sha256(verification_path)
    require(
        manifest_hash == EXPECTED_INPUT_HASHES[manifest_path.name],
        "unexpected deletion-v2 manifest hash",
    )
    require(
        verification_hash == EXPECTED_INPUT_HASHES[verification_path.name],
        "unexpected deletion-v2 verification hash",
    )
    manifest = load_json(manifest_path)
    verification = load_json(verification_path)
    require(
        manifest.get("schema") == "d6-residue-18-deletion-manifest-v2"
        and manifest.get("summary", {}).get("parents") == 911
        and manifest.get("summary", {}).get("unique_deletions") == 11_975,
        "deletion-v2 manifest boundary",
    )
    require(
        verification.get("schema") == "d6-residue-18-deletion-verification-v2"
        and verification.get("status") == "PASS"
        and verification.get("manifest", {}).get("sha256") == manifest_hash
        and all(verification.get("checks", {}).values()),
        "deletion-v2 independent verification boundary",
    )
    require(
        manifest.get("semantics", {}).get("nonedges")
        == "unconstrained and may also have distance one",
        "deletion-v2 optional-nonedge semantics",
    )
    records = manifest.get("unique_deletions")
    require(isinstance(records, list) and len(records) == 11_975, "deletion records")
    seen_ids: set[str] = set()
    for record in records:
        class_id = record.get("class_id")
        require(isinstance(class_id, str) and class_id not in seen_ids, "deletion class id")
        seen_ids.add(class_id)
        adjacency = tuple(record.get("adjacency", ()))
        validate_adjacency(adjacency, N)
        require(
            sum(row.bit_count() for row in adjacency) // 2 == record.get("edges"),
            "deletion edge count",
        )
    return manifest, {
        "kind": "independently_verified_d6_residue_18_deletions_v2",
        "manifest": {"path": str(manifest_path), "sha256": manifest_hash},
        "verification": {
            "path": str(verification_path),
            "sha256": verification_hash,
            "status": "PASS",
        },
        "unique_deletions": len(records),
    }


def build_dense_deletion_cover(parents: Sequence[dict], deletions: dict) -> tuple[list[dict], dict]:
    """Return the frozen maximum-edge/greedy cover and exact audit metadata."""

    parent_keys = [(str(parent["class"]), int(parent["index"])) for parent in parents]
    require(len(parent_keys) == len(set(parent_keys)) == EXPECTED_PARENT_COUNT, "parent keys")
    parent_set = set(parent_keys)
    occurrences: dict[tuple[str, int], list[tuple[int, int, dict]]] = defaultdict(list)
    records = deletions["unique_deletions"]
    for class_index, record in enumerate(records):
        for occurrence in record.get("occurrences", ()):
            parent = (occurrence.get("parent_class"), occurrence.get("parent_index"))
            if parent in parent_set:
                occurrences[parent].append((class_index, int(record["edges"]), occurrence))
    require(set(occurrences) == parent_set, "deletion manifest omits a v7 parent")

    dense_parent_records = []
    class_coverage: dict[int, set[tuple[str, int]]] = defaultdict(set)
    for parent in sorted(parent_set, key=parent_sort_key):
        options = occurrences[parent]
        require(
            len(options) == 19
            and {int(option[2]["deleted_vertex"]) for option in options} == set(range(19)),
            f"rooted deletion partition for {parent}",
        )
        maximum_edges = max(option[1] for option in options)
        dense_options = [option for option in options if option[1] == maximum_edges]
        for class_index, _edges, _occurrence in dense_options:
            class_coverage[class_index].add(parent)
        grouped: dict[int, list[int]] = defaultdict(list)
        for class_index, _edges, occurrence in dense_options:
            grouped[class_index].append(int(occurrence["deleted_vertex"]))
        dense_parent_records.append(
            {
                "parent": tagged_parent(*parent),
                "maximum_deletion_edges": maximum_edges,
                "dense_occurrences": len(dense_options),
                "dense_classes": [
                    {
                        "class_index": class_index,
                        "class_id": records[class_index]["class_id"],
                        "deleted_vertices": sorted(grouped[class_index]),
                    }
                    for class_index in sorted(grouped)
                ],
            }
        )

    uncovered = set(parent_set)
    selected_indices: list[int] = []
    greedy_steps = []
    while uncovered:
        candidates = [
            class_index
            for class_index, covered in class_coverage.items()
            if covered & uncovered
        ]
        require(bool(candidates), "dense deletion classes do not cover every parent")
        chosen = max(
            candidates,
            key=lambda class_index: (
                len(class_coverage[class_index] & uncovered),
                records[class_index]["class_id"],
            ),
        )
        newly_covered = sorted(class_coverage[chosen] & uncovered, key=parent_sort_key)
        selected_indices.append(chosen)
        greedy_steps.append(
            {
                "step": len(greedy_steps),
                "class_index": chosen,
                "class_id": records[chosen]["class_id"],
                "gain": len(newly_covered),
                "newly_covered_parents": [tagged_parent(*parent) for parent in newly_covered],
            }
        )
        uncovered.difference_update(newly_covered)

    require(len(selected_indices) == len(set(selected_indices)), "repeated greedy class")
    selected_graphs = []
    selected_records = []
    for class_index in selected_indices:
        record = records[class_index]
        covered = sorted(class_coverage[class_index], key=parent_sort_key)
        active_classes = {parent[0] for parent in covered}
        require(len(active_classes) == 1, "selected deletion spans v7 strata")
        population = "K7" if active_classes == {"K7"} else "K6"
        known_positive = known_realizable_deletion(record)
        selected_graphs.append(
            {
                "index": class_index,
                "population": population,
                "adjacency": tuple(map(int, record["adjacency"])),
            }
        )
        selected_records.append(
            {
                "class_index": class_index,
                "class_id": record["class_id"],
                "canonical_graph6": record["canonical_graph6"],
                "population": population,
                "edges": int(record["edges"]),
                "clique_number": int(record["clique_number"]),
                "active_parent_count": len(covered),
                "active_parents": [tagged_parent(*parent) for parent in covered],
                "known_realizable_positive_control": known_positive,
                "standard18_compatible": bool(record.get("standard18_compatible")),
            }
        )

    gain_histogram = {
        str(gain): count
        for gain, count in sorted(Counter(step["gain"] for step in greedy_steps).items())
    }
    population_counts = dict(Counter(graph["population"] for graph in selected_graphs))
    known_positives = [
        record for record in selected_records if record["known_realizable_positive_control"]
    ]
    class_ids = [records[index]["class_id"] for index in selected_indices]
    require(len(selected_graphs) == EXPECTED_SELECTED_COUNT, "selected class count drift")
    require(gain_histogram == EXPECTED_GREEDY_GAIN_HISTOGRAM, "greedy gain drift")
    require(population_counts == EXPECTED_SELECTED_POPULATIONS, "population drift")
    require(
        len(known_positives) == 1
        and known_positives[0]["class_index"] == EXPECTED_KNOWN_POSITIVE_CLASS_INDEX
        and known_positives[0]["class_id"] == EXPECTED_KNOWN_POSITIVE_CLASS_ID,
        "known-positive selection drift",
    )
    if len(EXPECTED_SELECTED_CLASS_INDICES_SHA256) == 64:
        require(
            stable_hash(selected_indices) == EXPECTED_SELECTED_CLASS_INDICES_SHA256,
            "selected class-index root drift",
        )
        require(
            stable_hash(class_ids) == EXPECTED_SELECTED_CLASS_IDS_SHA256,
            "selected class-id root drift",
        )
        require(
            stable_hash(greedy_steps) == EXPECTED_GREEDY_STEPS_SHA256,
            "greedy-step root drift",
        )

    # Engine execution is sorted by immutable deletion-manifest ordinal; the
    # greedy chronology remains separately committed above.
    selected_graphs.sort(key=lambda graph: graph["index"])
    metadata = {
        "algorithm": {
            "dense_filter": (
                "retain every deletion occurrence attaining its parent's "
                "maximum required-edge count"
            ),
            "greedy_primary_key": "maximum number of currently uncovered parents",
            "greedy_tie_break": "lexicographically largest class_id",
            "candidate_nonedges": "unconstrained and may also be unit distance",
        },
        "parents": {
            "count": len(parent_keys),
            "class_counts": dict(Counter(parent[0] for parent in parent_keys)),
            "ordered_tagged": [tagged_parent(*parent) for parent in parent_keys],
            "ordered_tagged_sha256": stable_hash(
                [tagged_parent(*parent) for parent in parent_keys]
            ),
            "dense_records": dense_parent_records,
            "dense_records_sha256": stable_hash(dense_parent_records),
        },
        "eligible_dense_classes": len(class_coverage),
        "selected_count": len(selected_indices),
        "selected_population_counts": population_counts,
        "selected_class_indices_greedy_order": selected_indices,
        "selected_class_indices_sha256": stable_hash(selected_indices),
        "selected_class_ids_sha256": stable_hash(class_ids),
        "selected_records": selected_records,
        "selected_records_sha256": stable_hash(selected_records),
        "greedy_gain_histogram": gain_histogram,
        "greedy_steps": greedy_steps,
        "greedy_steps_sha256": stable_hash(greedy_steps),
        "known_positive_selected": known_positives,
        "coverage_complete": True,
    }
    return selected_graphs, metadata


def build_actionable_campaign(
    selected_graphs: Sequence[dict], cover: dict, deletions: dict
) -> tuple[list[dict], dict]:
    """Retain the greedy cover and add a non-positive path for blocked parents.

    The known standard-realizable selected class is deliberately retained as
    a full-production-parameter positive control.  For every parent whose
    selected coverers are all known positive, choose a dense non-positive
    backup maximizing its number of active v7 parents, with the same
    lexicographically-largest-class-ID tie break as the base greedy cover.
    """

    deletion_records = deletions["unique_deletions"]
    selected_records = cover["selected_records"]
    selected_indices = {int(graph["index"]) for graph in selected_graphs}
    require(
        selected_indices
        == {int(record["class_index"]) for record in selected_records},
        "base selected graph/metadata mismatch",
    )

    selected_coverers: dict[str, list[dict]] = defaultdict(list)
    for record in selected_records:
        for parent in record["active_parents"]:
            selected_coverers[parent].append(record)
    blocked_parents = sorted(
        (
            parent
            for parent, records in selected_coverers.items()
            if records
            and all(record["known_realizable_positive_control"] for record in records)
        ),
        key=lambda tagged: parent_sort_key(
            (tagged.split(":", 1)[0], int(tagged.split(":", 1)[1]))
        ),
    )

    dense_by_parent = {
        record["parent"]: record for record in cover["parents"]["dense_records"]
    }
    dense_coverage: dict[int, set[str]] = defaultdict(set)
    for parent_record in cover["parents"]["dense_records"]:
        for option in parent_record["dense_classes"]:
            dense_coverage[int(option["class_index"])].add(parent_record["parent"])

    backup_records = []
    backup_indices: set[int] = set()
    for parent in blocked_parents:
        candidates = []
        for option in dense_by_parent[parent]["dense_classes"]:
            class_index = int(option["class_index"])
            record = deletion_records[class_index]
            if known_realizable_deletion(record):
                continue
            active_parents = sorted(
                dense_coverage[class_index],
                key=lambda tagged: parent_sort_key(
                    (tagged.split(":", 1)[0], int(tagged.split(":", 1)[1]))
                ),
            )
            candidates.append((len(active_parents), record["class_id"], class_index))
        require(bool(candidates), f"no non-positive dense backup for {parent}")
        _gain, _class_id, chosen = max(candidates)
        record = deletion_records[chosen]
        active_parents = sorted(
            dense_coverage[chosen],
            key=lambda tagged: parent_sort_key(
                (tagged.split(":", 1)[0], int(tagged.split(":", 1)[1]))
            ),
        )
        active_classes = {tagged.split(":", 1)[0] for tagged in active_parents}
        require(len(active_classes) == 1, "backup deletion spans v7 strata")
        population = "K7" if active_classes == {"K7"} else "K6"
        backup_indices.add(chosen)
        backup_records.append(
            {
                "blocked_parent": parent,
                "class_index": chosen,
                "class_id": record["class_id"],
                "canonical_graph6": record["canonical_graph6"],
                "population": population,
                "edges": int(record["edges"]),
                "clique_number": int(record["clique_number"]),
                "active_parent_count": len(active_parents),
                "active_parents": active_parents,
                "known_realizable_positive_control": False,
                "standard18_compatible": False,
                "selection_key": [len(active_parents), record["class_id"]],
            }
        )

    require(len(backup_indices) == len(backup_records), "duplicate actionable backup")
    require(not selected_indices & backup_indices, "backup already in greedy cover")
    frozen_backups = [
        {
            "blocked_parent": record["blocked_parent"],
            "class_index": record["class_index"],
            "class_id": record["class_id"],
        }
        for record in backup_records
    ]
    require(frozen_backups == EXPECTED_ACTIONABLE_BACKUPS, "actionable backup drift")

    production_graphs = [graph.copy() for graph in selected_graphs]
    for record in backup_records:
        source = deletion_records[record["class_index"]]
        production_graphs.append(
            {
                "index": record["class_index"],
                "population": record["population"],
                "adjacency": tuple(map(int, source["adjacency"])),
            }
        )
    production_graphs.sort(key=lambda graph: graph["index"])
    production_indices = [int(graph["index"]) for graph in production_graphs]
    production_populations = dict(
        Counter(graph["population"] for graph in production_graphs)
    )
    require(len(production_graphs) == EXPECTED_PRODUCTION_COUNT, "production count drift")
    require(
        production_populations == EXPECTED_PRODUCTION_POPULATIONS,
        "production population drift",
    )
    if len(EXPECTED_PRODUCTION_CLASS_INDICES_SHA256) == 64:
        require(
            stable_hash(production_indices) == EXPECTED_PRODUCTION_CLASS_INDICES_SHA256,
            "production class-index root drift",
        )
    return production_graphs, {
        "algorithm": (
            "retain the known-positive base class; for each parent having only "
            "known-positive selected coverers, add a non-positive dense class "
            "maximizing active-parent count then lexicographically largest class_id"
        ),
        "blocked_parents": blocked_parents,
        "backup_count": len(backup_records),
        "backup_records": backup_records,
        "backup_records_sha256": stable_hash(backup_records),
        "production_count": len(production_graphs),
        "production_population_counts": production_populations,
        "production_class_indices": production_indices,
        "production_class_indices_sha256": stable_hash(production_indices),
        "all_blocked_parents_have_nonpositive_backups": True,
        "known_positive_retained_as_full_parameter_control": True,
    }


def select_run_graphs(
    graphs: Sequence[dict],
    *,
    shards: int,
    shard: int,
    sample: int | None,
    sample_seed: int,
    indices: frozenset[int] | None,
) -> list[dict]:
    selected = [graph.copy() for graph in graphs]
    if indices is not None:
        available = {int(graph["index"]) for graph in selected}
        missing = sorted(indices - available)
        require(not missing, f"requested deletion class indices absent: {missing}")
        selected = [graph for graph in selected if int(graph["index"]) in indices]
    selected = [
        graph for position, graph in enumerate(selected) if position % shards == shard
    ]
    if sample is not None and sample < len(selected):
        random.Random(sample_seed).shuffle(selected)
        selected = selected[:sample]
    selected.sort(key=lambda graph: graph["index"])
    for ordinal, graph in enumerate(selected):
        graph["ordinal"] = ordinal
    return selected


def profile_placement_orders(graphs: Sequence[dict], max_orders: int) -> dict:
    """Inventory the wrapper's exact order generator without calling the kernel.

    The selected canonical deletion of parent 3950509 is a genuine exception:
    the current placement heuristic produces no complete order for it.  Every
    generated order for every other selected deletion uses one or two circle
    placements.  Recording the exception prevents a silent no-call result.
    """

    require(max_orders > 0, "placement-order cap must be positive")
    records = []
    circle_histogram: Counter[int] = Counter()
    minimum_circle_histogram: Counter[str] = Counter()
    no_order = []
    for graph in graphs:
        orders = engine.cdriver6.gen_orders(graph["adjacency"], N, kmax=max_orders)
        circles = [
            engine.ncircle(graph["adjacency"], seed, order)
            for seed, order in orders
        ]
        circle_histogram.update(circles)
        if not orders:
            no_order.append(int(graph["index"]))
            minimum_circle_histogram["unreachable"] += 1
        else:
            minimum_circle_histogram[str(min(circles))] += 1
        records.append(
            {
                "class_index": int(graph["index"]),
                "orders_available": len(orders),
                "circle_count_histogram": {
                    str(key): value for key, value in sorted(Counter(circles).items())
                },
            }
        )
    require(
        all(circle_count in {1, 2} for circle_count in circle_histogram),
        "a selected deletion unexpectedly has a zero- or >2-circle order",
    )
    full_cover = len(graphs) == EXPECTED_PRODUCTION_COUNT and {
        int(graph["index"]) for graph in graphs
    } >= {EXPECTED_NO_ORDER_CLASS_INDEX, EXPECTED_KNOWN_POSITIVE_CLASS_INDEX}
    if full_cover:
        require(
            no_order == [EXPECTED_NO_ORDER_CLASS_INDEX],
            "full-cover no-order exception drift",
        )
        require(
            dict(sorted(minimum_circle_histogram.items()))
            == EXPECTED_MINIMUM_CIRCLE_HISTOGRAM,
            "full-cover minimum-circle profile drift",
        )
    return {
        "order_cap": max_orders,
        "graphs_profiled": len(graphs),
        "graphs_with_orders": len(graphs) - len(no_order),
        "no_order_class_indices": no_order,
        "expected_full_cover_exception": {
            "class_index": EXPECTED_NO_ORDER_CLASS_INDEX,
            "class_id": EXPECTED_NO_ORDER_CLASS_ID,
            "active_parent": "K6_only:3950509",
            "meaning": "UNRESOLVED without a kernel call under this order generator",
        },
        "generated_circle_count_histogram": {
            str(key): value for key, value in sorted(circle_histogram.items())
        },
        "minimum_circle_count_by_graph_histogram": dict(
            sorted(minimum_circle_histogram.items())
        ),
        "all_generated_orders_use_one_or_two_circles": True,
        "records": records,
        "records_sha256": stable_hash(records),
    }


@contextmanager
def engine_order_18() -> Iterator[None]:
    """Scope the legacy runner's module-level order without leaking it."""

    previous = engine.N
    engine.N = N
    try:
        yield
    finally:
        engine.N = previous


def verify_positive_input_boundaries() -> tuple[dict, dict]:
    for name, expected in EXPECTED_INPUT_HASHES.items():
        require(engine.sha256(ROOT / name) == expected, f"positive input hash: {name}")
    standard = load_json(ROOT / "d6_standard18_geometry_report.json")
    standard_check = load_json(ROOT / "d6_standard18_geometry_verification.json")
    exact = load_json(ROOT / "d6_18_exact_reconstruction.json")
    exact_check = load_json(ROOT / "d6_18_exact_reconstruction_verification.json")
    require(
        standard.get("schema") == "d6-standard18-exact-geometry-v1"
        and standard.get("status") == "PASS"
        and standard_check.get("status") == "PASS",
        "standard18 positive boundary",
    )
    require(
        exact.get("schema") == "d6-18-exact-nonstandard-reconstruction-v1"
        and exact.get("status") == "PASS"
        and exact_check.get("status") == "PASS"
        and len(exact.get("configurations", ())) == 2,
        "nonstandard exact18 positive boundary",
    )
    return standard, exact


def run_n18_kernel_controls(deletions: dict, slices: int) -> dict:
    """Exercise negative, flexible-positive, and three exact positive inputs."""

    standard, exact = verify_positive_input_boundaries()
    base = engine.run_kernel_controls(slices)
    base["positive"]["description"] = (
        "K7 plus eleven distinct generic points on the common-unit circle "
        "of five seed vertices"
    )
    control_inputs = [
        (
            "selected_standard18_compatible_deletion",
            EXPECTED_KNOWN_POSITIVE_CLASS_INDEX,
            deletions["unique_deletions"][EXPECTED_KNOWN_POSITIVE_CLASS_INDEX][
                "adjacency"
            ],
        ),
        ("standard18_full_unit_graph", -101, standard["unit_graph"]["adjacency"]),
    ]
    control_inputs.extend(
        (
            f"nonstandard_exact18_{position + 1}",
            -102 - position,
            configuration["exact_geometry"]["unit_graph_adjacency"],
        )
        for position, configuration in enumerate(exact["configurations"])
    )
    exact_controls = []
    for ordinal, (name, index, adjacency_values) in enumerate(control_inputs):
        adjacency = tuple(map(int, adjacency_values))
        validate_adjacency(adjacency, N)
        result = engine.analyze_graph(
            {
                "ordinal": ordinal,
                "index": index,
                "population": "POSITIVE_CONTROL",
                "adjacency": adjacency,
            },
            max_orders=4,
            max_nodes=1_000,
            zero_circle_only=False,
            include_bulk_order=True,
            slices=slices,
        )
        require(result.status != "KILLED", f"kernel falsely killed {name}")
        result_record = asdict(result)
        # Controls are part of the checkpoint configuration hash.  Wall time
        # is intentionally excluded so an otherwise identical restart has the
        # same campaign identity.
        result_record.pop("elapsed_seconds")
        exact_controls.append(
            {
                "name": name,
                "required": "overall status is not KILLED",
                "status": result.status,
                "result": result_record,
            }
        )
    return {
        "status": "PASS",
        "graph_order": N,
        "flexible_circle_and_K8": base,
        "known_realizable_exact18": exact_controls,
        "positive_source_hashes": dict(sorted(EXPECTED_INPUT_HASHES.items())),
    }


def committed_source_boundary() -> dict:
    """Return a resume-stable, fail-closed executable-source boundary.

    Unrelated commits and porcelain changes must not perturb the campaign
    hash.  The boundary therefore names the newest commit touching any of the
    five executable sources and verifies every working file against both that
    commit's tree and HEAD's tree.
    """

    source_names = (
        Path(__file__).name,
        "run_d6_interval_residue.py",
        "cdriver6.py",
        "ckernel6.c",
        "ival.py",
    )

    def git(*arguments: str, binary: bool = False):
        return subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=not binary,
        ).stdout

    branch = str(git("branch", "--show-current")).strip()
    require(branch == "codex/dimension6", "production branch is not codex/dimension6")
    source_commit = str(
        git("rev-list", "-1", "HEAD", "--", *source_names)
    ).strip()
    require(len(source_commit) == 40, "could not resolve executable source commit")
    sources = {}
    for name in source_names:
        working_hash = engine.sha256(ROOT / name)
        head_blob = git("show", f"HEAD:{name}", binary=True)
        source_blob = git("show", f"{source_commit}:{name}", binary=True)
        head_hash = hashlib.sha256(head_blob).hexdigest()
        source_hash = hashlib.sha256(source_blob).hexdigest()
        require(
            working_hash == head_hash == source_hash,
            f"uncommitted or mixed-boundary production source: {name}",
        )
        sources[name] = working_hash
    return {
        "branch": branch,
        "source_commit": source_commit,
        "sources": dict(sorted(sources.items())),
        "all_sources_equal_head_and_source_commit_blobs": True,
        "identity_excludes_unrelated_head_and_porcelain_state": True,
    }


def compiler_version() -> str:
    try:
        return subprocess.run(
            ["cc", "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        ).stdout.splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError) as error:
        return f"unavailable: {type(error).__name__}: {error}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "d6_current_residue_manifest_v7.json"
    )
    parser.add_argument(
        "--verification",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v7_verification.json",
    )
    parser.add_argument(
        "--deletions",
        type=Path,
        default=ROOT / "d6_residue_18_deletions_v2.json",
    )
    parser.add_argument(
        "--deletion-verification",
        type=Path,
        default=ROOT / "d6_residue_18_deletions_v2_verification.json",
    )
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--orders", type=int, default=4)
    parser.add_argument("--cap", type=int, default=100_000)
    parser.add_argument("--slices", type=int, default=24)
    parser.add_argument("--sample", type=int)
    parser.add_argument("--sample-seed", type=int, default=618_263_179)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument(
        "--index",
        dest="indices",
        action="append",
        type=int,
        help="restrict to this deletion-manifest class ordinal; repeatable",
    )
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=ROOT / ".runs/d6_interval_18_cover_v7_cap100000",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / ".runs/d6_interval_18_cover_v7_cap100000_report.json",
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=ROOT / ".runs/d6_interval_18_cover_v7_cap100000_decisions.tsv",
    )
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=1)
    parser.add_argument("--retry-infra-errors", action="store_true")
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--outer-launch-command")
    args = parser.parse_args()
    if (
        args.workers < 1
        or args.orders < 1
        or args.cap < 1
        or args.slices < 1
        or args.checkpoint_every < 1
        or args.progress_every < 1
    ):
        parser.error("worker/search/checkpoint values must be positive")
    if args.shards < 1 or not 0 <= args.shard < args.shards:
        parser.error("require shards >= 1 and 0 <= shard < shards")
    if args.sample is not None and args.sample < 1:
        parser.error("sample must be positive")

    parents, parent_provenance = load_v7_parents(args.manifest, args.verification)
    deletions, deletion_provenance = load_deletion_boundary(
        args.deletions, args.deletion_verification
    )
    base_graphs, cover = build_dense_deletion_cover(parents, deletions)
    all_graphs, augmentation = build_actionable_campaign(
        base_graphs, cover, deletions
    )
    placement_profile = profile_placement_orders(all_graphs, args.orders)
    graphs = select_run_graphs(
        all_graphs,
        shards=args.shards,
        shard=args.shard,
        sample=args.sample,
        sample_seed=args.sample_seed,
        indices=frozenset(args.indices) if args.indices else None,
    )
    run_selection = {
        "graphs": len(graphs),
        "class_indices": [graph["index"] for graph in graphs],
        "class_indices_sha256": stable_hash([graph["index"] for graph in graphs]),
        "population_counts": dict(Counter(graph["population"] for graph in graphs)),
        "sampling": {
            "sample": args.sample,
            "sample_seed": args.sample_seed,
            "shards": args.shards,
            "shard": args.shard,
            "requested_indices": sorted(set(args.indices or ())),
        },
    }
    selection_summary = {
        "parent_boundary": parent_provenance,
        "deletion_boundary": deletion_provenance,
        "full_dense_cover": cover,
        "actionable_augmentation": augmentation,
        "placement_order_profile": placement_profile,
        "run_selection": run_selection,
    }
    if args.selection_only:
        print(json.dumps(selection_summary, indent=2, sort_keys=True))
        return

    source_boundary = committed_source_boundary()
    with engine_order_18():
        with engine._KERNEL_LOCK:
            engine.cdriver6._kernel()  # type: ignore[attr-defined]
        kernel_binary = Path(engine.cdriver6.HERE) / engine.cdriver6.LIBNAME
        controls = run_n18_kernel_controls(deletions, args.slices)
        search = {
            "orders": args.orders,
            "cap": args.cap,
            "slices": args.slices,
            # Every eligible selected graph needs one or two circle
            # placements, so zero-circle filtering would discard all work.
            "zero_circle_only": False,
            # Include position zero: it is itself a circle order at n=18.
            "include_bulk_order": True,
        }
        configuration = {
            "schema": 1,
            "graph_order": N,
            "sources": {
                **source_boundary["sources"],
                engine.cdriver6.LIBNAME: engine.sha256(kernel_binary),
            },
            "inputs": {
                "v7_manifest": parent_provenance["manifest"],
                "v7_verification": parent_provenance["verification"],
                "deletion_manifest": deletion_provenance["manifest"],
                "deletion_verification": deletion_provenance["verification"],
            },
            "selection": selection_summary,
            "search": search,
            "trust_assumptions": TRUST_ASSUMPTIONS,
            "environment": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "cpu_count": os.cpu_count(),
                "compiler": compiler_version(),
                "kernel_compile_command": (
                    f"cc -O2 -shared -o {engine.cdriver6.LIBNAME} ckernel6.c -lm"
                ),
            },
            "git": source_boundary,
            "launch": {
                "python_argv": [sys.executable, *sys.argv],
                "python_command": shlex.join([sys.executable, *sys.argv]),
                "outer_launch_command": args.outer_launch_command,
                "outer_launch_command_user_supplied": args.outer_launch_command is not None,
            },
            "kernel_controls": controls,
        }
        print(
            f"interval n=18 dense cover: {len(graphs)} classes, "
            f"{args.workers} native threads, up to {args.orders} orders, "
            f"cap {args.cap}; checkpoint {args.checkpoint_dir}",
            file=sys.stderr,
            flush=True,
        )
        report, has_infra = engine.run_campaign(
            graphs,
            configuration=configuration,
            workers=args.workers,
            checkpoint_every=args.checkpoint_every,
            progress_every=args.progress_every,
            retry_infra_errors=args.retry_infra_errors,
            checkpoint_dir=args.checkpoint_dir.resolve(),
            report_path=args.output.resolve(),
            decisions_path=args.decisions.resolve(),
        )

    positive_run_indices = {
        record["class_index"]
        for record in cover["known_positive_selected"]
        if record["class_index"] in set(run_selection["class_indices"])
    }
    falsely_killed = [
        result["index"]
        for result in report["results"]
        if result["index"] in positive_run_indices and result["status"] == "KILLED"
    ]
    if falsely_killed:
        raise RuntimeError(f"production search killed known positive classes: {falsely_killed}")
    counts = report["status_counts"]
    print(
        "complete: "
        + "; ".join(f"{status} {counts[status]}" for status in engine.RESULT_STATUSES),
        file=sys.stderr,
        flush=True,
    )
    if has_infra:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
