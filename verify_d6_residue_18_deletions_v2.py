#!/usr/bin/env python3
"""Independent checker for the dimension-six v5 18-deletion corpus.

The checker imports neither the v2 builder nor a geometric rejection kernel.
It reconstructs all 17,309 rooted deletions from the v5 parent adjacencies,
canonically labels them in a fresh batch, independently rebuilds the seven
standard18 support types and all compatibility pole pairs, and compares every
embedded unique record, occurrence, attachment mask, witness, count, and hash.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import shutil
import time
from collections import Counter, defaultdict
from pathlib import Path

import verify_d6_residue_18_deletions as reference


ROOT = Path(__file__).resolve().parent
BUILDER = "build_d6_residue_18_deletions_v2.py"
EXPECTED_BUILDER_SHA256 = (
    "ecf3b388af8e7216d099e1ad2b5eaaaf4151c6d145bc7e47fceef6bd2745c374"
)
INPUT_NAME = "d6_current_residue_manifest_v5.json"
INPUT_VERIFICATION_NAME = "d6_current_residue_manifest_v5_verification.json"
EXPECTED_FILES = {
    INPUT_NAME: "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5",
    INPUT_VERIFICATION_NAME: "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46",
    "verify_d6_residue_18_deletions.py": "5ac2baab623a49de603fe0d5d30581402e15a1753d48a02b8ac2ce3c4ec251d4",
}
EXPECTED_COUNTS = {
    "parents": 911,
    "deletion_occurrences": 17_309,
    "unique_deletions": 11_975,
    "unique_K7_parent_deletions": 2_288,
    "unique_K6_only_parent_deletions": 9_749,
    "cross_class_unique_overlap": 62,
    "standard_support_types": 7,
    "standard_support_labeled": 1_960,
    "standard_support_deletion_occurrences": 26,
    "standard_support_parent_graphs": 14,
    "standard_compatible_unique_deletions": 14,
    "standard_compatible_deletion_occurrences": 39,
    "standard_compatible_parent_graphs": 16,
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def object_sha256(value: object) -> str:
    data = json.dumps(
        value, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def histogram(values) -> dict[str, int]:
    return {str(key): count for key, count in sorted(Counter(values).items())}


def verify(manifest_path: Path, labelg: Path, expected_hash: str | None) -> dict:
    checks: dict[str, bool] = {}
    manifest_hash = file_sha256(manifest_path)
    if expected_hash is not None:
        require(manifest_hash == expected_hash, "manifest SHA-256 mismatch")
    require(file_sha256(ROOT / BUILDER) == EXPECTED_BUILDER_SHA256, "builder hash")
    for name, expected in EXPECTED_FILES.items():
        require(file_sha256(ROOT / name) == expected, f"pinned file changed: {name}")
    checks["source_hashes"] = True

    input_path = ROOT / INPUT_NAME
    input_check_path = ROOT / INPUT_VERIFICATION_NAME
    source = json.loads(input_path.read_text(encoding="utf-8"))
    source_check = json.loads(input_check_path.read_text(encoding="utf-8"))
    require(source.get("schema") == "d6-current-exact-residue-v5", "input schema")
    require(source.get("combined", {}).get("count") == 911, "input count")
    require(source_check.get("status") == "PASS", "input verification status")
    require(
        source_check.get("manifest", {}).get("sha256") == EXPECTED_FILES[INPUT_NAME],
        "input verification binding",
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(
        manifest.get("schema") == "d6-residue-18-deletion-manifest-v2",
        "manifest schema mismatch",
    )

    parents = []
    for class_name in ("K7", "K6_only"):
        for record in source["classes"][class_name]["graphs"]:
            index = record.get("index")
            require(type(index) is int, "noninteger parent index")
            parents.append(
                (
                    class_name,
                    index,
                    reference.rows_from_object(
                        record.get("adjacency"), 19, f"parent {index}"
                    ),
                )
            )
    require(len(parents) == EXPECTED_COUNTS["parents"], "parent count")

    occurrences = []
    packed_deletions = []
    for class_name, index, rows in parents:
        for deleted in range(19):
            deletion = reference.remove_vertex(rows, deleted)
            require(
                reference.no_independent_triple(deletion),
                "deletion has independent triple",
            )
            occurrences.append(
                {
                    "parent_class": class_name,
                    "parent_index": index,
                    "deleted_vertex": deleted,
                    "attachment_mask": reference.remove_bit(rows[deleted], deleted),
                    "attachment_degree": rows[deleted].bit_count(),
                }
            )
            packed_deletions.append(reference.graph6_pack(deletion))

    supports = reference.standard_support_records()
    packed_supports = [reference.graph6_pack(rows) for rows in supports]
    canonical = reference.canonical_batch(
        packed_deletions + packed_supports, labelg
    )
    deletion_codes = canonical[: len(packed_deletions)]
    support_codes = canonical[len(packed_deletions) :]
    expected_groups: dict[str, list[dict]] = defaultdict(list)
    for code, occurrence in zip(deletion_codes, occurrences):
        expected_groups[code].append(occurrence)
    support_histogram = Counter(support_codes)

    embedded = manifest.get("unique_deletions")
    require(isinstance(embedded, list), "unique_deletions is not a list")
    require(
        [record.get("canonical_graph6") for record in embedded]
        == sorted(expected_groups),
        "unique canonical order mismatch",
    )
    class_codes = {"K7": set(), "K6_only": set()}
    matched_occurrences = []
    matched_parents = set()
    compatible_occurrences = []
    compatible_parents = set()
    compatible_unique = 0
    support_code_set = set(support_histogram)
    standard_unit = reference.known_full_graph()
    for item in embedded:
        code = item.get("canonical_graph6")
        rows = reference.graph6_unpack(code)
        require(
            reference.rows_from_object(
                item.get("adjacency"), 18, "embedded deletion"
            )
            == rows,
            "embedded canonical adjacency mismatch",
        )
        require(
            item.get("class_id")
            == "u18-" + hashlib.sha256(code.encode("ascii")).hexdigest(),
            "stable class id",
        )
        group = expected_groups[code]
        require(item.get("occurrences") == group, "rooted occurrence mismatch")
        require(item.get("occurrence_count") == len(group), "occurrence count")
        parent_set = {
            (entry["parent_class"], entry["parent_index"]) for entry in group
        }
        require(item.get("parent_count") == len(parent_set), "parent count")
        classes = sorted({entry["parent_class"] for entry in group})
        require(item.get("parent_classes") == classes, "parent classes")
        for class_name in classes:
            class_codes[class_name].add(code)
        require(item.get("edges") == reference.number_of_edges(rows), "edge count")
        require(
            item.get("degree_sequence")
            == sorted((row.bit_count() for row in rows), reverse=True),
            "degree sequence",
        )
        has_k7 = reference.clique_exists(rows, 7)
        has_k6 = reference.clique_exists(rows, 6)
        require(item.get("contains_K7") is has_k7, "K7 flag")
        require(item.get("contains_K6") is has_k6, "K6 flag")
        require(item.get("clique_number") == (7 if has_k7 else 6), "clique number")

        if code in support_code_set:
            matched_occurrences.extend(group)
            matched_parents.update(parent_set)
            require(item.get("standard18_support_type") is not None, "support id")
        else:
            require(item.get("standard18_support_type") is None, "false support id")

        pole_pairs = reference.compatible_pole_pairs(rows)
        witnesses = item.get("standard18_pole_pair_witnesses")
        require(isinstance(witnesses, list), "standard witness shape")
        require(
            item.get("standard18_compatible") is bool(pole_pairs),
            "standard compatibility flag",
        )
        require(
            [witness.get("pole_pair") for witness in witnesses] == pole_pairs,
            "standard pole-pair completeness",
        )
        for witness in witnesses:
            permutation = witness.get("embedding_permutation")
            require(
                isinstance(permutation, list)
                and sorted(permutation) == list(range(18)),
                "standard embedding is not a permutation",
            )
            for first in range(18):
                for second in range(first):
                    if rows[first] & (1 << second):
                        require(
                            bool(
                                standard_unit[permutation[first]]
                                & (1 << permutation[second])
                            ),
                            "required edge missing in standard embedding",
                        )
        if pole_pairs:
            compatible_unique += 1
            compatible_occurrences.extend(group)
            compatible_parents.update(parent_set)
    checks["all_rooted_parent_deletions"] = True
    checks["canonical_unique_records"] = True

    standard = manifest.get("standard18")
    require(isinstance(standard, dict), "missing standard18 section")
    require(
        reference.rows_from_object(
            standard.get("full_unit_graph"), 18, "standard full graph"
        )
        == standard_unit,
        "standard full graph mismatch",
    )
    expected_types = []
    for ordinal, code in enumerate(
        sorted(
            support_histogram,
            key=lambda value: (
                reference.number_of_edges(reference.graph6_unpack(value)),
                value,
            ),
        ),
        start=1,
    ):
        rows = reference.graph6_unpack(code)
        expected_types.append(
            {
                "id": f"standard18-support-{ordinal}",
                "canonical_graph6": code,
                "adjacency": list(rows),
                "edges": reference.number_of_edges(rows),
                "labeled_supports": support_histogram[code],
            }
        )
    types = standard.get("minimal_support_types")
    require(types == expected_types, "standard support type mismatch")
    require(
        standard.get("matched_deletion_occurrences") == matched_occurrences,
        "standard matched occurrences",
    )
    require(
        standard.get("matched_parent_records")
        == [
            {"parent_class": class_name, "parent_index": index}
            for class_name, index in sorted(matched_parents)
        ],
        "standard matched parents",
    )
    require(
        standard.get("compatible_deletion_occurrences") == compatible_occurrences,
        "standard compatible occurrences",
    )
    require(
        standard.get("compatible_parent_records")
        == [
            {"parent_class": class_name, "parent_index": index}
            for class_name, index in sorted(compatible_parents)
        ],
        "standard compatible parents",
    )
    checks["standard18_reconstruction_and_embeddings"] = True

    observed = {
        "parents": len(parents),
        "deletion_occurrences": len(occurrences),
        "unique_deletions": len(expected_groups),
        "unique_K7_parent_deletions": len(class_codes["K7"]),
        "unique_K6_only_parent_deletions": len(class_codes["K6_only"]),
        "cross_class_unique_overlap": len(
            class_codes["K7"] & class_codes["K6_only"]
        ),
        "standard_support_types": len(support_histogram),
        "standard_support_labeled": len(supports),
        "standard_support_deletion_occurrences": len(matched_occurrences),
        "standard_support_parent_graphs": len(matched_parents),
        "standard_compatible_unique_deletions": compatible_unique,
        "standard_compatible_deletion_occurrences": len(compatible_occurrences),
        "standard_compatible_parent_graphs": len(compatible_parents),
    }
    require(observed == EXPECTED_COUNTS, "top-level exact count mismatch")
    summary = manifest.get("summary")
    require(isinstance(summary, dict), "summary is not an object")
    for key, value in observed.items():
        require(summary.get(key) == value, f"summary {key}")
    unique_codes = sorted(expected_groups)
    require(
        summary.get("parent_indices_sha256")
        == object_sha256([index for _, index, _ in parents]),
        "parent index hash",
    )
    require(
        summary.get("ordered_occurrences_sha256") == object_sha256(occurrences),
        "occurrence hash",
    )
    require(
        summary.get("ordered_deletion_canonical_graph6_sha256")
        == object_sha256(deletion_codes),
        "canonical occurrence hash",
    )
    require(
        summary.get("unique_canonical_graph6_sha256")
        == object_sha256(unique_codes),
        "unique code hash",
    )
    require(
        summary.get("unique_records_sha256") == object_sha256(embedded),
        "unique record hash",
    )
    require(
        summary.get("standard_support_types_sha256") == object_sha256(types),
        "support type hash",
    )
    require(
        summary.get("standard_support_occurrences_sha256")
        == object_sha256(matched_occurrences),
        "support occurrence hash",
    )
    require(
        summary.get("standard_compatible_occurrences_sha256")
        == object_sha256(compatible_occurrences),
        "compatible occurrence hash",
    )
    require(
        summary.get("occurrence_edge_histogram")
        == histogram(
            reference.number_of_edges(reference.graph6_unpack(code))
            for code in deletion_codes
        ),
        "occurrence edge histogram",
    )
    require(
        summary.get("unique_edge_histogram")
        == histogram(
            reference.number_of_edges(reference.graph6_unpack(code))
            for code in unique_codes
        ),
        "unique edge histogram",
    )
    require(
        summary.get("unique_clique_classes")
        == {
            "contains_K7": sum(
                reference.clique_exists(reference.graph6_unpack(code), 7)
                for code in unique_codes
            ),
            "K6_without_K7": sum(
                reference.clique_exists(reference.graph6_unpack(code), 6)
                and not reference.clique_exists(reference.graph6_unpack(code), 7)
                for code in unique_codes
            ),
            "no_K6": sum(
                not reference.clique_exists(reference.graph6_unpack(code), 6)
                for code in unique_codes
            ),
        },
        "unique clique histogram",
    )
    checks["summary_hashes_and_counts"] = True

    semantics = manifest.get("semantics")
    require(isinstance(semantics, dict), "semantics missing")
    require(semantics.get("current_rejections") == 0, "manifest claims rejection")
    require(
        semantics.get("realizability_claims") == 0,
        "manifest claims realization",
    )
    checks["no_geometric_overclaim"] = True

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require("build_d6_residue_18_deletions_v2" not in imports, "builder import")
    checks["import_independence"] = True

    return {
        "schema": "d6-residue-18-deletion-verification-v2",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "input": {"path": INPUT_NAME, "sha256": EXPECTED_FILES[INPUT_NAME]},
        "input_verification": {
            "path": INPUT_VERIFICATION_NAME,
            "sha256": EXPECTED_FILES[INPUT_VERIFICATION_NAME],
            "status": "PASS",
        },
        "builder_sha256": EXPECTED_BUILDER_SHA256,
        "verifier_source_sha256": file_sha256(Path(__file__).resolve()),
        "independent_reference_sha256": EXPECTED_FILES[
            "verify_d6_residue_18_deletions.py"
        ],
        "labelg": {"path": str(labelg), "sha256": file_sha256(labelg)},
        "counts": observed,
        "checks": checks,
        "trust": (
            "labelg is used only for isomorphism grouping; every parent "
            "deletion, rooted attachment, compatibility pair, and embedding "
            "is independently reconstructed"
        ),
    }


def atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "d6_residue_18_deletions_v2.json"
    )
    parser.add_argument("--expected-manifest-sha256")
    parser.add_argument("--labelg", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_residue_18_deletions_v2_verification.json",
    )
    args = parser.parse_args()
    executable = args.labelg or Path(shutil.which("labelg") or "")
    if not executable or not executable.is_file():
        raise SystemExit("labelg was not found")
    result = verify(
        args.manifest.resolve(),
        executable.resolve(),
        args.expected_manifest_sha256,
    )
    atomic_write(args.output.resolve(), result)
    print(
        json.dumps(
            {**result["counts"], "status": "PASS", "output": str(args.output)},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
