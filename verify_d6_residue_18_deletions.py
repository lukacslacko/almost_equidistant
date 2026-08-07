#!/usr/bin/env python3
"""Independent checker for the dimension-six v4 18-deletion manifest.

The checker imports neither the builder nor a geometric production kernel. It
reconstructs all 18,240 induced deletions from the hash-pinned v4 residue,
canonically labels them in one fresh invocation, verifies all 12,712 embedded
unique graphs and parent incidences, and independently reconstructs the seven
minimal support types carried by the known 18-point realization.

This verifies a corpus and its combinatorial metadata. It deliberately gives
no non-realizability credit to a canonical label or a standard-coordinate
support match.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
BUILDER = "build_d6_residue_18_deletions.py"
EXPECTED_BUILDER_SHA256 = (
    "f44ef42f8da3e6ba4b4c902b2db11468d781b8d9157eab3081a9df43674a08fd"
)
INPUT_NAME = "d6_current_residue_manifest_v4.json"
INPUT_VERIFICATION_NAME = "d6_current_residue_manifest_v4_verification.json"
EXPECTED_INPUT_SHA256 = (
    "6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4"
)
EXPECTED_INPUT_VERIFICATION_SHA256 = (
    "765eceb6782d131dae8c77c9b94de78735e23a070da051c8fffbd984790e1a41"
)
EXPECTED_COUNTS = {
    "parents": 960,
    "deletion_occurrences": 18_240,
    "unique_deletions": 12_712,
    "unique_K7_parent_deletions": 2_288,
    "unique_K6_only_parent_deletions": 10_486,
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


def rows_from_object(value: object, order: int, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == order
        and all(type(item) is int for item in value),
        f"{label}: malformed adjacency",
    )
    rows = tuple(value)
    for first in range(order):
        require(0 <= rows[first] < (1 << order), f"{label}: row out of range")
        require(not rows[first] & (1 << first), f"{label}: loop")
        for second in range(first):
            require(
                bool(rows[first] & (1 << second))
                == bool(rows[second] & (1 << first)),
                f"{label}: asymmetry",
            )
    return rows


def remove_vertex(rows: Sequence[int], vertex: int) -> tuple[int, ...]:
    answer = []
    below = (1 << vertex) - 1
    for old, neighbourhood in enumerate(rows):
        if old == vertex:
            continue
        answer.append(
            (neighbourhood & below) | ((neighbourhood >> (vertex + 1)) << vertex)
        )
    return tuple(answer)


def remove_bit(mask: int, vertex: int) -> int:
    below = (1 << vertex) - 1
    return (mask & below) | ((mask >> (vertex + 1)) << vertex)


def number_of_edges(rows: Sequence[int]) -> int:
    total = sum(row.bit_count() for row in rows)
    require(total % 2 == 0, "odd degree sum")
    return total // 2


def graph6_pack(rows: Sequence[int]) -> str:
    order = len(rows)
    require(order <= 62, "large graph6 order is unsupported")
    sequence = []
    for high in range(1, order):
        for low in range(high):
            sequence.append(int(bool(rows[low] & (1 << high))))
    while len(sequence) % 6:
        sequence.append(0)
    encoded = [chr(order + 63)]
    for start in range(0, len(sequence), 6):
        value = 0
        for bit in sequence[start : start + 6]:
            value = 2 * value + bit
        encoded.append(chr(value + 63))
    return "".join(encoded)


def graph6_unpack(code: str) -> tuple[int, ...]:
    require(bool(code) and not code.startswith(">>"), "bad graph6 header")
    order = ord(code[0]) - 63
    require(0 <= order <= 62, "bad graph6 order")
    bitstream = []
    for character in code[1:]:
        value = ord(character) - 63
        require(0 <= value < 64, "bad graph6 byte")
        bitstream.extend((value >> bit) & 1 for bit in (5, 4, 3, 2, 1, 0))
    required = order * (order - 1) // 2
    require(len(bitstream) >= required, "short graph6 payload")
    rows = [0] * order
    cursor = 0
    for high in range(1, order):
        for low in range(high):
            if bitstream[cursor]:
                rows[low] |= 1 << high
                rows[high] |= 1 << low
            cursor += 1
    answer = tuple(rows)
    require(graph6_pack(answer) == code, "graph6 record has trailing data")
    return answer


def canonical_batch(records: Sequence[str], executable: Path) -> list[str]:
    if not records:
        return []
    run = subprocess.run(
        [str(executable), "-q", "-g"],
        input="\n".join(records) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    require(run.returncode == 0, f"labelg failed: {run.stderr}")
    result = [line.strip() for line in run.stdout.splitlines() if line.strip()]
    require(len(result) == len(records), "labelg changed the record count")
    for record in result:
        graph6_unpack(record)
    return result


def no_independent_triple(rows: Sequence[int]) -> bool:
    order = len(rows)
    for first in range(order):
        for second in range(first):
            if rows[first] & (1 << second):
                continue
            for third in range(second):
                if not rows[first] & (1 << third) and not rows[second] & (1 << third):
                    return False
    return True


def clique_exists(rows: Sequence[int], size: int) -> bool:
    def recurse(candidates: list[int], chosen: int) -> bool:
        if chosen == size:
            return True
        if chosen + len(candidates) < size:
            return False
        while candidates:
            vertex = candidates.pop()
            following = [other for other in candidates if rows[vertex] & (1 << other)]
            if recurse(following, chosen + 1):
                return True
        return False

    return recurse(list(range(len(rows))), 0)


def restricted_graph(rows: Sequence[int], keep: Sequence[int]) -> tuple[int, ...]:
    answer = [0] * len(keep)
    for new_first, old_first in enumerate(keep):
        for new_second, old_second in enumerate(keep[:new_first]):
            if rows[old_first] & (1 << old_second):
                answer[new_first] |= 1 << new_second
                answer[new_second] |= 1 << new_first
    return tuple(answer)


def clebsch_parameters(rows: Sequence[int]) -> bool:
    if len(rows) != 16 or sorted(row.bit_count() for row in rows) != [5] * 16:
        return False
    for first in range(16):
        for second in range(first):
            common = (rows[first] & rows[second]).bit_count()
            expected = 0 if rows[first] & (1 << second) else 2
            if common != expected:
                return False
    return True


def compatible_pole_pairs(unit_rows: Sequence[int]) -> list[list[int]]:
    witnesses = []
    for high in range(18):
        for low in range(high):
            if unit_rows[high] & (1 << low):
                continue
            keep = [vertex for vertex in range(18) if vertex not in (low, high)]
            base_unit = restricted_graph(unit_rows, keep)
            complete = (1 << 16) - 1
            base_nonunit = tuple(
                complete & ~base_unit[vertex] & ~(1 << vertex)
                for vertex in range(16)
            )
            if clebsch_parameters(base_nonunit):
                witnesses.append([low, high])
    return witnesses


def known_full_graph() -> tuple[int, ...]:
    words = [word for word in range(32) if word.bit_count() % 2 == 1]
    rows = [0] * 18

    def join(first: int, second: int) -> None:
        rows[first] |= 1 << second
        rows[second] |= 1 << first

    for first in range(16):
        for second in range(first):
            if (words[first] ^ words[second]).bit_count() == 2:
                join(first, second)
    for pole in (16, 17):
        for base in range(16):
            join(pole, base)
    return tuple(rows)


def is_independent(rows: Sequence[int], mask: int) -> bool:
    chosen = [vertex for vertex in range(len(rows)) if mask & (1 << vertex)]
    return all(
        not rows[first] & (1 << second)
        for position, first in enumerate(chosen)
        for second in chosen[:position]
    )


def standard_support_records() -> list[tuple[int, ...]]:
    """Independent reconstruction using generic complement maximality checks."""

    unit = known_full_graph()
    base_mask = (1 << 16) - 1
    clebsch = tuple(
        base_mask & ~unit[vertex] & ~(1 << vertex) for vertex in range(16)
    )
    stable_sets = [
        mask for mask in range(1 << 16) if is_independent(clebsch, mask)
    ]
    answers = []
    for left in stable_sets:
        for right in stable_sets:
            if left & right:
                continue
            complement = list(clebsch) + [1 << 17, 1 << 16]
            for vertex in range(16):
                if left & (1 << vertex):
                    complement[16] |= 1 << vertex
                    complement[vertex] |= 1 << 16
                if right & (1 << vertex):
                    complement[17] |= 1 << vertex
                    complement[vertex] |= 1 << 17

            # Test maximal triangle-freeness directly. The starting graph is
            # triangle-free by the independent/disjoint construction.
            triangle_free = True
            for first in range(18):
                for second in range(first):
                    if complement[first] & (1 << second) and (
                        complement[first] & complement[second]
                    ):
                        triangle_free = False
                        break
                if not triangle_free:
                    break
            if not triangle_free:
                continue
            maximal = True
            for first in range(18):
                for second in range(first):
                    if complement[first] & (1 << second):
                        continue
                    if not (complement[first] & complement[second]):
                        maximal = False
                        break
                if not maximal:
                    break
            if not maximal:
                continue

            full = (1 << 18) - 1
            support = tuple(
                full & ~complement[vertex] & ~(1 << vertex)
                for vertex in range(18)
            )
            require(no_independent_triple(support), "bad reconstructed support")
            answers.append(support)
    return answers


def verify(manifest_path: Path, labelg: Path, expected_hash: str | None) -> dict:
    checks: dict[str, bool] = {}
    manifest_hash = file_sha256(manifest_path)
    if expected_hash is not None:
        require(manifest_hash == expected_hash, "manifest SHA-256 mismatch")
    input_path = ROOT / INPUT_NAME
    require(file_sha256(input_path) == EXPECTED_INPUT_SHA256, "input hash mismatch")
    input_verification_path = ROOT / INPUT_VERIFICATION_NAME
    require(
        file_sha256(input_verification_path)
        == EXPECTED_INPUT_VERIFICATION_SHA256,
        "input verification hash mismatch",
    )
    input_verification = json.loads(
        input_verification_path.read_text(encoding="utf-8")
    )
    require(input_verification.get("status") == "PASS", "input verification status")
    require(
        file_sha256(ROOT / BUILDER) == EXPECTED_BUILDER_SHA256,
        "builder source hash mismatch",
    )
    checks["source_hashes"] = True

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(
        manifest.get("schema") == "d6-residue-18-deletion-manifest-v1",
        "manifest schema mismatch",
    )
    source = json.loads(input_path.read_text(encoding="utf-8"))
    require(source.get("schema") == "d6-current-exact-residue-v4", "input schema")

    parents = []
    for class_name in ("K7", "K6_only"):
        for record in source["classes"][class_name]["graphs"]:
            index = record.get("index")
            require(type(index) is int, "noninteger parent index")
            parents.append(
                (
                    class_name,
                    index,
                    rows_from_object(record.get("adjacency"), 19, f"parent {index}"),
                )
            )
    require(len(parents) == EXPECTED_COUNTS["parents"], "parent count")

    occurrences = []
    packed_deletions = []
    for class_name, index, rows in parents:
        for deleted in range(19):
            deletion = remove_vertex(rows, deleted)
            require(no_independent_triple(deletion), "deletion has independent triple")
            occurrences.append(
                {
                    "parent_class": class_name,
                    "parent_index": index,
                    "deleted_vertex": deleted,
                    "attachment_mask": remove_bit(rows[deleted], deleted),
                    "attachment_degree": rows[deleted].bit_count(),
                }
            )
            packed_deletions.append(graph6_pack(deletion))

    supports = standard_support_records()
    packed_supports = [graph6_pack(rows) for rows in supports]
    canonical = canonical_batch(packed_deletions + packed_supports, labelg)
    deletion_codes = canonical[: len(packed_deletions)]
    support_codes = canonical[len(packed_deletions) :]

    expected_groups: dict[str, list[dict]] = defaultdict(list)
    for code, occurrence in zip(deletion_codes, occurrences):
        expected_groups[code].append(occurrence)
    support_histogram = Counter(support_codes)

    embedded = manifest.get("unique_deletions")
    require(isinstance(embedded, list), "unique_deletions is not a list")
    require(
        [item.get("canonical_graph6") for item in embedded]
        == sorted(expected_groups),
        "unique canonical order mismatch",
    )
    seen_occurrences = []
    class_codes = {"K7": set(), "K6_only": set()}
    matched_occurrences = []
    matched_parents = set()
    compatible_occurrences = []
    compatible_parents = set()
    compatible_unique = 0
    support_code_set = set(support_histogram)
    for item in embedded:
        code = item.get("canonical_graph6")
        rows = graph6_unpack(code)
        require(
            rows_from_object(item.get("adjacency"), 18, "embedded deletion") == rows,
            "embedded canonical adjacency mismatch",
        )
        require(
            item.get("class_id")
            == "u18-" + hashlib.sha256(code.encode("ascii")).hexdigest(),
            "stable class id",
        )
        expected_occurrences = expected_groups[code]
        require(item.get("occurrences") == expected_occurrences, "occurrence mismatch")
        require(item.get("occurrence_count") == len(expected_occurrences), "occurrence count")
        require(
            item.get("parent_count")
            == len(
                {
                    (entry["parent_class"], entry["parent_index"])
                    for entry in expected_occurrences
                }
            ),
            "parent count",
        )
        classes = sorted({entry["parent_class"] for entry in expected_occurrences})
        require(item.get("parent_classes") == classes, "parent classes")
        require(item.get("edges") == number_of_edges(rows), "edge count")
        require(
            item.get("degree_sequence")
            == sorted((row.bit_count() for row in rows), reverse=True),
            "degree sequence",
        )
        require(
            item.get("clique_number") == (7 if clique_exists(rows, 7) else 6),
            "clique number",
        )
        require(item.get("contains_K7") is clique_exists(rows, 7), "K7 flag")
        require(item.get("contains_K6") is clique_exists(rows, 6), "K6 flag")
        for class_name in classes:
            class_codes[class_name].add(code)
        seen_occurrences.extend(expected_occurrences)
        if code in support_code_set:
            matched_occurrences.extend(expected_occurrences)
            matched_parents.update(
                (entry["parent_class"], entry["parent_index"])
                for entry in expected_occurrences
            )
            require(item.get("standard18_support_type") is not None, "missing support id")
        else:
            require(item.get("standard18_support_type") is None, "false support id")
        pole_pairs = compatible_pole_pairs(rows)
        stored_witnesses = item.get("standard18_pole_pair_witnesses")
        require(isinstance(stored_witnesses, list), "standard witnesses shape")
        require(
            item.get("standard18_compatible") is bool(pole_pairs),
            "standard compatibility flag",
        )
        require(
            [witness.get("pole_pair") for witness in stored_witnesses] == pole_pairs,
            "standard pole-pair completeness",
        )
        standard_unit = known_full_graph()
        for witness in stored_witnesses:
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
            compatible_occurrences.extend(expected_occurrences)
            compatible_parents.update(
                (entry["parent_class"], entry["parent_index"])
                for entry in expected_occurrences
            )

    require(len(seen_occurrences) == len(occurrences), "lost occurrence")
    checks["all_parent_deletions"] = True
    checks["canonical_unique_records"] = True

    standard = manifest.get("standard18")
    require(isinstance(standard, dict), "missing standard18 section")
    require(
        rows_from_object(standard.get("full_unit_graph"), 18, "standard full")
        == known_full_graph(),
        "standard full graph mismatch",
    )
    types = standard.get("minimal_support_types")
    require(isinstance(types, list), "support types is not a list")
    expected_type_rows = []
    for ordinal, code in enumerate(
        sorted(support_histogram, key=lambda value: (number_of_edges(graph6_unpack(value)), value)),
        start=1,
    ):
        rows = graph6_unpack(code)
        expected_type_rows.append(
            {
                "id": f"standard18-support-{ordinal}",
                "canonical_graph6": code,
                "adjacency": list(rows),
                "edges": number_of_edges(rows),
                "labeled_supports": support_histogram[code],
            }
        )
    require(types == expected_type_rows, "standard support type mismatch")
    require(
        standard.get("matched_deletion_occurrences") == matched_occurrences,
        "standard occurrence mismatch",
    )
    require(
        standard.get("matched_parent_records")
        == [
            {"parent_class": class_name, "parent_index": index}
            for class_name, index in sorted(matched_parents)
        ],
        "standard parent mismatch",
    )
    require(
        standard.get("compatible_deletion_occurrences") == compatible_occurrences,
        "standard compatible occurrence mismatch",
    )
    require(
        standard.get("compatible_parent_records")
        == [
            {"parent_class": class_name, "parent_index": index}
            for class_name, index in sorted(compatible_parents)
        ],
        "standard compatible parent mismatch",
    )
    checks["standard_support_reconstruction"] = True

    observed = {
        "parents": len(parents),
        "deletion_occurrences": len(occurrences),
        "unique_deletions": len(expected_groups),
        "unique_K7_parent_deletions": len(class_codes["K7"]),
        "unique_K6_only_parent_deletions": len(class_codes["K6_only"]),
        "cross_class_unique_overlap": len(class_codes["K7"] & class_codes["K6_only"]),
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
    require(summary.get("ordered_occurrences_sha256") == object_sha256(occurrences), "occurrence hash")
    require(summary.get("ordered_deletion_canonical_graph6_sha256") == object_sha256(deletion_codes), "canonical hash")
    require(summary.get("unique_records_sha256") == object_sha256(embedded), "unique record hash")
    require(summary.get("standard_support_types_sha256") == object_sha256(types), "support type hash")
    require(
        summary.get("standard_compatible_occurrences_sha256")
        == object_sha256(compatible_occurrences),
        "standard compatible occurrence hash",
    )
    checks["summary_hashes_and_counts"] = True

    semantics = manifest.get("semantics")
    require(isinstance(semantics, dict), "semantics missing")
    require(semantics.get("current_rejections") == 0, "manifest claims a rejection")
    require(semantics.get("realizability_claims") == 0, "manifest claims a realization")
    checks["no_geometric_overclaim"] = True

    return {
        "schema": "d6-residue-18-deletion-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "input": {"path": INPUT_NAME, "sha256": EXPECTED_INPUT_SHA256},
        "input_verification": {
            "path": INPUT_VERIFICATION_NAME,
            "sha256": EXPECTED_INPUT_VERIFICATION_SHA256,
            "status": "PASS",
        },
        "builder_sha256": EXPECTED_BUILDER_SHA256,
        "labelg": {"path": str(labelg), "sha256": file_sha256(labelg)},
        "counts": observed,
        "checks": checks,
        "trust": (
            "labelg is used only to reproduce isomorphism grouping; every "
            "parent/deletion adjacency and occurrence is independently rebuilt"
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
        "--manifest", type=Path, default=ROOT / "d6_residue_18_deletions.json"
    )
    parser.add_argument("--expected-manifest-sha256")
    parser.add_argument("--labelg", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_residue_18_deletions_verification.json",
    )
    args = parser.parse_args()
    executable = args.labelg or Path(shutil.which("labelg") or "")
    if not executable or not executable.is_file():
        raise SystemExit("labelg was not found")
    result = verify(
        args.manifest.resolve(), executable.resolve(), args.expected_manifest_sha256
    )
    atomic_write(args.output, result)
    print(json.dumps({**result["counts"], "status": "PASS", "output": str(args.output)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
