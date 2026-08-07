#!/usr/bin/env python3
"""Build the canonical 18-vertex deletion corpus of the exact d=6 v4 residue.

Every 19-vertex graph in the v4 residue has 19 induced vertex deletions.  A
realization of the parent would realize every deletion, so a certified
non-realizable deletion is a reusable obstruction.  This builder performs no
geometric rejection: it reconstructs, canonically labels, deduplicates, and
classifies the deletion graphs, and recognizes the minimal unit supports
contained in the known halfcube-plus-two-apices realization.

Canonical labeling is operational metadata, not a mathematical rejection.
The exact parent/deletion incidences and the embedded adjacencies are retained
so that downstream checkers need not trust a canonical label to use a result.
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
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parent
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
    # Filled by the exact Clebsch-plus-apex-pair compatibility classifier.
    "standard_compatible_unique_deletions": 14,
    "standard_compatible_deletion_occurrences": 39,
    "standard_compatible_parent_graphs": 16,
}


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


def validate_rows(rows: object, vertices: int, label: str) -> tuple[int, ...]:
    if (
        not isinstance(rows, list)
        or len(rows) != vertices
        or any(type(row) is not int for row in rows)
    ):
        raise ValueError(f"{label} is not a {vertices}-row integer adjacency")
    bound = 1 << vertices
    result = tuple(rows)
    for vertex, row in enumerate(result):
        if row < 0 or row >= bound or row & (1 << vertex):
            raise ValueError(f"{label} has an invalid row {vertex}")
        for other in range(vertices):
            if bool(row & (1 << other)) != bool(result[other] & (1 << vertex)):
                raise ValueError(f"{label} is asymmetric at {vertex},{other}")
    return result


def edge_count(rows: Sequence[int]) -> int:
    degree_sum = sum(row.bit_count() for row in rows)
    if degree_sum & 1:
        raise ValueError("odd adjacency degree sum")
    return degree_sum // 2


def delete_vertex(rows: Sequence[int], deleted: int) -> tuple[int, ...]:
    n = len(rows)
    if not 0 <= deleted < n:
        raise ValueError("deleted vertex out of range")
    low_mask = (1 << deleted) - 1
    result = []
    for vertex, row in enumerate(rows):
        if vertex == deleted:
            continue
        low = row & low_mask
        high = row >> (deleted + 1)
        result.append(low | (high << deleted))
    return tuple(result)


def compress_mask(mask: int, deleted: int) -> int:
    low_mask = (1 << deleted) - 1
    return (mask & low_mask) | ((mask >> (deleted + 1)) << deleted)


def alpha_at_most_two(rows: Sequence[int]) -> bool:
    n = len(rows)
    full = (1 << n) - 1
    for first in range(n):
        nonneighbours = full & ~rows[first] & ~(1 << first)
        later = nonneighbours & ~((1 << (first + 1)) - 1)
        while later:
            bit = later & -later
            second = bit.bit_length() - 1
            third = nonneighbours & ~rows[second] & ~(1 << first) & ~(1 << second)
            if third:
                return False
            later ^= bit
    return True


def contains_clique(rows: Sequence[int], target: int) -> bool:
    def search(candidates: int, need: int) -> bool:
        if need == 0:
            return True
        if candidates.bit_count() < need:
            return False
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            if search(candidates & rows[vertex], need - 1):
                return True
            if candidates.bit_count() < need:
                break
        return False

    return search((1 << len(rows)) - 1, target)


def induced_rows(rows: Sequence[int], vertices: Sequence[int]) -> tuple[int, ...]:
    position = {old: new for new, old in enumerate(vertices)}
    result = [0] * len(vertices)
    for new, old in enumerate(vertices):
        neighbours = rows[old]
        for other_old, other_new in position.items():
            if neighbours & (1 << other_old):
                result[new] |= 1 << other_new
    return tuple(result)


def is_clebsch_graph(rows: Sequence[int]) -> bool:
    """Recognize the strongly regular Clebsch graph by exact parameters."""

    if len(rows) != 16 or any(row.bit_count() != 5 for row in rows):
        return False
    for first in range(16):
        for second in range(first):
            common = (rows[first] & rows[second]).bit_count()
            if rows[first] & (1 << second):
                if common != 0:
                    return False
            elif common != 2:
                return False
    return True


def graph_isomorphism(
    source: Sequence[int], target: Sequence[int]
) -> list[int] | None:
    """Return one exact source-to-target isomorphism by backtracking."""

    if len(source) != len(target):
        return None
    n = len(source)
    if sorted(row.bit_count() for row in source) != sorted(
        row.bit_count() for row in target
    ):
        return None
    mapping = [-1] * n
    used = 0

    def search(mapped_count: int) -> bool:
        nonlocal used
        if mapped_count == n:
            return True
        best = -1
        best_score = -1
        for vertex in range(n):
            if mapping[vertex] >= 0:
                continue
            score = sum(
                1 for other in range(n)
                if mapping[other] >= 0 and source[vertex] & (1 << other)
            )
            if score > best_score:
                best = vertex
                best_score = score
        degree = source[best].bit_count()
        for image in range(n):
            if used & (1 << image) or target[image].bit_count() != degree:
                continue
            compatible = True
            for other in range(n):
                other_image = mapping[other]
                if other_image < 0:
                    continue
                if bool(source[best] & (1 << other)) != bool(
                    target[image] & (1 << other_image)
                ):
                    compatible = False
                    break
            if not compatible:
                continue
            mapping[best] = image
            used |= 1 << image
            if search(mapped_count + 1):
                return True
            used ^= 1 << image
            mapping[best] = -1
        return False

    return mapping if search(0) else None


def standard_compatibility_witnesses(rows: Sequence[int]) -> list[dict]:
    """Return pole pairs proving E(rows) embeds in the full standard graph.

    The complement of the full standard unit graph is the disjoint union of a
    Clebsch graph and the edge between the two poles.  Since every deletion
    complement is triangle-free, a spanning required-edge embedding exists
    exactly when some candidate nonedge pair leaves an induced Clebsch graph.
    """

    if len(rows) != 18:
        raise ValueError("standard compatibility expects 18 vertices")
    witnesses: list[dict] = []
    target_unit = standard_eighteen_graph()
    complete = (1 << 16) - 1
    target_nonunit = tuple(
        complete & ~target_unit[vertex] & ~(1 << vertex) for vertex in range(16)
    )
    for first in range(18):
        for second in range(first):
            if rows[first] & (1 << second):
                continue
            base = [vertex for vertex in range(18) if vertex not in (second, first)]
            unit_base = induced_rows(rows, base)
            full = (1 << 16) - 1
            nonunit_base = tuple(
                full & ~unit_base[vertex] & ~(1 << vertex) for vertex in range(16)
            )
            if is_clebsch_graph(nonunit_base):
                base_isomorphism = graph_isomorphism(nonunit_base, target_nonunit)
                if base_isomorphism is None:
                    raise AssertionError("Clebsch parameter match lacked isomorphism")
                embedding = [-1] * 18
                for local, old in enumerate(base):
                    embedding[old] = base_isomorphism[local]
                embedding[second] = 16
                embedding[first] = 17
                for old in range(18):
                    for other in range(old):
                        if rows[old] & (1 << other):
                            if not target_unit[embedding[old]] & (1 << embedding[other]):
                                raise AssertionError("invalid standard embedding witness")
                witnesses.append(
                    {
                        "pole_pair": [second, first],
                        "embedding_permutation": embedding,
                    }
                )
    return witnesses


def graph6_encode(rows: Sequence[int]) -> str:
    n = len(rows)
    if not 0 <= n <= 62:
        raise ValueError("this compact graph6 encoder supports at most 62 vertices")
    values = []
    accumulator = 0
    used = 0
    for column in range(1, n):
        for row in range(column):
            accumulator = (accumulator << 1) | int(bool(rows[row] & (1 << column)))
            used += 1
            if used == 6:
                values.append(accumulator)
                accumulator = 0
                used = 0
    if used:
        values.append(accumulator << (6 - used))
    return chr(n + 63) + "".join(chr(value + 63) for value in values)


def graph6_decode(code: str) -> tuple[int, ...]:
    if not code or code.startswith(">>"):
        raise ValueError("expected a compact headerless graph6 record")
    n = ord(code[0]) - 63
    if not 0 <= n <= 62:
        raise ValueError("unsupported graph6 order")
    bits: list[int] = []
    for character in code[1:]:
        value = ord(character) - 63
        if not 0 <= value < 64:
            raise ValueError("invalid graph6 character")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    needed = n * (n - 1) // 2
    if len(bits) < needed:
        raise ValueError("truncated graph6 record")
    rows = [0] * n
    cursor = 0
    for column in range(1, n):
        for row in range(column):
            if bits[cursor]:
                rows[row] |= 1 << column
                rows[column] |= 1 << row
            cursor += 1
    result = tuple(rows)
    if graph6_encode(result) != code:
        raise ValueError("noncanonical graph6 padding or trailing data")
    return result


def canonicalize(codes: Sequence[str], labelg: Path) -> list[str]:
    if not codes:
        return []
    process = subprocess.run(
        [str(labelg), "-q", "-g"],
        input="\n".join(codes) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise RuntimeError(f"labelg failed: {process.stderr.strip()}")
    result = [line.strip() for line in process.stdout.splitlines() if line.strip()]
    if len(result) != len(codes):
        raise RuntimeError(
            f"labelg returned {len(result)} records for {len(codes)} inputs"
        )
    for code in result:
        graph6_decode(code)
    return result


def add_edge(rows: list[int], first: int, second: int) -> None:
    rows[first] |= 1 << second
    rows[second] |= 1 << first


def standard_eighteen_graph() -> tuple[int, ...]:
    words = [word for word in range(32) if word.bit_count() & 1]
    rows = [0] * 18
    for first, word in enumerate(words):
        for second, other in enumerate(words[:first]):
            if (word ^ other).bit_count() == 2:
                add_edge(rows, first, second)
    for apex in (16, 17):
        for base in range(16):
            add_edge(rows, apex, base)
    return tuple(rows)


def independent_masks(rows: Sequence[int]) -> list[int]:
    result = []
    for mask in range(1 << len(rows)):
        remaining = mask
        good = True
        while remaining:
            bit = remaining & -remaining
            vertex = bit.bit_length() - 1
            remaining ^= bit
            if rows[vertex] & remaining:
                good = False
                break
        if good:
            result.append(mask)
    return result


def standard_minimal_supports() -> list[tuple[int, ...]]:
    """Enumerate every labeled edge-minimal alpha<=2 support in the 18-set.

    The nonunit graph of the base is the Clebsch graph.  The full nonunit graph
    is Clebsch plus the edge between the two apices.  A triangle-free
    supergraph is obtained by assigning disjoint Clebsch-independent sets A,B
    to the two apex neighbourhoods.  It is maximal precisely when every
    unassigned base vertex has a Clebsch neighbour in both A and B.
    """

    full_unit = standard_eighteen_graph()
    base_full = (1 << 16) - 1
    clebsch = [0] * 16
    for vertex in range(16):
        clebsch[vertex] = base_full & ~full_unit[vertex] & ~(1 << vertex)
    independent = independent_masks(clebsch)
    supports: list[tuple[int, ...]] = []
    for first in independent:
        for second in independent:
            if first & second:
                continue
            outside = base_full & ~(first | second)
            remaining = outside
            maximal = True
            while remaining:
                bit = remaining & -remaining
                vertex = bit.bit_length() - 1
                remaining ^= bit
                if not (clebsch[vertex] & first) or not (clebsch[vertex] & second):
                    maximal = False
                    break
            if not maximal:
                continue
            rows = list(full_unit)
            for vertex in range(16):
                if first & (1 << vertex):
                    rows[16] &= ~(1 << vertex)
                    rows[vertex] &= ~(1 << 16)
                if second & (1 << vertex):
                    rows[17] &= ~(1 << vertex)
                    rows[vertex] &= ~(1 << 17)
            if not alpha_at_most_two(rows):
                raise AssertionError("constructed standard support has alpha at least three")
            supports.append(tuple(rows))
    return supports


def load_input(path: Path) -> tuple[dict, list[tuple[str, int, tuple[int, ...]]]]:
    if sha256(path) != EXPECTED_INPUT_SHA256:
        raise ValueError("the v4 residue manifest does not match its pinned hash")
    verification_path = path.with_name(INPUT_VERIFICATION_NAME)
    if sha256(verification_path) != EXPECTED_INPUT_VERIFICATION_SHA256:
        raise ValueError("the v4 independent verification hash does not match")
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    if verification.get("status") != "PASS":
        raise ValueError("the v4 independent verification did not pass")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "d6-current-exact-residue-v4":
        raise ValueError("unexpected v4 residue schema")
    parents = []
    for class_name in ("K7", "K6_only"):
        records = payload["classes"][class_name]["graphs"]
        for record in records:
            index = record.get("index")
            if type(index) is not int:
                raise ValueError("parent index is not an integer")
            rows = validate_rows(record.get("adjacency"), 19, f"parent {index}")
            parents.append((class_name, index, rows))
    if len(parents) != EXPECTED_COUNTS["parents"]:
        raise ValueError("wrong v4 parent count")
    return payload, parents


def _histogram(values: Iterable[int]) -> dict[str, int]:
    return {str(key): count for key, count in sorted(Counter(values).items())}


def build_manifest(input_path: Path, labelg: Path) -> dict:
    source, parents = load_input(input_path)
    labelg_hash = sha256(labelg)

    standard_supports = standard_minimal_supports()
    standard_inputs = [graph6_encode(rows) for rows in standard_supports]

    occurrences = []
    deletion_inputs = []
    for class_name, parent_index, rows in parents:
        for deleted in range(19):
            deletion = delete_vertex(rows, deleted)
            if not alpha_at_most_two(deletion):
                raise AssertionError("an induced deletion lost alpha<=2")
            occurrences.append(
                {
                    "parent_class": class_name,
                    "parent_index": parent_index,
                    "deleted_vertex": deleted,
                    "attachment_mask": compress_mask(rows[deleted], deleted),
                    "attachment_degree": rows[deleted].bit_count(),
                }
            )
            deletion_inputs.append(graph6_encode(deletion))

    all_canonical = canonicalize(standard_inputs + deletion_inputs, labelg)
    standard_canonical = all_canonical[: len(standard_inputs)]
    deletion_canonical = all_canonical[len(standard_inputs) :]

    standard_hist = Counter(standard_canonical)
    standard_types = []
    standard_code_to_id = {}
    for ordinal, code in enumerate(
        sorted(standard_hist, key=lambda item: (edge_count(graph6_decode(item)), item)),
        start=1,
    ):
        identifier = f"standard18-support-{ordinal}"
        standard_code_to_id[code] = identifier
        rows = graph6_decode(code)
        standard_types.append(
            {
                "id": identifier,
                "canonical_graph6": code,
                "adjacency": list(rows),
                "edges": edge_count(rows),
                "labeled_supports": standard_hist[code],
            }
        )

    grouped: dict[str, list[dict]] = defaultdict(list)
    for code, occurrence in zip(deletion_canonical, occurrences):
        grouped[code].append(occurrence)

    unique_records = []
    class_codes: dict[str, set[str]] = {"K7": set(), "K6_only": set()}
    standard_occurrences = []
    standard_parents = set()
    compatible_occurrences = []
    compatible_parents = set()
    compatible_unique = 0
    for code in sorted(grouped):
        rows = graph6_decode(code)
        group = grouped[code]
        classes = sorted({item["parent_class"] for item in group})
        for class_name in classes:
            class_codes[class_name].add(code)
        support_id = standard_code_to_id.get(code)
        if support_id is not None:
            standard_occurrences.extend(group)
            standard_parents.update(
                (item["parent_class"], item["parent_index"]) for item in group
            )
        compatibility = standard_compatibility_witnesses(rows)
        if compatibility:
            compatible_unique += 1
            compatible_occurrences.extend(group)
            compatible_parents.update(
                (item["parent_class"], item["parent_index"]) for item in group
            )
        unique_records.append(
            {
                "class_id": "u18-" + hashlib.sha256(code.encode("ascii")).hexdigest(),
                "canonical_graph6": code,
                "adjacency": list(rows),
                "edges": edge_count(rows),
                "degree_sequence": sorted(
                    (row.bit_count() for row in rows), reverse=True
                ),
                "clique_number": 7 if contains_clique(rows, 7) else 6,
                "contains_K7": contains_clique(rows, 7),
                "contains_K6": contains_clique(rows, 6),
                "parent_classes": classes,
                "occurrence_count": len(group),
                "parent_count": len(
                    {(item["parent_class"], item["parent_index"]) for item in group}
                ),
                "occurrences": group,
                "standard18_support_type": support_id,
                "standard18_compatible": bool(compatibility),
                "standard18_pole_pair_witnesses": compatibility,
            }
        )

    observed_counts = {
        "parents": len(parents),
        "deletion_occurrences": len(occurrences),
        "unique_deletions": len(unique_records),
        "unique_K7_parent_deletions": len(class_codes["K7"]),
        "unique_K6_only_parent_deletions": len(class_codes["K6_only"]),
        "cross_class_unique_overlap": len(
            class_codes["K7"] & class_codes["K6_only"]
        ),
        "standard_support_types": len(standard_types),
        "standard_support_labeled": len(standard_supports),
        "standard_support_deletion_occurrences": len(standard_occurrences),
        "standard_support_parent_graphs": len(standard_parents),
        "standard_compatible_unique_deletions": compatible_unique,
        "standard_compatible_deletion_occurrences": len(compatible_occurrences),
        "standard_compatible_parent_graphs": len(compatible_parents),
    }
    if observed_counts != EXPECTED_COUNTS:
        raise AssertionError(
            f"deletion corpus count drift: {observed_counts!r} != {EXPECTED_COUNTS!r}"
        )

    parent_indices = [index for _, index, _ in parents]
    unique_codes = [record["canonical_graph6"] for record in unique_records]
    summary = {
        **observed_counts,
        "parent_class_counts": {
            "K7": sum(class_name == "K7" for class_name, _, _ in parents),
            "K6_only": sum(
                class_name == "K6_only" for class_name, _, _ in parents
            ),
        },
        "occurrence_edge_histogram": _histogram(
            edge_count(graph6_decode(code)) for code in deletion_canonical
        ),
        "unique_edge_histogram": _histogram(
            record["edges"] for record in unique_records
        ),
        "unique_clique_classes": {
            "contains_K7": sum(record["contains_K7"] for record in unique_records),
            "K6_without_K7": sum(
                record["contains_K6"] and not record["contains_K7"]
                for record in unique_records
            ),
            "no_K6": sum(not record["contains_K6"] for record in unique_records),
        },
        "parent_indices_sha256": stable_hash(parent_indices),
        "ordered_occurrences_sha256": stable_hash(occurrences),
        "ordered_deletion_canonical_graph6_sha256": stable_hash(deletion_canonical),
        "unique_canonical_graph6_sha256": stable_hash(unique_codes),
        "unique_records_sha256": stable_hash(unique_records),
        "standard_support_types_sha256": stable_hash(standard_types),
        "standard_support_occurrences_sha256": stable_hash(standard_occurrences),
        "standard_compatible_occurrences_sha256": stable_hash(
            compatible_occurrences
        ),
    }

    return {
        "schema": "d6-residue-18-deletion-manifest-v1",
        "source": {
            "path": input_path.name,
            "sha256": EXPECTED_INPUT_SHA256,
            "schema": source["schema"],
            "combined_count": source["combined"]["count"],
            "ordered_indices_sha256": source["combined"][
                "ordered_indices_sha256"
            ],
        },
        "source_verification": {
            "path": INPUT_VERIFICATION_NAME,
            "sha256": EXPECTED_INPUT_VERIFICATION_SHA256,
            "status": "PASS",
        },
        "canonicalizer": {
            "program": "labelg",
            "path_at_build": str(labelg),
            "binary_sha256_at_build": labelg_hash,
            "command": ["labelg", "-q", "-g"],
            "role": (
                "deduplication and recognition metadata only; no geometric "
                "rejection follows from a canonical label"
            ),
        },
        "summary": summary,
        "standard18": {
            "full_unit_graph": list(standard_eighteen_graph()),
            "full_unit_graph_sha256": stable_hash(
                list(standard_eighteen_graph())
            ),
            "minimal_support_types": standard_types,
            "matched_deletion_occurrences": standard_occurrences,
            "matched_parent_records": [
                {"parent_class": class_name, "parent_index": parent_index}
                for class_name, parent_index in sorted(standard_parents)
            ],
            "compatible_deletion_occurrences": compatible_occurrences,
            "compatible_parent_records": [
                {"parent_class": class_name, "parent_index": parent_index}
                for class_name, parent_index in sorted(compatible_parents)
            ],
            "semantics": (
                "standard18_compatible is a spanning non-induced required-edge "
                "embedding into the full standard unit graph. It says that the "
                "known coordinates realize the deletion support. It neither "
                "proves uniqueness of that realization nor rejects the "
                "19-vertex parent."
            ),
        },
        "unique_deletions": unique_records,
        "semantics": {
            "edges": "required unit distances",
            "nonedges": "unconstrained and may also have distance one",
            "deletion_obstruction_rule": (
                "if a unique deletion graph is independently certified "
                "non-realizable by distinct points in R6 using only its required "
                "edges, every listed parent occurrence is non-realizable"
            ),
            "current_rejections": 0,
            "realizability_claims": 0,
            "canonical_isomorphism_claim_only": True,
        },
        "nonclaims": [
            "No deletion graph is rejected or claimed realizable by this manifest.",
            "A standard18 support match does not prove global uniqueness of its embedding.",
            "Failure to extend the standard coordinates does not reject another embedding.",
            "Canonical labeling is not a geometric certificate.",
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
    parser.add_argument("--input", type=Path, default=ROOT / INPUT_NAME)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_residue_18_deletions.json"
    )
    parser.add_argument("--labelg", type=Path)
    args = parser.parse_args()
    labelg = args.labelg or Path(shutil.which("labelg") or "")
    if not labelg or not labelg.is_file():
        raise SystemExit("labelg was not found; pass --labelg /absolute/path")
    manifest = build_manifest(args.input.resolve(), labelg.resolve())
    atomic_json(args.output, manifest)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output),
                **{
                    key: manifest["summary"][key]
                    for key in (
                        "parents",
                        "deletion_occurrences",
                        "unique_deletions",
                        "standard_support_deletion_occurrences",
                    )
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
