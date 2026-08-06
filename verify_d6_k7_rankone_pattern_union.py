#!/usr/bin/env python3
"""Independent checker for ``d6_k7_rankone_pattern_union.json``.

This checker imports neither the union builder nor either production search.
It separately parses the frozen selection, both verified reports, and every
tetrad TSV decision, then reconstructs the overlap and cover-profile sections.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
FILES = {
    "selection": ROOT / "d6_k7_rankone_tetrad_full_selection.json",
    "tetrad_report": ROOT / "d6_k7_rankone_tetrad_full_report.json",
    "tetrad_verification": (
        ROOT / "d6_k7_rankone_tetrad_full_verification_report.json"
    ),
    "tetrad_decisions": (
        ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"
    ),
    "pattern_report": ROOT / "d6_n14_pattern_954_containment_report.json",
    "pattern_verification": (
        ROOT / "d6_n14_pattern_954_containment_verification.json"
    ),
    "builder": ROOT / "build_d6_k7_rankone_pattern_union.py",
    "verifier": Path(__file__).resolve(),
}
EXPECTED_SELECTION_FILE = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
EXPECTED_SELECTION_LIST = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)
EXPECTED_PATTERN_REPORT = (
    "a2d4fb5a07dcea580315730b4155fc0c08ecd0439a9111e733d29e2da0b5be62"
)
EXPECTED_PATTERN_VERIFICATION = (
    "a3a353e4d6cab462697f9565b69845f4a8f296fd6f8085ed77d68897257eab32"
)
EXPECTED_GRAPHS = 12_839
COUNT_FIELDS = [
    "seeds",
    "covers",
    "enhanced_passing_covers",
    "strict_h_passing_covers",
    "prior_dual_failed_covers",
    "prior_dual_passing_covers",
    "prior_passing_no_saturating_clique_covers",
    "covers_with_near_clique",
    "covers_without_near_clique",
    "near_cliques_tested",
    "tetrad_failed_covers",
    "tetrad_passing_covers",
]
DECISION_FIELDS = [
    "ordinal",
    "index",
    "status",
    *COUNT_FIELDS,
    "prior_dual_rejected",
    "tetrad_rejected",
    "marginal_tetrad_rejected",
    "error_type",
]


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


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if type(value) is not dict:
        raise ValueError(f"{path.name} is not a JSON object")
    return value


def int_list(value: object, label: str) -> list[int]:
    if type(value) is not list or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    if len(value) != len(set(value)):
        raise ValueError(f"{label} has duplicate entries")
    return list(value)


def ordered(universe: Sequence[int], members: set[int]) -> list[int]:
    return [index for index in universe if index in members]


def parse_boolean(value: str, label: str) -> bool:
    if value not in {"True", "False"}:
        raise ValueError(f"invalid {label} at decision row")
    return value == "True"


def derive_profile(row: dict) -> dict:
    prior = row["prior_dual_passing_covers"]
    no_sat = row["prior_passing_no_saturating_clique_covers"]
    sat = prior - no_sat
    with_near = row["covers_with_near_clique"]
    without_near = row["covers_without_near_clique"]
    killed = row["tetrad_failed_covers"]
    resistant = with_near - killed
    if min(sat, no_sat, with_near, without_near, resistant) < 0:
        raise ValueError(f"negative derived count at index {row['index']}")
    if with_near + without_near != no_sat:
        raise ValueError(f"near-clique counts do not partition at {row['index']}")
    if row["tetrad_passing_covers"] != sat + without_near + resistant:
        raise ValueError(f"passing-cover arithmetic fails at {row['index']}")
    if no_sat == 0:
        structure = "saturating_only"
    elif sat == 0:
        structure = "no_saturating_only"
    else:
        structure = "mixed"
    return {
        "structure": structure,
        "prior_passing_covers": prior,
        "saturating_covers": sat,
        "no_saturating_covers": no_sat,
        "no_saturating_with_near_clique": with_near,
        "no_saturating_without_near_clique": without_near,
        "tetrad_killed_covers": killed,
        "tetrad_resistant_near_clique_covers": resistant,
        "tetrad_passing_covers": row["tetrad_passing_covers"],
        "has_saturating_cover": sat > 0,
        "has_no_near_clique_cover": without_near > 0,
        "has_tetrad_resistant_near_clique_cover": resistant > 0,
    }


def parse_decisions(selection: Sequence[int], report: dict) -> tuple[list[dict], dict]:
    path = FILES["tetrad_decisions"]
    artifacts = report["artifacts"]
    if artifacts["decisions_archive"] != path.name:
        raise ValueError("decision archive path disagrees with report")
    compressed_hash = sha256(path)
    if compressed_hash != artifacts["decisions_archive_sha256"]:
        raise ValueError("decision archive compressed hash mismatch")
    with gzip.open(path, "rb") as stream:
        raw = stream.read()
    raw_hash = hashlib.sha256(raw).hexdigest()
    if raw_hash != artifacts["decisions_uncompressed_sha256"]:
        raise ValueError("decision archive raw hash mismatch")
    reader = csv.DictReader(
        io.StringIO(raw.decode("utf-8"), newline=""), delimiter="\t"
    )
    if reader.fieldnames != DECISION_FIELDS:
        raise ValueError("decision archive header mismatch")
    source_rows = list(reader)
    if len(source_rows) != len(selection):
        raise ValueError("decision archive row count mismatch")
    rows = []
    statuses = Counter()
    totals = Counter()
    for ordinal, (source, expected) in enumerate(zip(source_rows, selection)):
        if int(source["ordinal"]) != ordinal or int(source["index"]) != expected:
            raise ValueError(f"decision archive order mismatch at {ordinal}")
        if source["status"] not in {"REJECTED", "SURVIVOR"}:
            raise ValueError(f"bad decision status at {expected}")
        if source["error_type"]:
            raise ValueError(f"decision infrastructure error at {expected}")
        counts = {name: int(source[name]) for name in COUNT_FIELDS}
        if any(value < 0 for value in counts.values()):
            raise ValueError(f"negative source count at {expected}")
        prior = parse_boolean(source["prior_dual_rejected"], "prior flag")
        rejected = parse_boolean(source["tetrad_rejected"], "tetrad flag")
        marginal = parse_boolean(
            source["marginal_tetrad_rejected"], "marginal flag"
        )
        if prior or rejected != marginal:
            raise ValueError(f"decision is not a marginal tetrad result at {expected}")
        if (source["status"] == "REJECTED") != marginal:
            raise ValueError(f"decision status/flag mismatch at {expected}")
        row = {
            "index": expected,
            "status": source["status"],
            **counts,
        }
        row["profile"] = derive_profile(row)
        rows.append(row)
        statuses[source["status"]] += 1
        totals.update(counts)
    summary = report["summary"]
    if dict(statuses) != summary["status_counts"]:
        raise ValueError("decision status counts disagree with report")
    if {name: totals[name] for name in COUNT_FIELDS} != summary["totals"]:
        raise ValueError("decision field totals disagree with report")
    return rows, {
        "rows": len(rows),
        "status_counts": dict(statuses),
        "uncompressed_sha256": raw_hash,
    }


def set_payload(indices: list[int]) -> dict:
    return {
        "count": len(indices),
        "indices": indices,
        "indices_sha256": stable_hash(indices),
    }


def distribution(rows: Sequence[dict], members: set[int]) -> dict:
    chosen = [row for row in rows if row["index"] in members]
    structures = Counter(row["profile"]["structure"] for row in chosen)
    return {
        "graphs": len(chosen),
        "structure": {
            name: structures.get(name, 0)
            for name in ("saturating_only", "mixed", "no_saturating_only")
        },
        "obstacle_flags": {
            flag: sum(row["profile"][flag] for row in chosen)
            for flag in (
                "has_saturating_cover",
                "has_no_near_clique_cover",
                "has_tetrad_resistant_near_clique_cover",
            )
        },
    }


def verify(manifest_path: Path, expected_manifest_sha256: str) -> dict:
    if sha256(manifest_path) != expected_manifest_sha256:
        raise ValueError("union manifest hash mismatch")
    manifest = load(manifest_path)
    if (
        manifest.get("schema") != "d6-k7-rankone-pattern-union-v1"
        or manifest.get("status") != "COMPLETE"
    ):
        raise ValueError("union manifest schema/status mismatch")
    source_hashes = manifest.get("source_hashes", {})
    fixed_hashes = {
        FILES["selection"].name: EXPECTED_SELECTION_FILE,
        FILES["pattern_report"].name: EXPECTED_PATTERN_REPORT,
        FILES["pattern_verification"].name: EXPECTED_PATTERN_VERIFICATION,
    }
    for name, expected in fixed_hashes.items():
        if source_hashes.get(name) != expected:
            raise ValueError(f"manifest has wrong pinned hash for {name}")
    for path in FILES.values():
        if source_hashes.get(path.name) != sha256(path):
            raise ValueError(f"manifest/source hash mismatch for {path.name}")
    if set(source_hashes) != {path.name for path in FILES.values()}:
        raise ValueError("manifest source-hash key set mismatch")

    selection_document = load(FILES["selection"])
    selection = int_list(selection_document["selected_indices"], "selection")
    if (
        selection_document.get("schema") != 1
        or selection_document.get("kind")
        != "d6_k7_rankone_tetrad_full_selection"
        or len(selection) != EXPECTED_GRAPHS
        or stable_hash(selection) != EXPECTED_SELECTION_LIST
        or selection_document["selected_indices_sha256"]
        != EXPECTED_SELECTION_LIST
    ):
        raise ValueError("frozen selection mismatch")
    universe = set(selection)

    tetrad_report = load(FILES["tetrad_report"])
    tetrad_verification = load(FILES["tetrad_verification"])
    report_hash = sha256(FILES["tetrad_report"])
    verification_hash = sha256(FILES["tetrad_verification"])
    summary = tetrad_report["summary"]
    tetrad_rejected_list = int_list(
        summary["marginal_tetrad_rejected_indices"], "tetrad rejections"
    )
    tetrad_rejected = set(tetrad_rejected_list)
    tconfiguration = tetrad_report["configuration"]
    tselection = tconfiguration["selection"]
    if (
        tconfiguration["limit"] is not None
        or tselection["selection_report_sha256"] != EXPECTED_SELECTION_FILE
        or tselection["selected_indices_sha256"] != EXPECTED_SELECTION_LIST
        or summary["graphs"] != len(selection)
        or summary["complete"] != len(selection)
        or summary["infra_errors"] != 0
        or summary["marginal_tetrad_rejected"] != len(tetrad_rejected)
        or tetrad_rejected_list != ordered(selection, tetrad_rejected)
        or not tetrad_rejected <= universe
    ):
        raise ValueError("tetrad production report gate failed")
    provenance = tetrad_verification["provenance"]
    vsummary = tetrad_verification["summary"]
    if (
        tetrad_verification["errors"] != []
        or provenance["tetrad_report_sha256"] != report_hash
        or provenance["selection_sha256"] != EXPECTED_SELECTION_FILE
        or provenance["tetrad_decisions_sha256"]
        != tetrad_report["artifacts"]["decisions_archive_sha256"]
        or vsummary["graphs"] != len(selection)
        or vsummary["status_counts"] != {"PASS": len(selection)}
        or vsummary["verified_rejections"] != len(tetrad_rejected)
    ):
        raise ValueError("independent tetrad verification gate failed")
    rows, decision_audit = parse_decisions(selection, tetrad_report)
    if [row["index"] for row in rows if row["status"] == "REJECTED"] != (
        tetrad_rejected_list
    ):
        raise ValueError("decision archive/report rejection mismatch")

    pattern_report = load(FILES["pattern_report"])
    pattern_verification = load(FILES["pattern_verification"])
    pattern_list = int_list(
        pattern_report["certified_hit_indices"], "pattern hits"
    )
    pattern = set(pattern_list)
    counts = pattern_report["counts"]
    pattern_residue = ordered(selection, universe - pattern)
    if (
        pattern_report["status"] != "COMPLETE"
        or pattern_report.get("required_edge_only") is not True
        or pattern_report["processed"] != len(selection)
        or counts != {
            "HIT": len(pattern),
            "INFRA_ERROR": 0,
            "NO_HIT": len(selection) - len(pattern),
            "TIMEOUT": 0,
        }
        or pattern_list != ordered(selection, pattern)
        or pattern_report["certified_hit_indices_sha256"]
        != stable_hash(pattern_list)
        or pattern_report["unresolved_indices"] != []
    ):
        raise ValueError("pattern production report gate failed")
    if (
        pattern_verification["status"] != "PASS"
        or pattern_verification.get("required_edge_only") is not True
        or int_list(
            pattern_verification["certified_hit_indices"],
            "verified pattern hits",
        )
        != pattern_list
        or int_list(
            pattern_verification["residue_indices"],
            "verified pattern residue",
        )
        != pattern_residue
        or pattern_verification["mapping_verification"]["all_injective"]
        is not True
        or pattern_verification["mapping_verification"]["mappings_checked"]
        != len(pattern)
        or pattern_verification["selection"]["indices_sha256"]
        != EXPECTED_SELECTION_LIST
        or pattern_verification["report"]["sha256"]
        != EXPECTED_PATTERN_REPORT
    ):
        raise ValueError("independent pattern verification gate failed")

    both = tetrad_rejected & pattern
    tetrad_only = tetrad_rejected - pattern
    pattern_only = pattern - tetrad_rejected
    union = tetrad_rejected | pattern
    residue = universe - union
    tetrad_survivors = universe - tetrad_rejected
    computed_sets = {
        "tetrad_rejected": ordered(selection, tetrad_rejected),
        "pattern_954_rejected": ordered(selection, pattern),
        "both": ordered(selection, both),
        "tetrad_only": ordered(selection, tetrad_only),
        "pattern_954_only": ordered(selection, pattern_only),
        "exact_union": ordered(selection, union),
        "exact_residue": ordered(selection, residue),
    }
    expected_sets = {
        name: set_payload(indices) for name, indices in computed_sets.items()
    }
    if manifest["sets"] != expected_sets:
        raise ValueError("union manifest set reconstruction mismatch")
    expected_overlap = {
        "tetrad_rejected": len(tetrad_rejected),
        "pattern_954_rejected": len(pattern),
        "both": len(both),
        "tetrad_only": len(tetrad_only),
        "pattern_954_only": len(pattern_only),
        "exact_union": len(union),
        "exact_residue": len(residue),
    }
    if manifest["overlap"] != expected_overlap:
        raise ValueError("union manifest overlap counts mismatch")

    profiles = [
        {"index": row["index"], **row["profile"]}
        for row in rows
        if row["index"] in residue
    ]
    expected_structure = {
        "all_selected": distribution(rows, universe),
        "tetrad_survivors": distribution(rows, tetrad_survivors),
        "exact_union_residue": distribution(rows, residue),
        "pattern_hits_among_tetrad_survivors": distribution(rows, pattern_only),
        "residue_profiles": profiles,
    }
    if manifest["cover_structure"] != expected_structure:
        raise ValueError("union manifest cover structure mismatch")
    expected_gates = {
        "selection": {
            "graphs": len(selection),
            "indices_sha256": stable_hash(selection),
        },
        "tetrad": {
            "report_sha256": report_hash,
            "verification_sha256": verification_hash,
            "verification_status_counts": vsummary["status_counts"],
            "verified_rejections": vsummary["verified_rejections"],
            "decision_archive_audit": decision_audit,
        },
        "pattern_954": {
            "report_sha256": EXPECTED_PATTERN_REPORT,
            "verification_sha256": EXPECTED_PATTERN_VERIFICATION,
            "mappings_rechecked": pattern_verification[
                "mapping_verification"
            ]["mappings_checked"],
            "all_mappings_injective": pattern_verification[
                "mapping_verification"
            ]["all_injective"],
        },
    }
    if manifest["gates"] != expected_gates:
        raise ValueError("manifest campaign gates mismatch")
    return {
        "schema": "d6-k7-rankone-pattern-union-verification-v1",
        "status": "PASS",
        "manifest": {
            "path": manifest_path.name,
            "sha256": expected_manifest_sha256,
        },
        "graphs": len(selection),
        "tetrad_rejected": len(tetrad_rejected),
        "pattern_954_rejected": len(pattern),
        "intersection": len(both),
        "exact_union": len(union),
        "exact_residue": len(residue),
        "residue_indices_sha256": stable_hash(computed_sets["exact_residue"]),
        "tetrad_decision_rows_checked": len(rows),
        "residue_profiles_checked": len(profiles),
        "verifier_source_sha256": sha256(Path(__file__)),
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "d6_k7_rankone_pattern_union.json",
    )
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_rankone_pattern_union_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.manifest, args.manifest_sha256)
    atomic_json(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
