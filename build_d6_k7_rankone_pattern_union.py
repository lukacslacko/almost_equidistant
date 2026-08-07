#!/usr/bin/env python3
"""Build the exact union of the K7 tetrad and pattern-954 rejections.

The two source campaigns attack the same frozen 12,839-graph K7 residue.
This builder accepts only hash-pinned, complete, independently verified
campaigns.  It parses every tetrad decision row, reconstructs all set
operations in the frozen selection order, and records the structure of the
remaining tetrad covers.  ``SURVIVOR`` never means realizable.
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
SELECTION = ROOT / "d6_k7_rankone_tetrad_full_selection.json"
TETRAD_REPORT = ROOT / "d6_k7_rankone_tetrad_full_report.json"
TETRAD_VERIFICATION = (
    ROOT / "d6_k7_rankone_tetrad_full_verification_report.json"
)
TETRAD_DECISIONS = ROOT / "d6_k7_rankone_tetrad_full_decisions.tsv.gz"
PATTERN_REPORT = ROOT / "d6_n14_pattern_954_containment_report.json"
PATTERN_VERIFICATION = (
    ROOT / "d6_n14_pattern_954_containment_verification.json"
)
OUTPUT = ROOT / "d6_k7_rankone_pattern_union.json"
VERIFIER = ROOT / "verify_d6_k7_rankone_pattern_union.py"

SELECTION_SHA256 = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
SELECTION_INDICES_SHA256 = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)
PATTERN_REPORT_SHA256 = (
    "a2d4fb5a07dcea580315730b4155fc0c08ecd0439a9111e733d29e2da0b5be62"
)
PATTERN_VERIFICATION_SHA256 = (
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
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path.name} is not a JSON object")
    return value


def require_file_hash(path: Path, expected: str) -> str:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(
            f"hash mismatch for {path.name}: {observed} != {expected}"
        )
    return observed


def integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise ValueError(f"{label} is not an integer list")
    if len(set(value)) != len(value):
        raise ValueError(f"{label} contains duplicates")
    return list(value)


def parse_bool(value: str, label: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"{label} is not a serialized Boolean: {value!r}")


def ordered(selection: Sequence[int], members: set[int]) -> list[int]:
    return [index for index in selection if index in members]


def set_record(indices: list[int], include: bool = True) -> dict:
    result = {
        "count": len(indices),
        "indices_sha256": stable_hash(indices),
    }
    if include:
        result["indices"] = indices
    return result


def cover_profile(row: dict) -> dict:
    prior = row["prior_dual_passing_covers"]
    no_saturating = row["prior_passing_no_saturating_clique_covers"]
    saturating = prior - no_saturating
    with_near = row["covers_with_near_clique"]
    without_near = row["covers_without_near_clique"]
    failed = row["tetrad_failed_covers"]
    near_resistant = with_near - failed
    if min(saturating, no_saturating, with_near, without_near, near_resistant) < 0:
        raise ValueError(f"negative derived cover count at index {row['index']}")
    if with_near + without_near != no_saturating:
        raise ValueError(f"near-clique partition mismatch at index {row['index']}")
    if row["tetrad_passing_covers"] != (
        saturating + without_near + near_resistant
    ):
        raise ValueError(f"tetrad passing-cover mismatch at index {row['index']}")
    if no_saturating == 0:
        structure = "saturating_only"
    elif saturating == 0:
        structure = "no_saturating_only"
    else:
        structure = "mixed"
    return {
        "structure": structure,
        "prior_passing_covers": prior,
        "saturating_covers": saturating,
        "no_saturating_covers": no_saturating,
        "no_saturating_with_near_clique": with_near,
        "no_saturating_without_near_clique": without_near,
        "tetrad_killed_covers": failed,
        "tetrad_resistant_near_clique_covers": near_resistant,
        "tetrad_passing_covers": row["tetrad_passing_covers"],
        "has_saturating_cover": saturating > 0,
        "has_no_near_clique_cover": without_near > 0,
        "has_tetrad_resistant_near_clique_cover": near_resistant > 0,
    }


def read_decisions(
    selection: Sequence[int], report: dict
) -> tuple[list[dict], dict, str]:
    artifacts = report.get("artifacts", {})
    recorded_path = artifacts.get("decisions_archive")
    if recorded_path != TETRAD_DECISIONS.name:
        raise ValueError(
            f"unexpected tetrad decision path in report: {recorded_path!r}"
        )
    compressed_hash = sha256(TETRAD_DECISIONS)
    if compressed_hash != artifacts.get("decisions_archive_sha256"):
        raise ValueError("tetrad decision compressed hash mismatch")
    with gzip.open(TETRAD_DECISIONS, "rb") as stream:
        raw = stream.read()
    raw_hash = hashlib.sha256(raw).hexdigest()
    if raw_hash != artifacts.get("decisions_uncompressed_sha256"):
        raise ValueError("tetrad decision uncompressed hash mismatch")
    reader = csv.DictReader(
        io.StringIO(raw.decode("utf-8"), newline=""), delimiter="\t"
    )
    if reader.fieldnames != DECISION_FIELDS:
        raise ValueError(f"unexpected tetrad decision schema: {reader.fieldnames}")
    raw_rows = list(reader)
    if len(raw_rows) != len(selection):
        raise ValueError("tetrad decision row count mismatch")

    rows = []
    totals = Counter()
    statuses = Counter()
    for ordinal, (raw_row, expected_index) in enumerate(zip(raw_rows, selection)):
        if int(raw_row["ordinal"]) != ordinal:
            raise ValueError(f"tetrad decision ordinal mismatch at {ordinal}")
        index = int(raw_row["index"])
        if index != expected_index:
            raise ValueError(f"tetrad decision index mismatch at {ordinal}")
        status = raw_row["status"]
        if status not in {"REJECTED", "SURVIVOR"} or raw_row["error_type"]:
            raise ValueError(f"bad tetrad decision status at index {index}")
        counts = {name: int(raw_row[name]) for name in COUNT_FIELDS}
        if any(value < 0 for value in counts.values()):
            raise ValueError(f"negative tetrad decision count at index {index}")
        prior_rejected = parse_bool(
            raw_row["prior_dual_rejected"], "prior_dual_rejected"
        )
        tetrad_rejected = parse_bool(
            raw_row["tetrad_rejected"], "tetrad_rejected"
        )
        marginal = parse_bool(
            raw_row["marginal_tetrad_rejected"],
            "marginal_tetrad_rejected",
        )
        if prior_rejected or tetrad_rejected != marginal:
            raise ValueError(f"non-marginal tetrad row at index {index}")
        if (status == "REJECTED") != marginal:
            raise ValueError(f"tetrad status/Boolean mismatch at index {index}")
        row = {
            "ordinal": ordinal,
            "index": index,
            "status": status,
            **counts,
        }
        profile = cover_profile(row)
        row["profile"] = profile
        rows.append(row)
        totals.update(counts)
        statuses[status] += 1

    report_summary = report.get("summary", {})
    if dict(statuses) != report_summary.get("status_counts"):
        raise ValueError("tetrad report status totals disagree with archive")
    if {name: totals[name] for name in COUNT_FIELDS} != report_summary.get(
        "totals"
    ):
        raise ValueError("tetrad report cover totals disagree with archive")
    return rows, {
        "rows": len(rows),
        "status_counts": dict(statuses),
        "uncompressed_sha256": raw_hash,
    }, compressed_hash


def distribution(rows: Sequence[dict], members: set[int]) -> dict:
    selected = [row for row in rows if row["index"] in members]
    structure = Counter(row["profile"]["structure"] for row in selected)
    flags = {
        name: sum(bool(row["profile"][name]) for row in selected)
        for name in (
            "has_saturating_cover",
            "has_no_near_clique_cover",
            "has_tetrad_resistant_near_clique_cover",
        )
    }
    return {
        "graphs": len(selected),
        "structure": {
            name: structure.get(name, 0)
            for name in ("saturating_only", "mixed", "no_saturating_only")
        },
        "obstacle_flags": flags,
    }


def validate_tetrad_gate(
    report: dict,
    report_hash: str,
    verification: dict,
    verification_hash: str,
    selection: Sequence[int],
) -> list[int]:
    configuration = report.get("configuration", {})
    report_selection = configuration.get("selection", {})
    summary = report.get("summary", {})
    if configuration.get("limit") is not None:
        raise ValueError("tetrad report is a limited run")
    if (
        report_selection.get("selection_report_sha256") != SELECTION_SHA256
        or report_selection.get("selected_indices_sha256")
        != SELECTION_INDICES_SHA256
        or report_selection.get("selected") != EXPECTED_GRAPHS
    ):
        raise ValueError("tetrad report selection provenance mismatch")
    rejected = integer_list(
        summary.get("marginal_tetrad_rejected_indices"),
        "tetrad rejected indices",
    )
    if (
        summary.get("schema") != 1
        or summary.get("graphs") != len(selection)
        or summary.get("complete") != len(selection)
        or summary.get("infra_errors") != 0
        or summary.get("marginal_tetrad_rejected") != len(rejected)
        or summary.get("status_counts", {}).get("REJECTED", 0) != len(rejected)
        or sum(summary.get("status_counts", {}).values()) != len(selection)
    ):
        raise ValueError("tetrad report is incomplete or internally inconsistent")
    if not set(rejected) <= set(selection):
        raise ValueError("tetrad report rejects an index outside the selection")
    if rejected != ordered(selection, set(rejected)):
        raise ValueError("tetrad rejected indices are not in selection order")

    provenance = verification.get("provenance", {})
    vsummary = verification.get("summary", {})
    if (
        verification.get("schema") != 1
        or verification.get("errors") != []
        or provenance.get("selection_sha256") != SELECTION_SHA256
        or provenance.get("tetrad_report_sha256") != report_hash
        or provenance.get("tetrad_decisions_sha256")
        != report.get("artifacts", {}).get("decisions_archive_sha256")
        or vsummary.get("graphs") != len(selection)
        or vsummary.get("status_counts") != {"PASS": len(selection)}
        or vsummary.get("verified_rejections") != len(rejected)
    ):
        raise ValueError("independent tetrad verification gate failed")
    if sha256(TETRAD_VERIFICATION) != verification_hash:
        raise ValueError("tetrad verification hash changed during build")
    return rejected


def validate_pattern_gate(
    report: dict, verification: dict, selection: Sequence[int]
) -> list[int]:
    hits = integer_list(
        report.get("certified_hit_indices"), "pattern certified hits"
    )
    counts = report.get("counts", {})
    if (
        report.get("schema") != 1
        or report.get("status") != "COMPLETE"
        or report.get("required_edge_only") is not True
        or report.get("processed") != len(selection)
        or counts != {
            "HIT": len(hits),
            "INFRA_ERROR": 0,
            "NO_HIT": len(selection) - len(hits),
            "TIMEOUT": 0,
        }
        or report.get("certified_hit_indices_sha256") != stable_hash(hits)
        or report.get("unresolved_indices") != []
        or report.get("unresolved_indices_sha256") != stable_hash([])
    ):
        raise ValueError("pattern report is incomplete or internally inconsistent")
    if hits != ordered(selection, set(hits)):
        raise ValueError("pattern hits are not exactly in selection order")

    verified_hits = integer_list(
        verification.get("certified_hit_indices"), "verified pattern hits"
    )
    verified_residue = integer_list(
        verification.get("residue_indices"), "verified pattern residue"
    )
    mapping = verification.get("mapping_verification", {})
    selection_record = verification.get("selection", {})
    report_record = verification.get("report", {})
    if (
        verification.get("schema") != 1
        or verification.get("status") != "PASS"
        or verification.get("required_edge_only") is not True
        or verified_hits != hits
        or verification.get("certified_hit_indices_sha256") != stable_hash(hits)
        or verified_residue != ordered(selection, set(selection) - set(hits))
        or verification.get("residue_indices_sha256")
        != stable_hash(verified_residue)
        or mapping.get("all_injective") is not True
        or mapping.get("mappings_checked") != len(hits)
        or selection_record.get("population") != len(selection)
        or selection_record.get("indices_sha256") != SELECTION_INDICES_SHA256
        or report_record.get("sha256") != PATTERN_REPORT_SHA256
        or report_record.get("processed") != len(selection)
        or report_record.get("counts") != counts
    ):
        raise ValueError("independent pattern verification gate failed")
    return hits


def build_manifest(
    tetrad_report_sha256: str, tetrad_verification_sha256: str
) -> dict:
    source_hashes = {
        SELECTION.name: require_file_hash(SELECTION, SELECTION_SHA256),
        TETRAD_REPORT.name: require_file_hash(
            TETRAD_REPORT, tetrad_report_sha256
        ),
        TETRAD_VERIFICATION.name: require_file_hash(
            TETRAD_VERIFICATION, tetrad_verification_sha256
        ),
        PATTERN_REPORT.name: require_file_hash(
            PATTERN_REPORT, PATTERN_REPORT_SHA256
        ),
        PATTERN_VERIFICATION.name: require_file_hash(
            PATTERN_VERIFICATION, PATTERN_VERIFICATION_SHA256
        ),
        Path(__file__).name: sha256(Path(__file__)),
        VERIFIER.name: sha256(VERIFIER),
    }
    selection_payload = load_json(SELECTION)
    selection = integer_list(
        selection_payload.get("selected_indices"), "frozen selection"
    )
    if (
        selection_payload.get("schema") != 1
        or selection_payload.get("kind")
        != "d6_k7_rankone_tetrad_full_selection"
        or len(selection) != EXPECTED_GRAPHS
        or selection_payload.get("selected_indices_sha256")
        != SELECTION_INDICES_SHA256
        or stable_hash(selection) != SELECTION_INDICES_SHA256
    ):
        raise ValueError("frozen tetrad selection mismatch")

    tetrad_report = load_json(TETRAD_REPORT)
    tetrad_verification = load_json(TETRAD_VERIFICATION)
    tetrad_rejected_list = validate_tetrad_gate(
        tetrad_report,
        tetrad_report_sha256,
        tetrad_verification,
        tetrad_verification_sha256,
        selection,
    )
    rows, decision_audit, decisions_hash = read_decisions(
        selection, tetrad_report
    )
    source_hashes[TETRAD_DECISIONS.name] = decisions_hash
    archive_rejected = [
        row["index"] for row in rows if row["status"] == "REJECTED"
    ]
    if archive_rejected != tetrad_rejected_list:
        raise ValueError("tetrad report and decision archive reject different sets")

    pattern_report = load_json(PATTERN_REPORT)
    pattern_verification = load_json(PATTERN_VERIFICATION)
    pattern_hits_list = validate_pattern_gate(
        pattern_report, pattern_verification, selection
    )

    tetrad_rejected = set(tetrad_rejected_list)
    pattern_hits = set(pattern_hits_list)
    both = tetrad_rejected & pattern_hits
    tetrad_only = tetrad_rejected - pattern_hits
    pattern_only = pattern_hits - tetrad_rejected
    union = tetrad_rejected | pattern_hits
    universe = set(selection)
    residue = universe - union
    tetrad_survivors = universe - tetrad_rejected

    ordered_sets = {
        "tetrad_rejected": ordered(selection, tetrad_rejected),
        "pattern_954_rejected": ordered(selection, pattern_hits),
        "both": ordered(selection, both),
        "tetrad_only": ordered(selection, tetrad_only),
        "pattern_954_only": ordered(selection, pattern_only),
        "exact_union": ordered(selection, union),
        "exact_residue": ordered(selection, residue),
    }
    residue_profiles = []
    for row in rows:
        if row["index"] in residue:
            residue_profiles.append({
                "index": row["index"],
                **row["profile"],
            })

    return {
        "schema": "d6-k7-rankone-pattern-union-v1",
        "status": "COMPLETE",
        "claim": (
            "The exact union consists only of independently verified tetrad "
            "rejections and injectively rechecked pattern-954 containments. "
            "The exact residue is their ordered complement in the frozen "
            "12,839-graph K7 selection; residue membership is not a "
            "realizability claim."
        ),
        "source_hashes": source_hashes,
        "gates": {
            "selection": {
                "graphs": len(selection),
                "indices_sha256": stable_hash(selection),
            },
            "tetrad": {
                "report_sha256": tetrad_report_sha256,
                "verification_sha256": tetrad_verification_sha256,
                "verification_status_counts": tetrad_verification["summary"][
                    "status_counts"
                ],
                "verified_rejections": tetrad_verification["summary"][
                    "verified_rejections"
                ],
                "decision_archive_audit": decision_audit,
            },
            "pattern_954": {
                "report_sha256": PATTERN_REPORT_SHA256,
                "verification_sha256": PATTERN_VERIFICATION_SHA256,
                "mappings_rechecked": pattern_verification[
                    "mapping_verification"
                ]["mappings_checked"],
                "all_mappings_injective": pattern_verification[
                    "mapping_verification"
                ]["all_injective"],
            },
        },
        "sets": {
            name: set_record(indices) for name, indices in ordered_sets.items()
        },
        "overlap": {
            "tetrad_rejected": len(tetrad_rejected),
            "pattern_954_rejected": len(pattern_hits),
            "both": len(both),
            "tetrad_only": len(tetrad_only),
            "pattern_954_only": len(pattern_only),
            "exact_union": len(union),
            "exact_residue": len(residue),
        },
        "cover_structure": {
            "all_selected": distribution(rows, universe),
            "tetrad_survivors": distribution(rows, tetrad_survivors),
            "exact_union_residue": distribution(rows, residue),
            "pattern_hits_among_tetrad_survivors": distribution(
                rows, pattern_only
            ),
            "residue_profiles": residue_profiles,
        },
        "trust_scope": {
            "exact": (
                "set arithmetic, decision-archive parsing, all tetrad "
                "rational identities accepted by the independent full "
                "verifier, and every positive pattern embedding rechecked "
                "edge by edge"
            ),
            "not_a_claim": (
                "absence of either certificate, or realizability of any "
                "SURVIVOR / exact-residue graph"
            ),
            "floating_point_used_for_union_membership": False,
        },
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tetrad-report-sha256", required=True)
    parser.add_argument("--tetrad-verification-sha256", required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    manifest = build_manifest(
        args.tetrad_report_sha256, args.tetrad_verification_sha256
    )
    atomic_json(args.output, manifest)
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "overlap": manifest["overlap"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
