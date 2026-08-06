#!/usr/bin/env python3
"""Independent archive checker for a full K7 positive-dual campaign.

This checker never imports or invokes the floating-point LP locator or the
production campaign runner.  It checks the selection binding, archive/report
hashes and coverage, every retained rational polynomial identity, and repeats
the exact seed/cover quantifiers for every graph claimed newly rejected.
Survivor rows remain explicitly non-claims.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import build_d6_k7_positive_dual_selection as selection_builder
import d6_k7_rank_reference as prior
import verify_d6_k7_positive_polynomial_dual as exact


ROOT = Path(__file__).resolve().parent
DECISION_FIELDS = (
    "ordinal", "index", "status", "applicable", "seeds", "covers",
    "enhanced_passing_covers", "strict_h_passing_covers",
    "dual_passing_covers", "strict_h_passing_cliques",
    "dual_failing_cliques", "strict_h_rejected", "dual_rejected",
    "marginal_dual_rejected", "error_type",
)
COUNT_FIELDS = (
    "seeds", "covers", "enhanced_passing_covers",
    "strict_h_passing_covers", "dual_passing_covers",
    "strict_h_passing_cliques", "dual_failing_cliques",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")).hexdigest()


def resolve(path_value: str, anchor: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    local = anchor.parent / path
    return local if local.exists() else ROOT / path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_and_verify_selection(input_path: Path, selection_path: Path) -> tuple[
    list[dict], dict
]:
    recorded = json.loads(selection_path.read_text(encoding="utf-8"))
    rebuilt = selection_builder.build_selection(input_path)
    require(recorded == rebuilt, "selection does not equal an independent rebuild")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graph_by_index = {int(graph["index"]): graph for graph in payload["graphs"]}
    selected = [graph_by_index[index] for index in recorded["selected_indices"]]
    return selected, recorded


def parse_decisions(content: bytes) -> list[dict]:
    try:
        text = content.decode("ascii")
    except UnicodeDecodeError as error:
        raise AssertionError("decision archive is not ASCII") from error
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    require(tuple(reader.fieldnames or ()) == DECISION_FIELDS, "decision header mismatch")
    rows: list[dict] = []
    for raw in reader:
        ordinal = int(raw["ordinal"])
        index = int(raw["index"])
        status = raw["status"]
        require(status in ("REJECTED", "SURVIVOR", "INFRA_ERROR"), "bad status")
        if status == "INFRA_ERROR":
            require(bool(raw["error_type"]), "INFRA_ERROR lacks error type")
            require(
                all(raw[field] == "" for field in (
                    "applicable", *COUNT_FIELDS, "strict_h_rejected",
                    "dual_rejected", "marginal_dual_rejected",
                )),
                "INFRA_ERROR has mathematical fields",
            )
            rows.append({
                "ordinal": ordinal, "index": index, "status": status,
                "analysis": None, "error_type": raw["error_type"],
            })
            continue
        require(raw["error_type"] == "", "complete row has error type")
        for field in (
            "applicable", "strict_h_rejected", "dual_rejected",
            "marginal_dual_rejected",
        ):
            require(raw[field] in ("0", "1"), f"non-Boolean field {field}")
        analysis = {
            "index": index,
            "applicable": bool(int(raw["applicable"])),
            **{field: int(raw[field]) for field in COUNT_FIELDS},
            "strict_h_rejected": bool(int(raw["strict_h_rejected"])),
            "dual_rejected": bool(int(raw["dual_rejected"])),
            "marginal_dual_rejected": bool(int(raw["marginal_dual_rejected"])),
            "dual_failure_witnesses": [],
        }
        require(
            all(analysis[field] >= 0 for field in COUNT_FIELDS),
            "negative analysis count",
        )
        require(
            analysis["marginal_dual_rejected"]
            == (analysis["dual_rejected"] and not analysis["strict_h_rejected"]),
            "inconsistent marginal decision",
        )
        require(
            (status == "REJECTED") == analysis["marginal_dual_rejected"],
            "decision/flag mismatch",
        )
        rows.append({
            "ordinal": ordinal, "index": index, "status": status,
            "analysis": analysis, "error_type": None,
        })
    return rows


def parse_certificates(content: bytes) -> dict[int, dict]:
    records: dict[int, dict] = {}
    for line_number, line in enumerate(content.splitlines(), 1):
        require(bool(line), f"blank certificate line {line_number}")
        value = json.loads(line)
        require(value.get("schema") == 1, "certificate-record schema mismatch")
        index = value.get("index")
        require(
            isinstance(index, int) and not isinstance(index, bool),
            "certificate record has bad index",
        )
        require(index not in records, "duplicate certificate graph record")
        witnesses = value.get("dual_failure_witnesses")
        require(isinstance(witnesses, list) and witnesses, "empty certificate record")
        records[index] = value
    return records


def verify_witness(graph: dict, witness: dict) -> None:
    adjacency = tuple(graph["adjacency"])
    nvertices = witness.get("nvertices")
    require(
        isinstance(nvertices, list)
        and len(nvertices) == len(set(nvertices))
        and all(
            isinstance(vertex, int) and 0 <= vertex < len(adjacency)
            for vertex in nvertices
        ),
        "invalid induced-vertex witness",
    )
    graph_n = prior.induced_graph(adjacency, nvertices)
    certificate = witness.get("certificate")
    require(isinstance(certificate, dict), "missing rational certificate")
    clique = certificate.get("clique")
    require(isinstance(clique, list), "certificate has no clique")
    clique_mask = sum(1 << vertex for vertex in clique)
    exact.verify_clique_certificate(tuple(graph_n), clique_mask, certificate)


def expected_summary(rows: list[dict]) -> dict:
    complete = [row for row in rows if row["analysis"] is not None]
    analyses = [row["analysis"] for row in complete]
    rejected = sorted(
        row["index"] for row in complete
        if row["analysis"]["marginal_dual_rejected"]
    )
    strict_rejected = sorted(
        row["index"] for row in complete
        if row["analysis"]["strict_h_rejected"]
    )
    return {
        "graphs": len(rows),
        "complete": len(complete),
        "infra_errors": len(rows) - len(complete),
        "marginal_dual_rejected": len(rejected),
        "survivors": len(complete) - len(rejected),
        "strict_h_rejected": len(strict_rejected),
        "marginal_dual_rejected_indices": rejected,
        "strict_h_rejected_indices": strict_rejected,
        "certificate_graphs": sum(
            bool(analysis["dual_failure_witnesses"]) for analysis in analyses
        ),
        "certificate_witnesses": sum(
            len(analysis["dual_failure_witnesses"]) for analysis in analyses
        ),
        "totals": {
            field: sum(analysis[field] for analysis in analyses)
            for field in COUNT_FIELDS
        },
    }


def verify_campaign(
    input_path: Path,
    selection_path: Path,
    report_path: Path,
    expected_report_sha256: str | None = None,
) -> dict:
    if expected_report_sha256 is not None:
        require(sha256(report_path) == expected_report_sha256, "report hash mismatch")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    require(report.get("schema") == 1, "report schema mismatch")
    configuration = report.get("configuration")
    require(isinstance(configuration, dict), "report has no configuration")
    require(
        report.get("config_sha256") == stable_hash(configuration),
        "configuration hash mismatch",
    )
    selected_full, selection = load_and_verify_selection(input_path, selection_path)
    require(
        configuration.get("selection", {}).get("input_sha256") == sha256(input_path)
        and configuration.get("selection", {}).get("selection_report_sha256")
        == sha256(selection_path)
        and configuration.get("full_selection_graphs") == len(selected_full)
        and configuration.get("full_selection_indices_sha256")
        == selection["selected_indices_sha256"],
        "campaign/selection binding mismatch",
    )
    for name, digest in configuration.get("dependencies_sha256", {}).items():
        require(sha256(ROOT / name) == digest, f"dependency hash mismatch: {name}")
    runner = ROOT / "run_d6_k7_positive_polynomial_dual_full.py"
    require(
        sha256(runner) == configuration.get("runner_source_sha256"),
        "runner source hash mismatch",
    )

    start = configuration.get("start")
    limit = configuration.get("limit")
    require(isinstance(start, int) and start >= 0, "invalid run start")
    stop = len(selected_full) if limit is None else start + int(limit)
    run_graphs = selected_full[start:stop]
    run_indices = [int(graph["index"]) for graph in run_graphs]
    require(
        len(run_graphs) == configuration.get("run_graphs")
        and stable_hash(run_indices) == configuration.get("run_indices_sha256"),
        "run stratum binding mismatch",
    )

    artifacts = report.get("artifacts", {})
    decisions_path = resolve(artifacts["decisions_archive"], report_path)
    certificates_path = resolve(artifacts["certificate_archive"], report_path)
    require(sha256(decisions_path) == artifacts["decisions_archive_sha256"],
            "decision archive hash mismatch")
    require(sha256(certificates_path) == artifacts["certificate_archive_sha256"],
            "certificate archive hash mismatch")
    with gzip.open(decisions_path, "rb") as stream:
        decision_bytes = stream.read()
    with gzip.open(certificates_path, "rb") as stream:
        certificate_bytes = stream.read()
    require(
        hashlib.sha256(decision_bytes).hexdigest()
        == artifacts["decisions_uncompressed_sha256"],
        "uncompressed decision hash mismatch",
    )
    require(
        hashlib.sha256(certificate_bytes).hexdigest()
        == artifacts["certificates_uncompressed_sha256"],
        "uncompressed certificate hash mismatch",
    )

    rows = parse_decisions(decision_bytes)
    require(len(rows) == artifacts["decision_rows"], "decision row count mismatch")
    require([row["ordinal"] for row in rows] == list(range(len(rows))),
            "decision ordinals are incomplete/out of order")
    require([row["index"] for row in rows] == run_indices,
            "decision graph coverage/order mismatch")
    graph_by_index = {int(graph["index"]): graph for graph in run_graphs}
    certificate_records = parse_certificates(certificate_bytes)
    require(
        set(certificate_records).issubset(graph_by_index),
        "certificate archive contains a graph outside this run",
    )
    for row in rows:
        record = certificate_records.get(row["index"])
        witnesses = [] if record is None else record["dual_failure_witnesses"]
        if record is not None:
            require(record["ordinal"] == row["ordinal"], "certificate ordinal mismatch")
        if row["analysis"] is None:
            require(not witnesses, "INFRA_ERROR has certificate witnesses")
            continue
        require(
            len(witnesses) == row["analysis"]["dual_failing_cliques"],
            "certificate witness count mismatch",
        )
        row["analysis"]["dual_failure_witnesses"] = witnesses
        for witness in witnesses:
            verify_witness(graph_by_index[row["index"]], witness)

    rebuilt_summary = expected_summary(rows)
    require(report.get("summary") == rebuilt_summary, "aggregate summary mismatch")
    require(rebuilt_summary["infra_errors"] == 0,
            "campaign has infrastructure errors and is not proof-complete")
    independently_rebuilt = 0
    for row in rows:
        if row["status"] != "REJECTED":
            continue
        exact.verify_graph(graph_by_index[row["index"]], row["analysis"])
        independently_rebuilt += 1
    return {
        "graphs": len(rows),
        "certificate_graphs": len(certificate_records),
        "certificate_witnesses": rebuilt_summary["certificate_witnesses"],
        "marginal_dual_rejected": rebuilt_summary["marginal_dual_rejected"],
        "rejected_graph_quantifiers_rebuilt": independently_rebuilt,
        "survivors_are_nonclaims": rebuilt_summary["survivors"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(".runs/d6_k7_rank_survivors.json")
    )
    parser.add_argument(
        "--selection", type=Path,
        default=ROOT / "d6_k7_positive_dual_selection.json",
    )
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k7_positive_dual_full_report.json",
    )
    parser.add_argument("--report-sha256")
    args = parser.parse_args()
    summary = verify_campaign(
        args.input, args.selection, args.report, args.report_sha256
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
