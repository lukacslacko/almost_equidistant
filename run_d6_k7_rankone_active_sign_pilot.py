#!/usr/bin/env python3
"""One-process cover-only pilot for the near-clique active/sign layer.

The deterministic sample contains 12 covers with a saturating clique and 12
without one, all surviving the completed degree-four campaign and all having
an applicable near-saturating clique.  This script reconstructs upstream
cover state but deliberately does not assert graph coverage: it may stop
partway through a graph's cover list, so the outer graph quantifier is not
complete.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

import d6_k7_rank_reference as prior
import d6_k7_rankone_active_sign as active_sign
import d6_k7_special_h_reference as strict_h
import run_d6_k7_rankone_tetrad_full as full
import verify_d6_k7_rankone_active_sign as active_sign_verify


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def survivor_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        return [
            row for row in csv.DictReader(stream, delimiter="\t")
            if row["status"] == "SURVIVOR"
        ]


SATURATING = "SATURATING_COVER"
NO_SATURATING = "NO_SATURATING_TETRAD_SURVIVOR"


def reconstructed_campaign_survivor_covers(
    graph: dict,
    prior_witnesses: list[dict],
    checkpoint: dict,
):
    """Yield applicable covers not rejected by any earlier/tetrad witness."""

    adj = tuple(int(value) for value in graph["adjacency"])
    prior.validate_graph(adj)
    old_witness_keys = {
        (tuple(item["seed"]), int(item["zmask"]))
        for item in prior_witnesses
    }
    tetrad_failure_keys = {
        (tuple(item["seed"]), int(item["zmask"]))
        for item in checkpoint["record"]["result"]["tetrad_failure_witnesses"]
    }
    support_solver = prior.SupportSolver()
    zero_forcing = prior.ZeroForcingSolver()
    clique_solver = prior.CliqueStructureSolver()

    for seed_mask in prior.clique_masks(adj, 7):
        seed, outside, defects, ladj, eligible = prior.seed_instance(adj, seed_mask)
        total_term_rank = prior.matching_size(defects)
        for zmask in prior.eligible_covers(ladj, eligible):
            analysis = prior.analyze_cover(
                adj,
                outside,
                defects,
                zmask,
                support_solver,
                zero_forcing,
                total_term_rank,
                clique_solver,
            )
            if analysis.enhanced_joint_failed:
                continue
            nvertices = [
                outside[index]
                for index in range(len(outside))
                if not (zmask & (1 << index))
            ]
            graph_n = prior.induced_graph(adj, nvertices)
            saturating = tuple(prior.clique_masks(
                graph_n, analysis.k_rank_upper
            )) if analysis.k_rank_upper else ()
            if any(
                strict_h.assess_saturating_clique(graph_n, mask).failed
                for mask in saturating
            ):
                continue
            key = (tuple(seed), zmask)
            if key in old_witness_keys:
                continue
            if saturating:
                near_size = analysis.k_rank_upper - 1
                near_cliques = tuple(prior.clique_masks(
                    graph_n, near_size
                )) if near_size > 0 else ()
                if not near_cliques:
                    raise ValueError("saturating clique has no near subclique")
                yield {
                    "stratum": SATURATING,
                    "seed": seed,
                    "zmask": zmask,
                    "nvertices": nvertices,
                    "graph_n": graph_n,
                    "k_rank_upper": analysis.k_rank_upper,
                    "near_cliques": near_cliques,
                }
                continue
            near_size = analysis.k_rank_upper - 1
            near_cliques = tuple(prior.clique_masks(
                graph_n, near_size
            )) if near_size > 0 else ()
            if not near_cliques or key in tetrad_failure_keys:
                continue
            yield {
                "stratum": NO_SATURATING,
                "seed": seed,
                "zmask": zmask,
                "nvertices": nvertices,
                "graph_n": graph_n,
                "k_rank_upper": analysis.k_rank_upper,
                "near_cliques": near_cliques,
            }


def assess_cover(graph_index: int, cover: dict) -> dict:
    clique_records = []
    cover_rejected = False
    for clique_mask in cover["near_cliques"]:
        variables, vertices, masks, targets = active_sign.graph_system(
            cover["graph_n"], clique_mask
        )
        assessment = active_sign.assess_system(
            variables, vertices, masks, targets
        )
        record = {
            "clique": [
                vertex for vertex in range(len(cover["graph_n"]))
                if clique_mask & (1 << vertex)
            ],
            "rejected": assessment.rejected,
            "reason": assessment.reason,
            "profile": assessment.profile,
            "diagonal_forced_active": list(
                assessment.diagonal_forced_active
            ),
            "forced_active": list(assessment.forced_active),
        }
        if assessment.rejected:
            certificate = active_sign.make_certificate(
                variables, vertices, masks, targets, assessment
            )
            verification = active_sign_verify.verify_certificate(certificate)
            record["certificate"] = certificate
            record["independent_verification"] = verification
            cover_rejected = True
        else:
            record["relaxed_active"] = list(assessment.relaxed_active or ())
            record["relaxed_signs"] = [
                list(item) for item in (assessment.relaxed_signs or ())
            ]
        clique_records.append(record)
    return {
        "graph_index": graph_index,
        "stratum": cover["stratum"],
        "seed": list(cover["seed"]),
        "zmask": int(cover["zmask"]),
        "nvertices": list(cover["nvertices"]),
        "k_rank_upper": int(cover["k_rank_upper"]),
        "near_cliques": clique_records,
        "near_clique_count": len(clique_records),
        "active_sign_rejected": cover_rejected,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-stratum", type=int, default=12)
    parser.add_argument(
        "--rank-input", type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    parser.add_argument(
        "--selection", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_selection.json"),
    )
    parser.add_argument(
        "--selection-sha256",
        default="86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479",
    )
    parser.add_argument(
        "--decisions", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_decisions.tsv.gz"),
    )
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=Path(".runs/d6_k7_rankone_tetrad_full"),
    )
    parser.add_argument(
        "--full-report", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_report.json"),
    )
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_rankone_active_sign_pilot_report.json"),
    )
    args = parser.parse_args()
    if not 1 <= args.per_stratum <= 12:
        raise ValueError("this pilot requires 1 <= per-stratum <= 12")

    started = time.perf_counter()
    report_payload = json.loads(args.full_report.read_text())
    expected_config = report_payload["configuration"]["config_sha256"]
    selected, selection_provenance, _, prior_certificates = full.load_inputs(
        args.rank_input, args.selection, args.selection_sha256
    )
    selected_by_index = {int(graph["index"]): graph for graph in selected}
    rows = survivor_rows(args.decisions)

    cover_records = []
    selected_by_stratum: Counter[str] = Counter()
    quotas = {
        SATURATING: args.per_stratum,
        NO_SATURATING: args.per_stratum,
    }
    inspected_graphs = []
    for row in rows:
        if all(selected_by_stratum[name] >= quota for name, quota in quotas.items()):
            break
        graph_index = int(row["index"])
        graph = selected_by_index[graph_index]
        ordinal = int(row["ordinal"])
        checkpoint_path = args.checkpoint_dir / f"graph_{ordinal:05d}_{graph_index}.json"
        checkpoint = json.loads(checkpoint_path.read_text())
        if checkpoint.get("config_sha256") != expected_config:
            raise ValueError("checkpoint/full-report configuration mismatch")
        record = checkpoint.get("record", {})
        if (
            record.get("status") != "SURVIVOR"
            or int(record.get("index")) != graph_index
            or int(record.get("ordinal")) != ordinal
        ):
            raise ValueError("decision/checkpoint survivor mismatch")
        inspected_graphs.append(graph_index)
        for cover in reconstructed_campaign_survivor_covers(
            graph,
            prior_certificates.get(graph_index, []),
            checkpoint,
        ):
            stratum = cover["stratum"]
            if selected_by_stratum[stratum] >= quotas[stratum]:
                continue
            cover_records.append(assess_cover(graph_index, cover))
            selected_by_stratum[stratum] += 1
            if all(
                selected_by_stratum[name] >= quota
                for name, quota in quotas.items()
            ):
                break

    expected_records = 2 * args.per_stratum
    if len(cover_records) != expected_records:
        raise ValueError(
            f"only reconstructed {len(cover_records)} applicable covers, "
            f"expected {expected_records}: {dict(selected_by_stratum)}"
        )

    pair_profile: Counter[str] = Counter()
    near_cliques = rejected_near_cliques = verified_certificates = 0
    for cover in cover_records:
        for clique in cover["near_cliques"]:
            near_cliques += 1
            pair_profile.update(clique["profile"])
            if clique["rejected"]:
                rejected_near_cliques += 1
                verified_certificates += 1
    rejected_covers = sum(
        int(cover["active_sign_rejected"]) for cover in cover_records
    )
    stratum_coverage = {}
    for stratum in (SATURATING, NO_SATURATING):
        chosen = [cover for cover in cover_records if cover["stratum"] == stratum]
        stratum_cliques = [
            clique for cover in chosen for clique in cover["near_cliques"]
        ]
        stratum_profile: Counter[str] = Counter()
        for clique in stratum_cliques:
            stratum_profile.update(clique["profile"])
        stratum_coverage[stratum] = {
            "covers_tested": len(chosen),
            "covers_rejected": sum(
                int(cover["active_sign_rejected"]) for cover in chosen
            ),
            "near_cliques_tested": len(stratum_cliques),
            "near_cliques_rejected": sum(
                int(clique["rejected"]) for clique in stratum_cliques
            ),
            "near_cliques_with_diagonal_forcing": sum(
                bool(clique["diagonal_forced_active"])
                for clique in stratum_cliques
            ),
            "diagonal_forced_active_incidences": sum(
                len(clique["diagonal_forced_active"])
                for clique in stratum_cliques
            ),
            "pair_classification_totals": dict(sorted(
                stratum_profile.items()
            )),
        }
    output = {
        "schema": 1,
        "kind": "d6_k7_rankone_active_sign_cover_only_pilot",
        "description": (
            "Deterministic 12+12 stratified applicable covers among completed "
            "degree-four campaign survivors; no graph quantifier claim."
        ),
        "provenance": {
            "rank_and_selection": selection_provenance,
            "full_report": {
                "path": str(args.full_report),
                "sha256": sha256(args.full_report),
                "config_sha256": expected_config,
            },
            "decisions": {
                "path": str(args.decisions),
                "sha256": sha256(args.decisions),
            },
            "checkpoint_directory": str(args.checkpoint_dir),
            "locator_source_sha256": sha256(ROOT / "d6_k7_rankone_active_sign.py"),
            "verifier_source_sha256": sha256(ROOT / "verify_d6_k7_rankone_active_sign.py"),
            "pilot_source_sha256": sha256(Path(__file__).resolve()),
        },
        "selection": {
            "rule": (
                "survivor rows and covers in degree-four decision order; "
                "reconstruct earlier filters; take the first quota from both "
                "saturating and no-saturating applicable strata"
            ),
            "per_stratum": args.per_stratum,
            "total_limit": expected_records,
            "inspected_graph_indices": inspected_graphs,
        },
        "coverage": {
            "covers_tested": len(cover_records),
            "covers_rejected": rejected_covers,
            "covers_surviving": len(cover_records) - rejected_covers,
            "near_cliques_tested": near_cliques,
            "near_cliques_rejected": rejected_near_cliques,
            "near_cliques_with_diagonal_forcing": sum(
                bool(clique["diagonal_forced_active"])
                for cover in cover_records
                for clique in cover["near_cliques"]
            ),
            "diagonal_forced_active_incidences": sum(
                len(clique["diagonal_forced_active"])
                for cover in cover_records
                for clique in cover["near_cliques"]
            ),
            "independently_verified_abstract_certificates": verified_certificates,
            "pair_classification_totals": dict(sorted(pair_profile.items())),
            "by_stratum": stratum_coverage,
            "graph_coverage": {
                "status": "NOT_ASSESSED",
                "certified_rejections": 0,
                "reason": (
                    "the 24-cover truncation can stop inside a graph; full "
                    "K7 seed and eligible-cover quantifiers were not checked"
                ),
            },
        },
        "records": cover_records,
        "runtime": {
            "processes": 1,
            "wall_seconds": time.perf_counter() - started,
        },
        "claim": (
            "Each rejected record is an abstract post-cover near-clique "
            "contradiction checked independently. This pilot makes no graph "
            "rejection claim."
        ),
    }
    args.report.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output["coverage"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
