#!/usr/bin/env python3
"""Run the first-core exact K5/rank-two Schur pentad pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import d6_k7_ranktwo_pentad as pentad


ROOT = Path(__file__).resolve().parent
EXPECTED_BASIS_PILOT_SHA256 = (
    "2c2dd18ded8a0daf505700001e64eb31fdcbf987f6292f91e58a3c07c8dc8339"
)
EXPECTED_BASIS_VERIFICATION_SHA256 = (
    "789f53f235599b9b4ce0c4b8b288f7fb1ed767d9038d340dff861d18cfee2d44"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        ).stdout.strip()

    try:
        commit = git("rev-parse", "HEAD")
        branch = git("branch", "--show-current")
        status = git("status", "--porcelain=v1", "--untracked-files=all")
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}
    return {
        "available": True,
        "commit": commit,
        "branch": branch,
        "dirty": bool(status),
        "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
        "porcelain_lines": status.splitlines(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--basis-pilot", type=Path,
        default=Path("d6_k7_arbitrary_basis_extension_pilot_report.json"),
    )
    parser.add_argument(
        "--basis-verification", type=Path,
        default=Path("d6_k7_arbitrary_basis_extension_pilot_verification.json"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_ranktwo_pentad_pilot_report.json"),
    )
    args = parser.parse_args()
    started = time.perf_counter()

    basis_hash = sha256(args.basis_pilot)
    verification_hash = sha256(args.basis_verification)
    if basis_hash != EXPECTED_BASIS_PILOT_SHA256:
        raise ValueError("basis pilot report hash mismatch")
    if verification_hash != EXPECTED_BASIS_VERIFICATION_SHA256:
        raise ValueError("basis pilot verification hash mismatch")
    basis_report = json.loads(args.basis_pilot.read_text(encoding="utf-8"))
    basis_verification = json.loads(
        args.basis_verification.read_text(encoding="utf-8")
    )
    if (
        basis_report.get("status") != "COMPLETE"
        or basis_verification.get("status") != "PASS"
        or basis_verification["report"]["sha256"] != basis_hash
    ):
        raise ValueError("upstream exact pilot/verification gate failed")

    # Deliberately one core only.  No second cover is touched until this exact
    # construction and its independent symbolic checker are in place.
    cover = basis_report["covers"][0]
    graph_n = tuple(int(row) for row in cover["graph_n"])
    clique = tuple(int(vertex) for vertex in cover["fixed_maximum_clique"])
    clique_mask = sum(1 << vertex for vertex in clique)
    if len(clique) != 5 or int(cover["rank_upper"]) != 7:
        raise ValueError("first exact cover is not a K5/rank-seven instance")
    system = pentad.rank_two_system(graph_n, clique_mask)
    certificate = pentad.find_coefficientwise_certificate(graph_n, clique_mask)
    if certificate is None:
        raise ValueError("first core has no coefficientwise pentad certificate")
    positive = pentad.verify_coefficientwise_certificate(
        graph_n, clique_mask, certificate
    )

    degree_counts = Counter(
        pentad.polynomial_degree(equation)
        for equation in system.pentad_equations if equation
    )
    term_counts = [len(equation) for equation in system.pentad_equations]
    source_paths = [
        Path(__file__).resolve(),
        ROOT / "d6_k7_ranktwo_pentad.py",
        ROOT / "d6_k7_positive_polynomial_dual.py",
        ROOT / "d6_k7_rankone_tetrad_pilot.py",
    ]
    report = {
        "schema": 1,
        "kind": "d6_k7_ranktwo_pentad_first_core_pilot",
        "status": "REJECTED",
        "claim": (
            "For this one exact no-near cover, a required K5 leaves a PSD "
            "Schur complement of rank at most two. The stored signed pentad "
            "is identically zero for rank-two Gram data but expands to a "
            "nonzero polynomial with every coefficient positive in the "
            "strictly positive Sherman variables, an exact contradiction."
        ),
        "command": " ".join(sys.argv),
        "created_at": utc_now(),
        "elapsed_seconds": time.perf_counter() - started,
        "upstream": {
            "basis_pilot": {
                "path": str(args.basis_pilot), "sha256": basis_hash,
            },
            "basis_verification": {
                "path": str(args.basis_verification),
                "sha256": verification_hash,
                "status": basis_verification["status"],
            },
        },
        "source_hashes": {path.name: sha256(path) for path in source_paths},
        "git": git_provenance(),
        "machine": {
            "platform": platform.platform(),
            "python": sys.version,
            "implementation": platform.python_implementation(),
        },
        "scope": {
            "cores_tested": 1,
            "covers_tested": 1,
            "graph_index": int(cover["graph_index"]),
            "seed": cover["seed"],
            "zmask": int(cover["zmask"]),
            "nvertices": cover["nvertices"],
            "graph_n": cover["graph_n"],
            "rank_upper": int(cover["rank_upper"]),
            "clique": list(clique),
        },
        "system": {
            "remainder": list(system.remainder),
            "off_diagonal_pairs": len(system.pair_labels),
            "five_vertex_subsets": len(system.pentad_labels),
            "nonzero_pentads": sum(bool(item) for item in system.pentad_equations),
            "zero_pentads": sum(not item for item in system.pentad_equations),
            "pentad_degree_counts": {
                str(degree): count for degree, count in sorted(degree_counts.items())
            },
            "pentad_term_count_min": min(term_counts),
            "pentad_term_count_max": max(term_counts),
        },
        "certificate": certificate,
        "certificate_summary": {
            "pentad_index": int(certificate["pentad_index"]),
            "pentad_vertices": certificate["pentad_vertices"],
            "sign": int(certificate["sign"]),
            "positive_terms": len(positive),
            "polynomial_degree": pentad.polynomial_degree(positive),
            "coefficient_min": pentad.fraction_json(min(positive.values())),
            "coefficient_max": pentad.fraction_json(max(positive.values())),
            "monomial_degree_counts": {
                str(degree): count for degree, count in sorted(Counter(
                    sum(exponent) for exponent in positive
                ).items())
            },
        },
        "trust_scope": {
            "exact": (
                "Integer graph data and exact rational polynomial expansion."
            ),
            "not_claimed": (
                "This report tests one cover only; it does not by itself "
                "claim rejection of the other 52 covers or the full residue."
            ),
        },
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "status": report["status"],
        "certificate_summary": report["certificate_summary"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
