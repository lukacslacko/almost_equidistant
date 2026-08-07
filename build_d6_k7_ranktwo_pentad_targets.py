#!/usr/bin/env python3
"""Build the exact 53-cover target manifest for the full pentad scan."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import run_d6_k7_arbitrary_basis_extension_pilot as loader


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rank-input", type=Path,
        default=Path(".runs/d6_k7_rank_survivors.json"),
    )
    parser.add_argument(
        "--selection", type=Path,
        default=Path("d6_k7_rankone_tetrad_full_selection.json"),
    )
    parser.add_argument(
        "--union", type=Path,
        default=Path("d6_k7_rankone_pattern_union.json"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_ranktwo_pentad_targets.json"),
    )
    args = parser.parse_args()
    graph_by_index, graph_indices, witnesses, provenance = (
        loader.load_campaign_inputs(
            args.rank_input, args.selection, args.union
        )
    )
    targets = loader.extract_no_near_covers(
        graph_by_index, graph_indices, witnesses
    )
    if len(targets) != 53 or len({item["graph_index"] for item in targets}) != 50:
        raise ValueError("target manifest is not 53 covers on 50 graphs")
    if any(
        len(item["graph_n"]) != 12
        or int(item["rank_upper"]) != 7
        or int(item["maximum_clique_size"]) != 5
        or int(item["zmask"]) != 0
        for item in targets
    ):
        raise ValueError("target manifest does not share n=12/U=7/omega=5/z=0")
    keys = [
        [item["graph_index"], item["seed"], item["zmask"]]
        for item in targets
    ]
    payload = {
        "schema": 1,
        "kind": "d6_k7_ranktwo_pentad_exact_target_manifest",
        "status": "COMPLETE",
        "claim": (
            "This is the exact replayed list of prior-passing K7 covers in "
            "the 258-graph union residue that have neither a saturating nor "
            "near-saturating clique. It makes no realizability claim."
        ),
        "provenance": provenance,
        "source_hashes": {
            Path(__file__).name: sha256(Path(__file__).resolve()),
            "run_d6_k7_arbitrary_basis_extension_pilot.py": sha256(
                Path(__file__).resolve().parent
                / "run_d6_k7_arbitrary_basis_extension_pilot.py"
            ),
        },
        "population": {
            "graphs": len(graph_indices),
            "covers": len(targets),
            "cover_keys_sha256": stable_hash(keys),
            "shared_parameters": {
                "n": 12, "rank_upper": 7, "maximum_clique": 5,
                "remainder": 7, "pentads_per_cover": 21, "zmask": 0,
            },
        },
        "targets": targets,
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "graphs": len(graph_indices),
        "covers": len(targets),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
