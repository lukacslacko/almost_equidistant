#!/usr/bin/env python3
"""Independent checker for near-clique active/sign certificates.

This checker intentionally does not import the certificate locator.  It
reconstructs the Venn regions, forced-active set, zero contradictions, and
signed-graph parity directly from the serialized abstract system.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from itertools import combinations
from pathlib import Path


POSITIVE = "FORCED_POSITIVE"
NEGATIVE = "FORCED_NEGATIVE"
ZERO = "IDENTICALLY_ZERO"
VARIABLE = "FLEXIBLE"


def _pair_kind(first: int, second: int, target: int, width: int) -> tuple[str, int | None]:
    universe = (1 << width) - 1
    common = first & second
    first_only = first & ~second & universe
    second_only = second & ~first & universe
    outside = universe & ~(first | second)
    if target == 0:
        if common == 0:
            return (POSITIVE, 1) if first_only and second_only else (ZERO, 0)
        if first_only == 0 or second_only == 0:
            return NEGATIVE, -1
        return VARIABLE, None
    if common == 0 or outside == 0:
        return POSITIVE, 1
    return VARIABLE, None


def _parse_system(
    certificate: dict,
) -> tuple[
    int,
    tuple[int, ...],
    tuple[int, ...],
    dict[tuple[int, int], int],
]:
    if certificate.get("schema") != 1:
        raise ValueError("unsupported certificate schema")
    if certificate.get("kind") != "rank_one_near_clique_active_sign_contradiction":
        raise ValueError("wrong certificate kind")
    width = int(certificate["variables"])
    if width < 0:
        raise ValueError("negative variable count")
    vertices = tuple(int(value) for value in certificate["vertices"])
    masks = tuple(int(value) for value in certificate["masks"])
    if len(vertices) != len(masks) or len(set(vertices)) != len(vertices):
        raise ValueError("malformed vertices or masks")
    if any(mask < 0 or mask >> width for mask in masks):
        raise ValueError("mask outside variable universe")
    targets: dict[tuple[int, int], int] = {}
    for item in certificate["targets"]:
        if not isinstance(item, list) or len(item) != 3:
            raise ValueError("malformed target record")
        left, right, target = map(int, item)
        if left >= right or target not in (0, 1):
            raise ValueError("target pair must be ordered and binary")
        if (left, right) in targets:
            raise ValueError("duplicate target pair")
        targets[left, right] = target
    expected = set(combinations(sorted(vertices), 2))
    if set(targets) != expected:
        raise ValueError("target pair set is incomplete or extraneous")
    return width, vertices, masks, targets


def _balanced(vertices: tuple[int, ...], forced: list[tuple[int, int, int]]) -> bool:
    neighbourhood = {vertex: [] for vertex in vertices}
    for left, right, sign in forced:
        neighbourhood[left].append((right, sign))
        neighbourhood[right].append((left, sign))
    assigned: dict[int, int] = {}
    for root in vertices:
        if root in assigned:
            continue
        assigned[root] = 1
        queue = deque([root])
        while queue:
            left = queue.popleft()
            for right, sign in neighbourhood[left]:
                required = assigned[left] * sign
                if right not in assigned:
                    assigned[right] = required
                    queue.append(right)
                elif assigned[right] != required:
                    return False
    return True


def verify_certificate(certificate: dict) -> dict:
    """Raise ``ValueError`` unless ``certificate`` proves the stated conflict."""

    width, vertices, masks, targets = _parse_system(certificate)
    mask_of = dict(zip(vertices, masks, strict=True))
    records: list[tuple[int, int, str, int | None]] = []
    for left, right in combinations(sorted(vertices), 2):
        kind, sign = _pair_kind(mask_of[left], mask_of[right], targets[left, right], width)
        records.append((left, right, kind, sign))

    profile_counter = Counter(kind for _, _, kind, _ in records)
    profile = {
        kind: profile_counter.get(kind, 0)
        for kind in (POSITIVE, NEGATIVE, ZERO, VARIABLE)
    }
    forced = [
        (left, right, int(sign))
        for left, right, _, sign in records
        if sign in (-1, 1)
    ]
    universe = (1 << width) - 1
    diagonal_forced_active = sorted(
        vertex
        for vertex, mask in zip(vertices, masks, strict=True)
        if mask in (0, universe)
    )
    forced_active = sorted(
        {vertex for left, right, _ in forced for vertex in (left, right)}
        | set(diagonal_forced_active)
    )
    active_set = set(forced_active)
    zero_conflicts = [
        (left, right)
        for left, right, kind, _ in records
        if kind == ZERO and left in active_set and right in active_set
    ]
    parity_conflict = not _balanced(tuple(forced_active), forced)
    reason = "ZERO_INSIDE_FORCED_ACTIVE" if zero_conflicts else (
        "FORCED_SIGN_PARITY" if parity_conflict else None
    )
    if reason is None:
        raise ValueError("serialized system has no active/sign contradiction")

    assessment = certificate.get("assessment")
    if not isinstance(assessment, dict) or assessment.get("rejected") is not True:
        raise ValueError("certificate assessment is not a rejection")
    if assessment.get("reason") != reason:
        raise ValueError("recorded rejection reason does not match reconstruction")
    if assessment.get("profile") != profile:
        raise ValueError("recorded classification profile does not match")
    if assessment.get("diagonal_forced_active") != diagonal_forced_active:
        raise ValueError("recorded diagonal-forced active set does not match")
    if assessment.get("forced_active") != forced_active:
        raise ValueError("recorded forced-active set does not match")
    if assessment.get("relaxed_active") is not None or assessment.get("relaxed_signs") is not None:
        raise ValueError("rejected certificate cannot contain a relaxed survivor")
    witness = assessment.get("witness")
    if not isinstance(witness, dict):
        raise ValueError("missing contradiction witness")

    if reason == "ZERO_INSIDE_FORCED_ACTIVE":
        zero_pair = witness.get("zero_pair")
        if zero_pair not in [list(pair) for pair in zero_conflicts]:
            raise ValueError("witnessed pair is not an identically-zero active pair")
        expected_forcing = [
            [left, right, sign]
            for left, right, sign in forced
            if zero_pair[0] in (left, right) or zero_pair[1] in (left, right)
        ]
        if witness.get("forcing_pairs") != expected_forcing:
            raise ValueError("zero witness forcing pairs do not match")
        expected_diagonal = [
            vertex for vertex in zero_pair
            if vertex in set(diagonal_forced_active)
        ]
        if witness.get("diagonal_forced_endpoints") != expected_diagonal:
            raise ValueError("zero witness diagonal activation does not match")
    else:
        if witness.get("forced_pairs") != [list(pair) for pair in forced]:
            raise ValueError("parity witness forced pairs do not match")
        conflict_edge = witness.get("conflict_edge")
        if not isinstance(conflict_edge, list) or len(conflict_edge) != 2:
            raise ValueError("malformed parity conflict edge")
        edge = tuple(sorted(map(int, conflict_edge)))
        sign_by_edge = {(left, right): sign for left, right, sign in forced}
        if edge not in sign_by_edge or int(witness.get("edge_sign")) != sign_by_edge[edge]:
            raise ValueError("parity conflict edge/sign is not forced")
        left_sign = int(witness.get("assigned_left_sign"))
        right_sign = int(witness.get("assigned_right_sign"))
        if left_sign not in (-1, 1) or right_sign not in (-1, 1):
            raise ValueError("invalid parity witness vertex sign")
        if left_sign * right_sign == int(witness["edge_sign"]):
            raise ValueError("recorded edge is not contradictory to recorded signs")

    return {
        "status": "VERIFIED",
        "reason": reason,
        "variables": width,
        "remainder_vertices": len(vertices),
        "diagonal_forced_active_vertices": len(diagonal_forced_active),
        "forced_active_vertices": len(forced_active),
        "forced_pairs": len(forced),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.certificate.read_text())
    if "certificate" in payload:
        payload = payload["certificate"]
    print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
