#!/usr/bin/env python3
"""Exact discrete active/sign consequences of a rank-one Schur complement.

This module deliberately stops before polynomial feasibility.  Given a
required clique of size one below the normalized-Gram rank bound, it classifies
the sign of every known off-diagonal Schur numerator using only which positive
variables occur in the four Venn regions of the two clique-neighbour masks.

The output is either a small, exact contradiction certificate or a labelled
survivor of this *relaxed* discrete layer.  A survivor is not a realization.
In particular, ``FLEXIBLE`` means sign-variable and zero-capable; it is never
silently treated as nonzero.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Mapping, Sequence


FORCED_POSITIVE = "FORCED_POSITIVE"
FORCED_NEGATIVE = "FORCED_NEGATIVE"
IDENTICALLY_ZERO = "IDENTICALLY_ZERO"
FLEXIBLE = "FLEXIBLE"


@dataclass(frozen=True)
class PairClassification:
    """Exact qualitative classification of one Schur numerator."""

    left: int
    right: int
    target: int
    intersection: int
    left_only: int
    right_only: int
    outside_union: int
    classification: str
    forced_sign: int | None


@dataclass(frozen=True)
class ActiveSignAssessment:
    """Result for one near-clique system.

    ``relaxed_active`` and ``relaxed_signs`` describe one assignment satisfying
    just the discrete sign/zero abstraction.  They make no claim about a common
    positive choice of the clique variables or the diagonal Schur entries.
    """

    rejected: bool
    reason: str | None
    profile: dict[str, int]
    diagonal_forced_active: tuple[int, ...]
    forced_active: tuple[int, ...]
    relaxed_active: tuple[int, ...] | None
    relaxed_signs: tuple[tuple[int, int], ...] | None
    witness: dict | None


def _validate_mask(mask: int, variables: int, name: str) -> None:
    if variables < 0:
        raise ValueError("variables must be nonnegative")
    if mask < 0 or mask >> variables:
        raise ValueError(f"{name} has a bit outside the {variables} variables")


def classify_pair(
    left_mask: int,
    right_mask: int,
    target: int,
    variables: int,
    *,
    left: int = 0,
    right: int = 1,
) -> PairClassification:
    """Classify ``g=T(k-p_left^T H p_right)`` for positive variables.

    ``target`` is the known normalized-Gram off-diagonal entry: one for a
    required edge and zero for a nonedge of the post-cover exact support graph.
    The latter zero is valid only after the K7/zero-factor cover quantifiers
    have established it; this routine does not reinterpret original candidate
    nonedges as prescribed non-unit distances.
    """

    _validate_mask(left_mask, variables, "left mask")
    _validate_mask(right_mask, variables, "right mask")
    if target not in (0, 1):
        raise ValueError("target must be zero or one")
    if left == right:
        raise ValueError("pair labels must be distinct")

    full = (1 << variables) - 1
    intersection = left_mask & right_mask
    left_only = left_mask & ~right_mask & full
    right_only = right_mask & ~left_mask & full
    outside_union = full & ~(left_mask | right_mask)

    if target == 0:
        if not intersection:
            if left_only and right_only:
                kind, sign = FORCED_POSITIVE, 1
            else:
                kind, sign = IDENTICALLY_ZERO, 0
        elif not left_only or not right_only:
            kind, sign = FORCED_NEGATIVE, -1
        else:
            kind, sign = FLEXIBLE, None
    elif not intersection or not outside_union:
        kind, sign = FORCED_POSITIVE, 1
    else:
        kind, sign = FLEXIBLE, None

    return PairClassification(
        left=left,
        right=right,
        target=target,
        intersection=intersection,
        left_only=left_only,
        right_only=right_only,
        outside_union=outside_union,
        classification=kind,
        forced_sign=sign,
    )


def _normalise_targets(
    vertices: Sequence[int], targets: Mapping[tuple[int, int], int]
) -> dict[tuple[int, int], int]:
    labels = tuple(int(vertex) for vertex in vertices)
    if len(set(labels)) != len(labels):
        raise ValueError("vertex labels must be distinct")
    expected = {tuple(sorted(pair)) for pair in combinations(labels, 2)}
    normalised: dict[tuple[int, int], int] = {}
    for raw_pair, raw_target in targets.items():
        if len(raw_pair) != 2:
            raise ValueError("target keys must be vertex pairs")
        left, right = map(int, raw_pair)
        if left == right:
            raise ValueError("target keys cannot be loops")
        pair = tuple(sorted((left, right)))
        if pair in normalised:
            raise ValueError(f"duplicate target pair {pair}")
        target = int(raw_target)
        if target not in (0, 1):
            raise ValueError("targets must be zero or one")
        normalised[pair] = target
    if set(normalised) != expected:
        missing = sorted(expected - set(normalised))
        extra = sorted(set(normalised) - expected)
        raise ValueError(f"target pairs mismatch: missing={missing}, extra={extra}")
    return normalised


def classify_system(
    variables: int,
    vertices: Sequence[int],
    masks: Sequence[int],
    targets: Mapping[tuple[int, int], int],
) -> tuple[PairClassification, ...]:
    """Classify all pairs of one abstract near-clique system."""

    labels = tuple(int(vertex) for vertex in vertices)
    if len(labels) != len(masks):
        raise ValueError("one mask is required for every remainder vertex")
    mask_by_vertex = dict(zip(labels, map(int, masks), strict=True))
    for vertex, mask in mask_by_vertex.items():
        _validate_mask(mask, variables, f"mask for vertex {vertex}")
    normalised = _normalise_targets(labels, targets)
    ordered = tuple(sorted(labels))
    return tuple(
        classify_pair(
            mask_by_vertex[left],
            mask_by_vertex[right],
            normalised[(left, right)],
            variables,
            left=left,
            right=right,
        )
        for left, right in combinations(ordered, 2)
    )


def _signed_colouring(
    vertices: Sequence[int],
    forced_pairs: Sequence[PairClassification],
) -> tuple[dict[int, int] | None, dict | None]:
    adjacency: dict[int, list[tuple[int, int]]] = {
        vertex: [] for vertex in vertices
    }
    for pair in forced_pairs:
        assert pair.forced_sign in (-1, 1)
        adjacency[pair.left].append((pair.right, pair.forced_sign))
        adjacency[pair.right].append((pair.left, pair.forced_sign))
    for neighbours in adjacency.values():
        neighbours.sort()

    signs: dict[int, int] = {}
    for root in sorted(vertices):
        if root in signs:
            continue
        signs[root] = 1
        queue = deque([root])
        while queue:
            left = queue.popleft()
            for right, edge_sign in adjacency[left]:
                required = signs[left] * edge_sign
                if right not in signs:
                    signs[right] = required
                    queue.append(right)
                elif signs[right] != required:
                    return None, {
                        "conflict_edge": [left, right],
                        "edge_sign": edge_sign,
                        "assigned_left_sign": signs[left],
                        "assigned_right_sign": signs[right],
                    }
    return signs, None


def assess_system(
    variables: int,
    vertices: Sequence[int],
    masks: Sequence[int],
    targets: Mapping[tuple[int, int], int],
) -> ActiveSignAssessment:
    """Apply the complete pairwise active/sign abstraction.

    Every forced-nonzero pair forces both endpoints active.  The diagonal
    Schur condition also forces vertices with empty or full masks active:
    inactivity would require ``T*w(A)-w(A)^2-T>0``, which is respectively
    ``-T`` or ``-1`` for those two masks.  An identically zero pair between
    forced-active endpoints is impossible.  Otherwise the remaining forced
    signs must be a balanced signed graph.
    """

    pairs = classify_system(variables, vertices, masks, targets)
    profile = dict(Counter(pair.classification for pair in pairs))
    profile = {
        kind: profile.get(kind, 0)
        for kind in (
            FORCED_POSITIVE,
            FORCED_NEGATIVE,
            IDENTICALLY_ZERO,
            FLEXIBLE,
        )
    }
    forced_pairs = tuple(
        pair for pair in pairs if pair.forced_sign in (-1, 1)
    )
    full = (1 << variables) - 1
    diagonal_forced_active = tuple(sorted(
        int(vertex)
        for vertex, mask in zip(vertices, masks, strict=True)
        if int(mask) in (0, full)
    ))
    forced_active = tuple(sorted({
        endpoint
        for pair in forced_pairs
        for endpoint in (pair.left, pair.right)
    } | set(diagonal_forced_active)))
    active_set = set(forced_active)

    zero_conflict = next(
        (
            pair for pair in pairs
            if pair.classification == IDENTICALLY_ZERO
            and pair.left in active_set
            and pair.right in active_set
        ),
        None,
    )
    if zero_conflict is not None:
        return ActiveSignAssessment(
            rejected=True,
            reason="ZERO_INSIDE_FORCED_ACTIVE",
            profile=profile,
            diagonal_forced_active=diagonal_forced_active,
            forced_active=forced_active,
            relaxed_active=None,
            relaxed_signs=None,
            witness={
                "zero_pair": [zero_conflict.left, zero_conflict.right],
                "forcing_pairs": [
                    [pair.left, pair.right, pair.forced_sign]
                    for pair in forced_pairs
                    if zero_conflict.left in (pair.left, pair.right)
                    or zero_conflict.right in (pair.left, pair.right)
                ],
                "diagonal_forced_endpoints": [
                    vertex
                    for vertex in (zero_conflict.left, zero_conflict.right)
                    if vertex in set(diagonal_forced_active)
                ],
            },
        )

    signs, conflict = _signed_colouring(forced_active, forced_pairs)
    if conflict is not None:
        return ActiveSignAssessment(
            rejected=True,
            reason="FORCED_SIGN_PARITY",
            profile=profile,
            diagonal_forced_active=diagonal_forced_active,
            forced_active=forced_active,
            relaxed_active=None,
            relaxed_signs=None,
            witness={
                **conflict,
                "forced_pairs": [
                    [pair.left, pair.right, pair.forced_sign]
                    for pair in forced_pairs
                ],
            },
        )

    assert signs is not None
    return ActiveSignAssessment(
        rejected=False,
        reason=None,
        profile=profile,
        diagonal_forced_active=diagonal_forced_active,
        forced_active=forced_active,
        relaxed_active=forced_active,
        relaxed_signs=tuple(sorted(signs.items())),
        witness=None,
    )


def branch_is_pairwise_compatible(
    pairs: Sequence[PairClassification],
    active: Sequence[int],
    signs: Mapping[int, int] | None = None,
    required_active: Sequence[int] = (),
) -> bool:
    """Check one active/sign branch in the deliberately relaxed abstraction.

    This helper is useful for controls and downstream branch generation.  A
    flexible pair is accepted both at zero and at either nonzero sign because
    the pairwise abstraction does not assert simultaneous realizability.
    """

    active_set = set(map(int, active))
    if not set(map(int, required_active)) <= active_set:
        return False
    sign_map = {} if signs is None else {int(k): int(v) for k, v in signs.items()}
    if set(sign_map) - active_set or any(value not in (-1, 1) for value in sign_map.values()):
        return False
    for pair in pairs:
        both = pair.left in active_set and pair.right in active_set
        if not both:
            if pair.forced_sign in (-1, 1):
                return False
            continue
        if pair.classification == IDENTICALLY_ZERO:
            return False
        if pair.forced_sign in (-1, 1):
            if pair.left not in sign_map or pair.right not in sign_map:
                return False
            if sign_map[pair.left] * sign_map[pair.right] != pair.forced_sign:
                return False
    return True


def graph_system(
    adj: Sequence[int], clique_mask: int
) -> tuple[int, tuple[int, ...], tuple[int, ...], dict[tuple[int, int], int]]:
    """Extract masks and exact zero/one targets from a post-cover graph."""

    nvertices = len(adj)
    full = (1 << nvertices) - 1
    for vertex, row in enumerate(adj):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError("invalid adjacency row")
        if any(
            bool(row & (1 << other)) != bool(adj[other] & (1 << vertex))
            for other in range(vertex)
        ):
            raise ValueError("asymmetric adjacency")
    clique = tuple(vertex for vertex in range(nvertices) if clique_mask & (1 << vertex))
    if not clique or clique_mask >> nvertices:
        raise ValueError("invalid or empty clique mask")
    if any(not (adj[left] & (1 << right)) for left, right in combinations(clique, 2)):
        raise ValueError("clique mask is not a required clique")
    remainder = tuple(vertex for vertex in range(nvertices) if not (clique_mask & (1 << vertex)))
    masks = tuple(
        sum(
            1 << coordinate
            for coordinate, basis in enumerate(clique)
            if adj[vertex] & (1 << basis)
        )
        for vertex in remainder
    )
    targets = {
        (left, right): int(bool(adj[left] & (1 << right)))
        for left, right in combinations(remainder, 2)
    }
    return len(clique), remainder, masks, targets


def assess_graph(adj: Sequence[int], clique_mask: int) -> ActiveSignAssessment:
    """Assess one near clique in an already justified exact support graph."""

    variables, vertices, masks, targets = graph_system(adj, clique_mask)
    return assess_system(variables, vertices, masks, targets)


def make_certificate(
    variables: int,
    vertices: Sequence[int],
    masks: Sequence[int],
    targets: Mapping[tuple[int, int], int],
    assessment: ActiveSignAssessment,
) -> dict:
    """Serialize a rejected system for the independent verifier."""

    if not assessment.rejected:
        raise ValueError("only contradictions are certificates")
    normalised = _normalise_targets(vertices, targets)
    return {
        "schema": 1,
        "kind": "rank_one_near_clique_active_sign_contradiction",
        "variables": int(variables),
        "vertices": list(map(int, vertices)),
        "masks": list(map(int, masks)),
        "targets": [
            [left, right, target]
            for (left, right), target in sorted(normalised.items())
        ],
        # Force tuple-valued dataclass fields into their on-disk JSON shape so
        # the in-memory certificate and its serialized form verify identically.
        "assessment": json.loads(json.dumps(asdict(assessment))),
    }


def _load_system(path: Path) -> tuple[int, list[int], list[int], dict[tuple[int, int], int]]:
    payload = json.loads(path.read_text())
    return (
        int(payload["variables"]),
        [int(value) for value in payload["vertices"]],
        [int(value) for value in payload["masks"]],
        {
            (int(left), int(right)): int(target)
            for left, right, target in payload["targets"]
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("system", type=Path, help="JSON abstract system")
    parser.add_argument("--certificate", type=Path)
    args = parser.parse_args()
    variables, vertices, masks, targets = _load_system(args.system)
    assessment = assess_system(variables, vertices, masks, targets)
    output: dict = {"assessment": asdict(assessment)}
    if assessment.rejected:
        output["certificate"] = make_certificate(
            variables, vertices, masks, targets, assessment
        )
    encoded = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.certificate is None:
        print(encoded, end="")
    else:
        args.certificate.write_text(encoded)


if __name__ == "__main__":
    main()
