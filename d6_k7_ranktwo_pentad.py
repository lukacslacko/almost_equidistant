#!/usr/bin/env python3
"""Exact rank-two Schur pentads over a required K5.

If a normalized Gram matrix has rank at most seven and C is a required K5,
then the Schur complement after C is positive semidefinite of rank at most
two.  The known off-diagonal entries of any rank-two Gram matrix satisfy one
degree-five pentad on every five indices.  This module substitutes the exact
Sherman numerator polynomials into those pentads and recognizes the simplest
positive-polynomial certificate: one signed pentad whose every coefficient
is strictly positive.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping, Sequence

import d6_k7_positive_polynomial_dual as schur_dual
import d6_k7_rankone_tetrad_pilot as polynomial


Q = Fraction
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Q]
Edge = tuple[int, int]


# Local positions a,b,c,d,e are 0,1,2,3,4.  The signs and five edges are the
# exact expansion obtained by eliminating the first squared norm from two
# residual rank-one tetrads; see d6_k7_ranktwo_pentad.md.
PENTAD_TERMS: tuple[tuple[int, tuple[Edge, ...]], ...] = (
    (+1, ((0, 1), (0, 2), (1, 3), (2, 4), (3, 4))),
    (-1, ((0, 1), (0, 2), (1, 4), (2, 3), (3, 4))),
    (-1, ((0, 1), (0, 3), (1, 2), (2, 4), (3, 4))),
    (+1, ((0, 1), (0, 3), (1, 4), (2, 3), (2, 4))),
    (+1, ((0, 1), (0, 4), (1, 2), (2, 3), (3, 4))),
    (-1, ((0, 1), (0, 4), (1, 3), (2, 3), (2, 4))),
    (+1, ((0, 2), (0, 3), (1, 2), (1, 4), (3, 4))),
    (-1, ((0, 2), (0, 3), (1, 3), (1, 4), (2, 4))),
    (-1, ((0, 2), (0, 4), (1, 2), (1, 3), (3, 4))),
    (+1, ((0, 2), (0, 4), (1, 3), (1, 4), (2, 3))),
    (+1, ((0, 3), (0, 4), (1, 2), (1, 3), (2, 4))),
    (-1, ((0, 3), (0, 4), (1, 2), (1, 4), (2, 3))),
)


def _pair(first: int, second: int) -> Edge:
    return (first, second) if first < second else (second, first)


def pentad_polynomial(
    vertices: Sequence[int],
    pair_polynomials: Mapping[Edge, Polynomial],
    variables: int,
) -> Polynomial:
    """Expand the 12-term pentad for five ordered distinct vertices."""

    vertices = tuple(int(vertex) for vertex in vertices)
    if len(vertices) != 5 or len(set(vertices)) != 5:
        raise ValueError("a pentad needs five distinct vertices")
    one = {(0,) * variables: Q(1)}
    output: Polynomial = {}
    for sign, local_edges in PENTAD_TERMS:
        term = one
        for left, right in local_edges:
            edge = _pair(vertices[left], vertices[right])
            if edge not in pair_polynomials:
                raise ValueError(f"missing off-diagonal pair {edge}")
            term = polynomial.multiply(term, pair_polynomials[edge])
        for exponent, coefficient in term.items():
            polynomial.add_term(output, exponent, Q(sign) * coefficient)
    return output


@dataclass(frozen=True)
class RankTwoSystem:
    clique: tuple[int, ...]
    remainder: tuple[int, ...]
    masks: tuple[tuple[int, ...], ...]
    pair_labels: tuple[Edge, ...]
    pair_polynomials: tuple[Polynomial, ...]
    pentad_labels: tuple[tuple[int, ...], ...]
    pentad_equations: tuple[Polynomial, ...]


def rank_two_system(adj: Sequence[int], clique_mask: int) -> RankTwoSystem:
    """Construct all off-diagonal rank-two pentads after a required K5."""

    schur = schur_dual.schur_polynomial_system(adj, clique_mask)
    if len(schur.clique) != 5:
        raise ValueError("rank-two K7 pilot requires a five-vertex clique")
    # The saturating equation is f_yz=-g_yz; g_yz=T*S_yz.
    pair_polynomials = tuple(
        {exponent: -value for exponent, value in equation.items()}
        for equation in schur.equations
    )
    pair_by_label = dict(zip(schur.equation_pairs, pair_polynomials))
    labels = tuple(itertools.combinations(schur.remainder, 5))
    equations = tuple(
        pentad_polynomial(label, pair_by_label, len(schur.clique))
        for label in labels
    )
    return RankTwoSystem(
        clique=schur.clique,
        remainder=schur.remainder,
        masks=schur.masks,
        pair_labels=schur.equation_pairs,
        pair_polynomials=pair_polynomials,
        pentad_labels=labels,
        pentad_equations=equations,
    )


def fraction_json(value: Q) -> list[int]:
    return [value.numerator, value.denominator]


def parse_fraction(value: Sequence[int]) -> Q:
    if len(value) != 2:
        raise ValueError("fraction must have two entries")
    return Q(int(value[0]), int(value[1]))


def polynomial_degree(polynomial_value: Polynomial) -> int:
    return max(map(sum, polynomial_value), default=-1)


def verify_coefficientwise_certificate(
    adj: Sequence[int], clique_mask: int, certificate: dict
) -> Polynomial:
    """Re-expand one signed pentad and check strict coefficient positivity."""

    system = rank_two_system(adj, clique_mask)
    if certificate.get("schema") != 1:
        raise ValueError("unexpected certificate schema")
    if certificate.get("kind") != "rank_two_schur_pentad_coefficientwise_positive":
        raise ValueError("unexpected certificate kind")
    if certificate.get("clique") != list(system.clique):
        raise ValueError("certificate clique mismatch")
    if certificate.get("remainder") != list(system.remainder):
        raise ValueError("certificate remainder mismatch")
    if certificate.get("basis_neighbour_masks") != [
        list(mask) for mask in system.masks
    ]:
        raise ValueError("certificate basis masks mismatch")
    if certificate.get("pair_labels") != [
        list(pair) for pair in system.pair_labels
    ]:
        raise ValueError("certificate pair labels mismatch")
    index = int(certificate["pentad_index"])
    if not 0 <= index < len(system.pentad_equations):
        raise ValueError("pentad index out of range")
    if certificate.get("pentad_vertices") != list(system.pentad_labels[index]):
        raise ValueError("pentad vertex label mismatch")
    sign = int(certificate["sign"])
    if sign not in (-1, 1):
        raise ValueError("pentad sign must be plus or minus one")
    signed = {
        exponent: Q(sign) * value
        for exponent, value in system.pentad_equations[index].items()
    }
    if not signed or any(value <= 0 for value in signed.values()):
        raise ValueError("signed pentad is not coefficientwise strictly positive")
    serialized: Polynomial = {}
    for item in certificate.get("positive_polynomial", []):
        exponent = tuple(int(value) for value in item["monomial"])
        if len(exponent) != len(system.clique) or any(value < 0 for value in exponent):
            raise ValueError("invalid positive monomial")
        coefficient = parse_fraction(item["coefficient"])
        if coefficient <= 0 or exponent in serialized:
            raise ValueError("invalid/repeated positive coefficient")
        serialized[exponent] = coefficient
    if serialized != signed:
        raise ValueError("serialized positive polynomial mismatch")
    if int(certificate["polynomial_degree"]) != polynomial_degree(signed):
        raise ValueError("pentad polynomial degree mismatch")
    return signed


def find_coefficientwise_certificate(
    adj: Sequence[int], clique_mask: int
) -> dict | None:
    """Find a single pentad with one-sign nonzero coefficients."""

    system = rank_two_system(adj, clique_mask)
    for index, equation in enumerate(system.pentad_equations):
        if not equation:
            continue
        values = tuple(equation.values())
        if all(value > 0 for value in values):
            sign = 1
        elif all(value < 0 for value in values):
            sign = -1
        else:
            continue
        signed = {
            exponent: Q(sign) * value for exponent, value in equation.items()
        }
        certificate = {
            "schema": 1,
            "kind": "rank_two_schur_pentad_coefficientwise_positive",
            "clique": list(system.clique),
            "remainder": list(system.remainder),
            "basis_neighbour_masks": [list(mask) for mask in system.masks],
            "pair_labels": [list(pair) for pair in system.pair_labels],
            "pentad_index": index,
            "pentad_vertices": list(system.pentad_labels[index]),
            "sign": sign,
            "polynomial_degree": polynomial_degree(signed),
            "positive_polynomial": [
                {
                    "monomial": list(exponent),
                    "coefficient": fraction_json(value),
                }
                for exponent, value in sorted(signed.items())
            ],
        }
        verify_coefficientwise_certificate(adj, clique_mask, certificate)
        return certificate
    return None
