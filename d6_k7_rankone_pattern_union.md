# Exact K7 rank-one-tetrad / pattern-954 union

## Certified result

The independently verified degree-four rank-one tetrad campaign and the
independently edge-checked pattern-954 containment campaign act on the same
ordered 12,839-graph K7 selection.  Their exact overlap is:

| set | graphs |
|---|---:|
| tetrad rejected | 11,902 |
| pattern-954 rejected | 3,403 |
| rejected by both | 2,724 |
| tetrad only | 9,178 |
| pattern-954 only | 679 |
| exact union | **12,581** |
| ordered complement | **258** |

The exact K7 residue index hash is

```text
55ab329dbe3ae4dbfa10378ff023a168fc6ab306ffc7af14bbd9790a68108a09
```

Residue membership is not a realizability claim.  It means only that neither
of these two exact certificates rejects the graph.

## Verification gates

The union builder refuses to run unless all of the following hashes and
complete-status checks agree:

```text
frozen selection
  86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479
tetrad production report
  ae2075ac83abcdfc0b7d42f63b9515c4e48b40cf9c976178aa366d151a78996a
tetrad independent verification
  1812524c835fd6635b9815c0f0c25e9a49d99ff5d3b042e4b194ac5d7898dc97
pattern-954 production report
  a2d4fb5a07dcea580315730b4155fc0c08ecd0439a9111e733d29e2da0b5be62
pattern-954 independent verification
  a3a353e4d6cab462697f9565b69845f4a8f296fd6f8085ed77d68897257eab32
```

The tetrad verifier replayed all 12,839 graph quantifiers, 1,131 inherited
degree-one identities, and 35,160 new tetrad identities over exact rational
arithmetic.  Its status counts are `PASS: 12,839` with an empty error list.
The pattern verifier rechecked all 3,403 positive mappings for injectivity
and all 54 required edges.  Pattern nonedges impose no distance conditions.

The union builder then parses all 12,839 rows of the compressed tetrad
decision archive, checks its compressed and uncompressed hashes, reconstructs
the production totals, and performs set arithmetic in frozen selection order.
The independent union checker imports neither the builder nor either search
runner and reconstructs the same sets and cover profiles.

Final union artifacts:

```text
d6_k7_rankone_pattern_union.json
  1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599
d6_k7_rankone_pattern_union_verification.json
  7b9d5e35bf787e5d05aa486f914416598c31836264a479fbbc1af117f1dfcc14
```

The independent result is `PASS`, checks all 12,839 decision rows and all
258 residue profiles, and obtains the same residue hash shown above.

## Structure of the exact residue

The production tetrad layer applied rank-one equations only to prior-passing
covers without a saturating clique and with a near-saturating clique.  The
258 remaining graphs split as follows:

| aggregate cover structure | graphs |
|---|---:|
| only saturating prior-passing covers | 23 |
| both saturating and non-saturating covers | 157 |
| only non-saturating covers | 78 |

Obstacle flags can overlap:

| surviving mechanism present | graphs |
|---|---:|
| at least one saturating cover | 180 |
| at least one cover with no near-saturating clique | 50 |
| at least one tetrad-resistant near-clique cover | 135 |

Every residue index and its aggregate cover counts are explicit in the union
manifest.  These figures identify the next bounded pilots: saturating-clique
special-inverse algebra, maximum-clique extensions for the 50 no-near cases,
and stronger rank-one consequences for the 135 tetrad-resistant cases.

## Deterministic reproduction

```text
python3 build_d6_k7_rankone_pattern_union.py \
  --tetrad-report-sha256 \
    ae2075ac83abcdfc0b7d42f63b9515c4e48b40cf9c976178aa366d151a78996a \
  --tetrad-verification-sha256 \
    1812524c835fd6635b9815c0f0c25e9a49d99ff5d3b042e4b194ac5d7898dc97

python3 verify_d6_k7_rankone_pattern_union.py \
  --manifest-sha256 \
    1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599

python3 -m unittest -v test_d6_k7_rankone_pattern_union.py
```

The rebuild is byte-identical, the independent checker returns `PASS`, and
all seven focused tests pass.

## Trust scope

The mathematical rejection union trusts the two independently verified
source certificates, exact JSON/TSV parsing, SHA-256 binding, Python
arbitrary-precision integer and `Fraction` arithmetic, and the documented
Euclidean reductions behind the source campaigns.  The numerical LP locator
used to discover tetrad identities is outside the trust boundary.  No
floating-point rank or failed search is used as a rejection.
