# Full K6 normal-inertia campaign

The exact normal-coordinate inertia rule proved in
`d6_k6_normal_coordinates.md` was run over every one of the 1,097 survivors
of the frozen K6 Lorentz/support/bipartite-rank campaign.  It rejects 107
graphs and leaves 990 K6-only candidates.

## Production result

The production source boundary was commit
`9b2c10ad4a89858a8b7ed942d282070218ab0b37` on branch
`codex/dimension6`.  The exact command was

```text
caffeinate -dimsu python3 d6_k6_normal_inertia_full_runner.py \
  --workers 11 --chunksize 1 --checkpoint-every 32
```

On macOS 14.5 arm64 with Python 3.11.15, the run completed in 0.670 seconds
without a resume or infrastructure error.

| quantity | exact value |
|---|---:|
| input graphs | 1,097 |
| graphs rejected | 107 |
| graphs surviving | 990 |
| K6 seeds checked | 35,925 |
| impossible K6 seeds | 210 |
| `Z0` subsets considered | 54,114 |
| support-matchable `Z0` subsets | 54,114 |
| bipartite components checked | 295,234 |
| failing `Z0` choices | 18,399 |

The complete 107-index rejection set is stored in the report and every
per-graph decision is stored in the deterministic gzip archive.  Decompressing
the archive gives exactly 1,097 newline-delimited JSON records with SHA-256

```text
c2f954dd46a8959cafbe40206e4cc7cc16a0df9598623d976ed9b09c1894dab2.
```

## Independent verification

The independent verifier does not import the production evaluator or runner.
It reconstructs the candidate validation, K6 clique enumeration, defect
masks, disjoint-defect graph, zero-factor subsets, support matching, and
bipartite components separately.  Instead of rational symmetric elimination,
it obtains every inertia from an exact SymPy integer characteristic
polynomial, square-free factorization, and exact Sturm root counts.

The command

```text
caffeinate -dimsu python3 d6_k6_normal_inertia_full_verifier.py \
  --workers 11 --chunksize 1
```

recomputed all 1,097 decisions in 1.366 seconds.  It matched:

- the complete 107-graph rejection set;
- all nine scalar decision fields for every graph;
- the first impossible K6 seed for every rejected graph;
- all archive, checkpoint, selection, input, and source hashes.

The known realizable 18-point construction independently passed all 32 of its
K6 seeds.  The verification status is `PASS`.

## Artifact hashes

```text
d6_k6_normal_inertia_full_runner.py
  a3ff3859686b503dd7c5066dcc77d07d15b57ea8eb169df28efc182ea5db5c70
d6_k6_normal_inertia_full_verifier.py
  7abd8d8c23d3a8f893eaf40aed42e1b5b953b7eb4392c82fcd9fd164ad16caa6
d6_k6_normal_inertia_full_test.py
  83cc285448fdc913304a4ffe5514469bc84ca26f72e5af0821734f2c8be22e15
d6_k6_normal_inertia_full_decisions.jsonl.gz
  614bbb64989b5834ed0dbe1c12e7718583e53cfd68483668ffbbad44608b3828
d6_k6_normal_inertia_full_checkpoint.json
  f4cb37faab90ed65b08b0e0ae2015a226de577caad5c74d23b530377d089962d
d6_k6_normal_inertia_full_report.json
  8fd3d9d10038040f42407b25ea007f59f0e84e6eb0da8bec2694d1831e1991c5
d6_k6_normal_inertia_full_verification_report.json
  61889185deb0e773bf899cdbb82c96467cc38a4688560a50f53c7761a7615e6b
```

The production report pins the frozen 1,098-graph input and earlier report;
the exact survivor selection removes only their already certified index
`461363`.  Candidate nonedges remain unconstrained throughout, and allowed
defect coordinates may be zero.
