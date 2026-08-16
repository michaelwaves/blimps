# `gonogo` — the control run, in phases

A runnable harness for `bench/repeat-binder-gonogo`: validate a binder-scoring
pipeline against repeat-protein complexes whose outcomes are already known,
before trusting it on a novel design.

The feasibility assessment this produced is in `docs/gonogo-feasibility.md`.
The orchestration it implements is `docs/design-loop-sketch.md`.

## Run it

```bash
.venv/bin/python -m gonogo preflight     # free: what can run, what it costs, what blocks it
.venv/bin/python -m gonogo inputs        # free: resolve 23 complexes, re-derive the case set's claims
.venv/bin/python -m gonogo --scorer boltz2 run          # the real run — needs a deployed GPU tool
.venv/bin/python -m gonogo verify                       # the bench's own verifier
```

Harness self-test, with meaningless numbers, kept away from the submission path:

```bash
.venv/bin/python -m gonogo --scorer synthetic --out demo/selftest run --allow-synthetic
.venv/bin/python -m gonogo --out demo/selftest verify    # 5 passed
```

`preflight` exits non-zero when the scorer is unavailable. That is the check to
run before a run counts, because the bench verifier scores a `not_run` G2 as a
failure rather than an abstention.

## The four rules, and where each is enforced

| Rule from `task.md` | Enforced in |
|---|---|
| every score carries a confidence from the model that produced it | `scorers.Score` — there is no constructor that omits it |
| every prediction records tool and location | same |
| a gate that could not run is `not_run` with a reason | `gates.Verdict` — every branch returns evidence |
| one score convention, higher-is-better | `Boltz2Scorer.score_definition`, quoted into the submission notes |

## Layout

```
registry.py   cached, provenanced retrieval: PDB FASTA, UniProt features,
              author-numbered residues from mmCIF, terminal-tag stripping
caseset.py    cases.yaml -> concrete two-chain complexes; Tier A ladder
              construction; decoy construction that excludes cross-reactive pairs
scorers.py    Boltz2Scorer (real, Modal) and SyntheticScorer (self-test only)
gates.py      G1-G5 arithmetic; each returns pass / fail / not_run + evidence
cli.py        the phases, the Referee's veto, the submission writer
trajectory.py the Scribe: one JSONL line per step, dead ends included
data/         tier_a_kd.yaml (63 K_D values, line-cited) and baselines.yaml
```

## Adding a scorer

Subclass `Scorer`, set `key` / `ran_on` / `confidence_metric` /
`score_definition`, implement `preflight()` (raise `Unavailable` with a reason a
reader can act on) and `score()` (yield `Score` with replicates so a band
exists). Register it in `SCORERS`. Nothing else in the harness needs to change.

A scorer that cannot produce replicates cannot be used for G2 — the gate's
verdict is a margin measured against the model's own uncertainty band, and the
harness will report `not_run` rather than invent one.
