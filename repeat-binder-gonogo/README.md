# repeat-binder-gonogo — a control run that stopped itself

Validate a binder-scoring pipeline against repeat-protein complexes whose
outcomes are already known, and report whether it is fit to be trusted on novel
designs. Nothing here is a design — these are controls, chosen because they are
the cheapest systems on which a scoring model can be caught being wrong.

**Result of the run in `bench/repeat-binder-gonogo/results/`:** the scorer
(Boltz-2, on Modal) passed the decoy gate and failed the specificity gate. The
run stopped there rather than spending on the remaining tiers.

| Gate | | Verdict |
|---|---|---|
| **G5** | decoy rejection | **pass** — all 9 decoys below the weakest true positive (B2 = 66.748) |
| **G2** | specificity: DARPin 4m3 vs mouse/human cathepsin B | **fail** — margin 38.941 ≤ band 52.671 |
| G1 | interface recovery | `not_run` — vetoed by G2 |
| G3 | affinity ranking | `not_run` — vetoed by G2 |
| G4 | calibration | `not_run` — vetoed by G2 |

Bench verifier: **4 of 5 checks pass**, reward **0.0** — `verifier/test.sh` caps
the reward at zero on any G2 failure, by design.

## Why G2 failed, in one paragraph

B1 (DARPin 4m3 × *mouse* cathepsin B, a 65.7 nM binder) scored 78.688. B4 (the
same DARPin × *human* cathepsin B, 83% identical, no binding at 1000× molar
excess) scored 39.747. The ranking is correct. The problem is reproducibility:
B4's three seeds came back **14.618, 37.333, 67.289** — one of them scores the
true non-binder above B1's own weakest replicate. The gate compares the margin
against the model's own seed spread and refused to certify specificity on a
signal that unstable.

That is the harness working as intended. A pipeline that cannot tell these two
apart *repeatably* cannot support a specificity claim on a novel design, and the
point of a go/no-go set is to find that out for $3.66 instead of downstream.

⚠️ **Known limitation, documented rather than patched:** the band is the sample
range (`max − min`), whose expectation *grows* with the number of seeds — so
more evidence can never un-stick this gate. Two of three defensible band
definitions still say fail, so the verdict stands, but the estimator should be
re-declared before the next run. See §4 of
[`docs/gonogo-feasibility.md`](docs/gonogo-feasibility.md).

## Run it

Requires Python 3.12+, `pyyaml`, `pytest`, and — for the real scorer — a Modal
account with `boltz2-prediction` deployed via proto-tools.

```bash
python -m gonogo preflight   # free: what can run, what it costs, what blocks it
python -m gonogo inputs      # free: resolve 23 complexes, re-derive the case set's claims
python -m gonogo verify      # free: the bench's own verifier against the submission
```

The two free phases are worth running first — `inputs` re-derives every factual
claim in `cases.yaml` from primary sources (the 83% identity, the I65/Q66
selectivity residues, the construct boundaries, both Tier A pocket mappings) and
prints pass/fail per claim. That is where nearly every error was caught, at zero
cost.

The real run needs a GPU and bills your own Modal account:

```bash
python -m gonogo --scorer boltz2 run
```

Harness self-test, with **meaningless numbers**, kept away from the submission
path — every prediction is stamped `synthetic: true` and the CLI refuses to
write it into the bench results dir:

```bash
python -m gonogo --scorer synthetic --out demo/selftest run --allow-synthetic
python -m gonogo --out demo/selftest verify    # 5 passed
```

## Layout

The harness resolves the task directory relative to its own package location
(`gonogo/../bench/repeat-binder-gonogo`), so **these two must keep their relative
positions** or nothing resolves:

```
gonogo/                      the harness
  registry.py                cached, provenanced retrieval: PDB FASTA, UniProt
                             features, author-numbered residues, tag stripping
  caseset.py                 cases.yaml -> concrete two-chain complexes
  scorers.py                 Boltz2Scorer (real, Modal) + SyntheticScorer (self-test)
  gates.py                   G1-G5 arithmetic; each returns pass/fail/not_run + evidence
  cli.py                     the phases, the Referee's veto, the submission writer
  trajectory.py              the Scribe: one JSONL line per step, dead ends included
  data/tier_a_kd.yaml        63 K_D values extracted from PMC12202725 Table S1
  .cache/                    committed on purpose — lets `inputs` run offline

bench/repeat-binder-gonogo/  the task, unmodified
  task.md, cases.yaml, REFERENCES.md
  verifier/                  the bench's own pytest suite + reward script
  results/                   the submission, the trajectory, the per-complex score cache

demo/selftest/               self-test output (synthetic numbers, not results)
docs/gonogo-feasibility.md   what was checked, what broke, what it cost
smoke_test.py                one representative tool per patched proto-tools app
```

## Provenance

Every score in `results/predictions.json` carries its confidence (ipTM), the
tool that produced it (`boltz2-prediction`), where it ran (`modal`), its
per-seed replicates and its avg PAE. Every gate carries evidence, and a gate
that could not run says so with a reason rather than being omitted.
`results/trajectory.jsonl` logs one line per step including the dead ends — the
four earlier attempts that aborted on a Modal billing block are still in there.

Measurements are line-cited in `bench/repeat-binder-gonogo/REFERENCES.md`:
Tier B from PMC12732864, Tier A and its physics baselines from PMC12202725.
