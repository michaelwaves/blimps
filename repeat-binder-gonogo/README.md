# repeat-binder-gonogo

Validate a binder-scoring pipeline against repeat-protein complexes whose
outcomes are already known, and report whether it is fit to be trusted on novel
designs. The controls are then applied to a matched, de novo GvpA lattice-binder
experiment so every design claim stays inside the evidence they license.

Scorer under test: **Boltz-2**, run on Modal via proto-tools.

## GvpA design pipeline

`gvpa_rfd3.py` generates eight hotspot-constrained and eight free binder
backbones against a five-chain GvpA lattice patch. `gvpa_full_pipeline.py`
applies a geometry/clash referee, designs four ProteinMPNN sequences per
backbone, and screens distinct clash-free finalists with Boltz-2.

The matched experiment produced 64 sequences. Constrained backbones contact
6–10 nominated hotspot residues while free backbones contact none; the
constrained winner's three-seed mean binder-target ipTM is 0.2605 versus 0.0951
for the free control. Because G1 below fails on repeat proteins, these are
triage results and testable sequences—not validated binders or poses.

```bash
python gvpa_rfd3.py --arm constrained --n-designs 8
python gvpa_rfd3.py --arm free --n-designs 8
python gvpa_full_pipeline.py mpnn --n-per-backbone 4
python gvpa_full_pipeline.py compile
python gvpa_full_pipeline.py boltz --top-per-arm 2
```

The repository commits the 16 materialized RFdiffusion backbones, compact
candidate/result tables, and the browser demo. Raw Modal response payloads are
excluded. The new stages require `numpy`, `biotite`, and `proto-tools` in
addition to the benchmark dependencies below.

## Result

| Gate | | Verdict |
|---|---|---|
| **G5** | decoy rejection | **pass** — all 9 decoys below the weakest true positive (B2 = 66.75) |
| **G2** | specificity: DARPin 4m3 × mouse vs human cathepsin B | **pass** — margin 34.67 > band 12.08, Welch t = 6.08 |
| **G1** | interface recovery (DockQ vs deposit) | **fail** — 2/8 clear DockQ ≥ 0.49; median iRMSD 16.95 Å |
| **G3** | affinity ranking | **fail** — per-pocket rho ranges from 0.667 to −0.245 |
| **G4** | calibration | **pass** — AUC 0.775, but read it narrowly |

Bench verifier: **5/5**.

## The finding

**Boltz-2 fails on every genuine repeat protein in the set, and the case set's
own scaffold label hides it.**

| scaffold | DockQ ≥ 0.49 |
|---|---|
| non-repeat (affitin ×2) | **2/2** — 0.887, 0.819 |
| repeat (DARPin ×3, repebody, dArmRP, alphaRep) | **0/6** — 0.022 down to 0.009 |

Perfect separation, Fisher exact p = 0.036. `cases.yaml` declares
`scaffold_class: repeat proteins (ankyrin / armadillo / LRR / affitin /
alphaRep)`, but affitins are engineered from Sac7d — a single ~66-residue
Sulfolobus domain with no repeat architecture. Both successes are the two cases
that are not repeat proteins.

A mechanism consistent with the confidence data: a repeat protein presents a
*periodic* surface, so docking solutions are near-degenerate in register. The
model can build a correct solenoid and place it on the wrong repeat — confident
geometry, wrong answer. That is the observed signature, ipTM 0.59–0.66 against
DockQ ≈ 0.01 on C_1SVX, C_8RCI and C_9SPO. This is a hypothesis the run supports
but **did not test**; separating wrong-epitope from wrong-register needs a
per-repeat contact comparison that was not made.

The consequence for the stated purpose: this set exists to qualify a pipeline
for repeat-protein binder design. G2's pass says the scorer can tell two
83%-identical orthologs apart. It does not say the scorer can place a DARPin,
and G1 says it cannot.

## G2, and why the band was replaced

The first run reported G2 as **fail**, using a band it resolved silently in code
as the sample range (`max − min`) — 52.67 against a margin of 38.94.

Twelve seeds per case later (`results/g2_seed_depth.json`):

| n | margin | range band | | 95% CI half-width | | Welch t |
|---|---|---|---|---|---|---|
| 3 | 38.94 | 52.67 | fail | 66.24 | fail | 2.53 |
| 4 | 39.15 | 52.67 | fail | 36.47 | **pass** | 3.41 |
| 12 | 34.67 | 52.67 | fail | 12.14 | **pass** | 6.08 |

The range band is **frozen**. Both of B4's extremes landed in the first three
seeds, so nine further folds moved it by nothing; it would read 52.67 at n=100.
It is set by two draws and discards the rest. The CI half-width falls 66 → 12.

So the estimator failed the gate, not the scorer. The fix was **pre-registered
before it was applied** — `cases.yaml` now pins the band as a 95% CI half-width
on the difference of means, Welch, with a floor of 5 seeds per case below which
G2 returns `not_run` rather than `fail`. The superseded n=3 submission is
preserved unmodified at `results/predictions.n3-range-band.json`.

Choosing the estimator *is* choosing whether the gate is purchasable: under a CI
band this was answerable for ~$0.50, under the range band for no amount.

## Constraints — read these before trusting any number

`results/predictions.json` carries a `constraints` block. The load-bearing one:

> **ipTM is used as a proxy for interaction strength.** That holds for G2 and G5,
> where only ordering matters. It does **not** hold for G3 — ipTM is a structural
> confidence with no thermodynamic content, and the Tier A ladder spans 1–1000 nM
> (~4 kcal/mol) from single-residue substitutions in a 10-mer. The published
> baselines it is measured against compute actual energies. G3 tests that
> assumption; it does not measure affinity.

Also not modelled, each with why it matters: glycosylation, disulfide bonds (six
in mature cathepsin B — a plausible cause of the Tier C failures, undiagnosed),
pH (ground truth is pH-dependent: 65.7 nM at pH 6 vs 108.3 at pH 7), and
stoichiometry beyond 1:1 (4CJ2 is 2:2 in the crystal).

Sequence edits applied: purification tags stripped from 6 of 14 deposited binder
chains; mature chains throughout; decoys built by composition-preserving shuffle
and non-cognate pairing, excluding cross-reactive binders.

## Uncertainty

Every prediction carries mean, per-seed replicates, sample range and avg PAE —
no bare point estimates. Beyond that, the `uncertainty` block names where
confidence is *known to mislead*:

**G4's pass should be read narrowly.** AUC 0.775 scores confidence against
binder-*call* correctness on Tier B and the decoys. The DockQ numbers measure
something G4 never looked at, and there ipTM is uninformative. Confidence here
is informative about binding calls and not about structures.

Seeds were bought unevenly on purpose: 12 for B1/B4 where G2 turns on them, 3
for cases no gate turns on.

## Run it

Requires Python 3.12+, `pyyaml`, `pytest`; `DockQ` for G1; a Modal account with
`boltz2-prediction` deployed for the GPU phases.

```bash
python -m gonogo preflight   # free: what can run, what it costs, what blocks it
python -m gonogo inputs      # free: resolve 23 complexes, re-derive the case set's claims
python -m gonogo verify      # free: the bench's own verifier against the submission
```

`inputs` re-derives every factual claim in `cases.yaml` from primary sources —
the 83% identity, the I65/Q66 selectivity residues, the construct boundaries,
both Tier A pocket mappings — and prints pass/fail per claim. That phase caught
nearly every error in the case set, at zero cost.

GPU phases, billed to your own Modal account:

```bash
python -m gonogo --scorer boltz2 run      # G5, then G2, with the Referee's veto live
python g2_seed_depth.py --to 12           # extend B1/B4; reuses cached seeds
python g1_dockq.py                        # G1 via local DockQ; reuses cached folds
python finalize_submission.py             # assemble; re-runs no model
```

Harness self-test with **meaningless numbers**, kept off the submission path —
every prediction stamped `synthetic: true`, and the CLI refuses to write it into
the bench results directory:

```bash
python -m gonogo --scorer synthetic --out demo/selftest run --allow-synthetic
python -m gonogo --out demo/selftest verify    # 5 passed
```

## Layout

The harness resolves the task directory relative to its own package location
(`gonogo/../bench/repeat-binder-gonogo`), so **those two must keep their relative
positions**.

```
gonogo/                      the harness
  registry.py                cached, provenanced retrieval: PDB FASTA, UniProt,
                             author-numbered residues, purification-tag stripping
  caseset.py                 cases.yaml -> concrete two-chain complexes
  scorers.py                 Boltz2Scorer (Modal) + SyntheticScorer (self-test)
  gates.py                   G1-G5; reads G2's pre-registered band from cases.yaml
  cli.py                     phases, the Referee's veto, the submission writer
  trajectory.py              the Scribe: one JSONL line per step, dead ends included
  data/tier_a_kd.yaml        63 K_D values from PMC12202725 Table S1
  .cache/                    committed on purpose — lets `inputs` run offline

bench/repeat-binder-gonogo/  the task, unmodified except the G2 pre-registration
  results/predictions.json               the submission
  results/predictions.n3-range-band.json the superseded one, kept
  results/g2_seed_depth.json             the band diagnostic
  results/g1_dockq.json                  per-complex DockQ
  results/trajectory.jsonl               every step, including the dead ends

demo/selftest/               self-test output (synthetic, not results)
demo/full-sweep/             --ignore-veto sweep. NOT A VALID SUBMISSION: it
                             violates the task's stop-at-failed-G2 rule and says so
docs/gonogo-feasibility.md   what was checked, what broke, what it cost
g1_dockq.py, g2_seed_depth.py, finalize_submission.py
smoke_test.py                one representative tool per patched proto-tools app
```

## Provenance and dead ends

`results/trajectory.jsonl` keeps everything, including four runs aborted on a
Modal billing block, one deliberate veto override, and a Tier A sweep killed at
30/63 by an upstream connection reset. None of it is pruned — a record that only
shows the path that worked is not a record.

Measurements are line-cited in `bench/repeat-binder-gonogo/REFERENCES.md`:
Tier B from PMC12732864, Tier A and its baselines from PMC12202725.
