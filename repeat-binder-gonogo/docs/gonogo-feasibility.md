# Is the design loop feasible? Checked against `bench/repeat-binder-gonogo/task.md`

Companion to `design-loop-sketch.md`. That document proposes how agents compose
design loops; this one runs the proposal against the one task in the repo that
has a verifier attached, and reports where it holds and where it breaks.

Everything below was executed, not reasoned about. The harness is
`gonogo/`, runnable now:

```
.venv/bin/python -m gonogo preflight                 # what can run, and what it costs
.venv/bin/python -m gonogo inputs                    # resolve cases, re-derive their claims
.venv/bin/python -m gonogo --scorer boltz2 run       # the real run, gates in cost order
.venv/bin/python -m gonogo verify                    # the bench's own verifier
```

---

## 1. Headline

**The loop's shape survives contact with the task. The run completed. It was
the scorer that failed the gate, not the harness that failed to reach it.**

The GPU block described in earlier drafts of this document is resolved:
`boltz2-prediction` is deployed to `proto-env` in workspace `ghosos`, and
`gonogo preflight` reports the scorer AVAILABLE. The control run executed on
2026-08-15 at 20:00, 15 complexes at three seeds each, ~$3.66:

| Gate | Status | Why |
|---|---|---|
| **G5** decoy rejection | **pass** | all 9 decoys scored below the weakest true positive (B2 = 66.748) |
| **G2** specificity | **fail** | margin 38.941 ≤ band 52.671 |
| G1 / G3 / G4 | `not_run` | the Referee's veto — the cost policy stops the run when G2 fails |

The verifier scores this **0.0**, and that number needs reading carefully. Four
of its five checks pass; the fifth is `test_g2_separates_mouse_from_human_cathepsin_b`,
and `test.sh` caps the whole run's reward at zero on any G2 failure. A 0.0 here
does not mean "nothing worked" — it means the one gate the exercise exists for
did not clear, which is exactly the signal the cap is designed to send.

Stated plainly: Boltz-2 ranks the true binder above the true non-binder
(78.688 vs 39.747), but it does not do so *reliably*. B4's three seeds are
14.618, 37.333 and 67.289 — one of them scores the human ortholog, a case with
no binding at 1000× molar excess, above B1's own weakest replicate. The gate
refused to certify specificity on a signal that unstable, which is the correct
verdict and the useful one (§4).

---

## 2. What was verified on CPU, for free, before any budget question

The whole input-resolution phase runs in-process through proto-tools and RCSB
with no billing. It re-derives the case set's own claims instead of trusting
them. All five checks pass:

| Claim in `cases.yaml` | Re-derived | Result |
|---|---|---|
| mouse / human cathepsin B are 83% identical | global alignment of both mature chains (UniProt P10605 / P07858, residues 80–333) | **83.07%** |
| selectivity residues mouse I65/Q66, human S65/M66 | mature-chain positions 65–66 | **mouse IQ, human SM** |
| 9S60 is the mouse complex | deposited target chain vs UniProt mature chain | **byte-identical** |
| Tier A pocket positions (ref [2] L56) | residues 155/159/162/194/197/198/201 of 5AEI chain A | **read `QWSQQEW`**, the Arg-binder pocket |
| same, other scaffold | residues 364/368/371/403/406/407/410 of 6SA8 chain A | **read `QWSQQEW`** |

That last pair matters more than it looks. The literature specifies pocket
positions in author numbering, which is not the index into the SEQRES sequence.
Had they been applied as sequence offsets, Tier A would have been built on seven
wrong residues and produced confident, fully-provenanced, meaningless numbers.

**The Tier A `status: TODO` is closed.** All 63 per-variant K_D values across
five pockets were extracted from Table S1 of ref [2] into
`gonogo/data/tier_a_kd.yaml`, line-cited to `PMC12202725#L547-L629`, with
`null` (not zero) where the paper publishes no error margin. The published
physics baselines are in `gonogo/data/baselines.yaml`, per pocket and per input
structure, with the statistic each one is reported in. G3 no longer lacks ground
truth; it lacks model runs.

---

## 3. Gate-by-gate feasibility

| Gate | Can the loop execute it? | The real obstacle |
|---|---|---|
| **G5** decoy rejection | Yes, once any scorer exists | Not compute — *decoy construction*. See below. |
| **G2** specificity | Yes; 6 complexes, ~23 GPU-min at 3 seeds | GPU access only |
| **G1** interface recovery | **No, as written** | DockQ is absent from the proto-tools catalogue. The pass condition is written in DockQ ≥ 0.49 and this stack cannot produce that number. |
| **G3** affinity ranking | Yes, and now with ground truth | Resolution, not availability. See below. |
| **G4** calibration | Yes, with one caveat | Needs an a-priori decision threshold, and can be legitimately undefined. |

### G5: the decoy recipe in the task has a trap

The task says to build non-cognate pairs "drawn from cases already in the set".
Done naively, that produces DARPin 81 × mouse cathepsin B — and `cases.yaml`
itself states that DARPins 81 and 8h6 are **cross-reactive**. That decoy is
plausibly a real binder, and a pipeline penalised for scoring it highly is being
penalised for being right. The harness therefore draws non-cognate partners only
from unrelated systems (anti-MBP off7, anti-lysozyme affitin H4, anti-RBD
repebody A6 against cathepsin B, and 4m3 against MBP and lysozyme) and states
the exclusion in code. Nine decoys, deterministic under a fixed seed.

Second, quieter issue: G5 compares decoys to "the weakest true positive", and the
Tier B positives are B1 (mouse) plus B2/B3 (whose deposited targets, 5MBL and
5MBM, are **human** cathepsin B). The weakest-positive floor is therefore set
across two orthologs. That is defensible but should be deliberate, not
accidental.

### G1: the metric does not exist in this stack

`search_tools("dockq")` returns nothing; the catalogue's structure-alignment
tools are TMalign, USalign, FoldMason, Foldseek and PyMOL RMSD. Three honest
options, in order of preference:

1. `uv pip install DockQ` and call it locally on predicted-vs-deposited pairs —
   CPU-only, closes the gate as written.
2. Rewrite G1 around USalign multimer TM-score, and state the substitution.
3. Report G1 `not_run` with the reason, which is what the harness does today.

What is not acceptable is silently reporting a different metric under the DockQ
threshold, so the harness refuses to.

### G3: ground truth is closed, but the gate may be underpowered

The ladder varies **one residue of a 10-mer peptide at one pocket position**. It
does not vary repeat count or peptide length, contrary to the rationale
`cases.yaml` gives for choosing the system. Asking an ipTM-style confidence to
resolve a 2-fold K_D difference caused by a single side chain is asking for
resolution the metric is not known to have — ref [2] says so itself, at L15,
which is why its authors benchmarked physics methods rather than a folding
model. Expect G3 to fail on ipTM even if the pipeline is sound; that is
information about the *score*, not about the pipeline.

Also note G3's two clauses can disagree. On the Ile pocket the published methods
score rho = 0.143 / 0.253 / 0.291 — a model can beat every published baseline
there and still miss the absolute rho ≥ 0.5 floor. The harness reports the two
clauses separately.

### G4: the threshold must be fixed before the answers are seen

Correctness needs a decision threshold, and choosing one after seeing which
predictions were right is circular. The harness defines it as *the highest decoy
score from G5* — set by the cheapest gate, before anything else is looked at.

G4 can also be legitimately undefined: if every labelled prediction lands on the
correct side, AUC has one class and no value. The harness reports `not_run` with
that reason rather than a perfect score. The self-test run below actually
**fails** G4 (AUC 0.550), which is the useful demonstration that the gate is not
decorative.

---

## 4. The finding that should change the design loop

`design-loop-sketch.md` puts validation before design and gives the Referee a
veto. Correct. But the bench's verifier adds a constraint the sketch does not
have:

> A `not_run` G2 is scored as a **failure**, not as an abstention. The verifier
> caps the whole run's reward at 0.

So "report honestly that we could not run it" and "run it and fail" are worth
the same. The loop must therefore treat *scorer availability* as a **pre-run
condition**, checked and resolved before the run is considered started —
`gonogo preflight` exists for exactly this and exits non-zero when the scorer is
unavailable. An agent team that discovers mid-run that it has no GPU has already
lost the run; an agent team that discovers it in preflight has lost nothing.

This is the concrete version of the Economist role: not "justify expensive
calls" but "prove the expensive call is *possible* before the run counts".

Three smaller things the run surfaced, each now enforced in code rather than in
prose:

- **Construct boundaries are part of the gate.** B4 must use the same construct
  boundaries as B1 — otherwise the comparison measures the construct, not the
  ortholog. Verified identical, in `inputs`.
- **Purification tags are an edit, not a detail.** Six of the fourteen deposited
  binder chains carry `MRGSHHHHHHGS` or `LEVLFQGPLEHHHHHH`. Stripping them
  changes what gets folded, so each removal is recorded on the chain. An early
  version of the stripper truncated DARPin 81 to four residues; the run would
  have completed and produced numbers.
- **The uncertainty band needs a definition, not a value.** G2's verdict is a
  margin against a band. The harness defines it as the spread across independent
  seeds, and says so in the submission: it measures run-to-run instability, not
  model error. A single-seed scorer cannot produce one, so the harness reports
  G2 `not_run` rather than inventing a band.

### The band definition is the harness's own weakest point

Now that the gate has actually issued a `fail`, the definition above deserves
the same scrutiny the harness applies to everything else.

`gates.g2_specificity` takes the band to be the **sample range**, `max − min`,
of the wider of the two cases. The range has a property that matters here: its
expectation grows monotonically with the number of seeds. Every additional
replicate can only hold the range constant or widen it. So under the harness as
written, **G2 gets strictly harder to pass the more evidence you collect** —
which is backwards for a test of whether a margin is real.

The obvious response is to swap in a standard error. That response has to be
resisted, or at least quarantined, because the swap would be made *after* seeing
the verdict, and one of the three candidate definitions flips the result:

| Band definition | Value at n=3 | G2 |
|---|---|---|
| sample range — pre-registered, what the harness uses | 52.67 | **fail** |
| 95% CI half-width on the difference of means (Welch, t≈4.30 at df≈2.1) | 66.2 | **fail** |
| bare SE of the difference | 15.4 | pass |

Two of the three say fail, so **the reported verdict is not an artifact of the
range choice** — which is the reason it is safe to leave the code alone. The
third is the least defensible of the three (an SE is a standard error of a mean,
not an interval you compare a margin against) and it is the only one that would
have rescued the run. Changing the definition now would be exactly the post-hoc
tuning that this task guards against elsewhere, where it fixes the G4 decision
threshold at the highest decoy score precisely "so that it cannot be tuned after
the fact".

The honest split, and what this document recommends:

1. **Leave the shipped verdict as it stands.** The range band is what was
   declared, and it says fail. Two of three defensible bands agree.
2. **Fix the definition for the next run, before that run is scored** — a CI
   half-width on the difference of means, which *shrinks* as `1/√n` and is
   therefore answerable by more seeds. Declare it in `cases.yaml` first.
3. **Re-run B1/B4 at n≥10 only under the new, pre-declared band.** Under the
   current one the extra seeds cannot help, and buying GPU for a foregone
   conclusion is the Economist's job to refuse.

The general lesson for `design-loop-sketch.md`: a gate that compares a margin to
an uncertainty band is only as good as the estimator behind the band, and that
estimator has to be pinned down at the same time as the threshold — before any
numbers are seen. The harness pinned the threshold and left the estimator
implicit. That is the gap this control run found in itself.

---

## 5. Cost, and the phase plan

From `gonogo preflight`, at three seeds per complex, ~60 GPU-s per complex per
seed on an A100-class GPU at ~$4/GPU-hour:

| Phase | What | Complexes | Cost | Gate |
|---|---|---|---|---|
| **P0** | resolve inputs, verify claims, close the K_D TODO, build decoys | 23 | **$0** | — |
| **P1** | G5 then G2 | 15 | **~$4** | if G2 fails, stop and report |
| **P2** | Tier C sweep, then G1/G4 | +8 | ~$2 | only bought if G2 passes |
| **P3** | Tier A ladder, one scaffold | +63 | ~$13 | only worth buying with an affinity-capable score |
| **P4** | Tier A on the second scaffold (the robustness check ref [2] asks for) | +63 | ~$13 | optional |

P0 is where nearly all the errors were caught, and it cost nothing. That is the
argument for the phase split, in one line.

---

## 6. Demos, and what each one proves

**Demo 1 — the real state of the world.**

```
.venv/bin/python -m gonogo --scorer boltz2 run && .venv/bin/python -m gonogo verify
```

Resolves 23 complexes, buys G5 (12 complexes, ~$2.73), passes it, buys G2
(3 complexes, ~$0.93), fails it, and then **stops** — G1/G3/G4 are never
purchased. Writes a schema-valid submission and scores 0.0 on the verifier,
4 of 5 checks passing.

This is the demonstration that the Referee's veto is real, and it is a better
demo than the one this document originally planned for. The veto is not firing
on a missing tool; it is firing on a model that produced numbers, ranked them in
the right order, and still could not be trusted — B4's seeds span 14.6 to 67.3.
Nothing downstream was bought on top of that. No number was invented to make the
file look finished.

**Demo 2 — the harness self-test.**

```
.venv/bin/python -m gonogo --scorer synthetic --out demo/selftest run --allow-synthetic
.venv/bin/python -m gonogo --out demo/selftest verify        # 5 passed
```

Exercises the gate arithmetic, the submission schema, the decoy naming, the
band/margin logic and the verifier wiring end to end — **with numbers that are a
hash of the input sequences and mean nothing**. Every prediction is stamped
`synthetic: true`, the scorer refuses to run without `--allow-synthetic`, and the
CLI refuses to write its output into the bench submission path. It proves the
plumbing, and it is labelled loudly enough that it cannot be mistaken for a
result.

---

## 7. What to fix in `cases.yaml`

Feeding back to whoever owns the case set:

1. Tier A's rationale says repeat count and peptide length "can be varied
   independently". The ground truth varies neither — it is a single-residue
   peptide substitution series. The rationale should match the data.
2. `task.md` says G2 is "four complexes total"; Tier B defines six cases
   (B1–B6). Either the extra two negatives are in scope or they are not.
3. The G5 recipe should exclude cross-reactive binders from non-cognate pairing,
   or it manufactures false decoys.
4. B6 (DARPin E3_5) cites no structure. Its sequence is available from its own
   deposit, 1MJ0, which the harness now uses.
5. G1's pass condition is written in a metric this toolchain cannot compute.

None of these is fatal. All five are the kind of thing a control run is supposed
to surface before a design run inherits them.

--------
REFERENCES
[1] Zarić, M. et al. "Structural and Proteomic Analysis of the Mouse Cathepsin
    B-DARPin 4m3 Complex Reveals Species-Specific Binding Determinants."
    *International Journal of Molecular Sciences* (2025). doi:10.3390/ijms262411910
    https://paperclip.gxl.ai/citations/papers/PMC12732864#L10,L19,L21,L22
[2] Ayyildiz, M., Noske, J., Gisdon, F. J., Kynast, J. P. & Höcker, B.
    "Evaluation of Physics-Based Protein Design Methods for Predicting Single
    Residue Effects on Peptide Binding Specificities." *Journal of Computational
    Chemistry* (2025). doi:10.1002/jcc.70160
    https://paperclip.gxl.ai/citations/papers/PMC12202725#L15,L41,L42,L55,L56,L57,L547-L629
