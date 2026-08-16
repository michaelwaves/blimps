# re:AGENT submission — repeat-binder-gonogo + GvpA design run

Two connected pieces of work. A go/no-go validation that asks whether a
binder-scoring pipeline can be trusted, and a design run against a real target
that inherits the answer.

```
trajectory.jsonl    174 records, every phase, every dead end
manifest.json       record counts, validation checks, artifact hashes
artifacts/          every file a claim in the trajectory points at
```

## The trajectory is segmented, not pruned

174 records across four phases. The control run alone is **seventeen separate
invocations** — `t` restarts on each, and thirty lines are byte-identical
repeats of per-attempt setup.

Those repeats were kept. They are evidence that setup re-ran after each abort,
and the rubric checks honest reporting as a pass/fail gate before it scores
anything. Instead of merging or deleting, every record carries a `run_id`, a
`phase`, and a globally monotonic `seq`, so a reader can follow one ordered
narrative or filter to a single attempt without anything having been removed.

| outcome | count |
|---|---|
| ok | 144 |
| veto | 13 |
| gap | 9 |
| blocked | 8 |

The 30 non-`ok` records are the point, not the noise: four runs aborted on a
Modal billing block, one deliberate veto override, a Tier A sweep killed at
30/63 by an upstream connection reset, three declined deploys, and two of my own
analysis errors caught and corrected.

The `gvpa-design` records are marked `logged: retrospective` — that work ran
through standalone scripts outside the scribe, and the records were
reconstructed from the artifacts each step wrote. They are not live
instrumentation and are labelled so they cannot be mistaken for it.

## Part 1 — the control run

Boltz-2 against repeat-protein complexes with published outcomes.

| Gate | Verdict | What passing would license |
|---|---|---|
| **G5** decoy rejection | **pass** | that scores carry signal at all |
| **G2** specificity | **pass** | specificity claims — hits target, not paralog |
| **G1** interface recovery | **fail** | epitope-directed design |
| **G3** affinity ranking | **fail** | ranking and affinity maturation |
| **G4** calibration | **pass** | triage — where to spend wet-lab budget |

Bench verifier 5/5. Three findings worth the space:

**G2 was failed by its estimator, not the model.** The band was the sample range
(`max − min`), which is set by two draws and discards the rest. B4's extremes
both landed in the first three seeds, so twelve seeds moved it by nothing — it
reads 52.67 at n=3 and at n=12. A CI half-width falls 66 → 12 and crosses the
margin at n=4. At twelve seeds the separation is Welch t = 6.08. The estimator
was then re-declared in `cases.yaml` **before** being applied, with a 5-seed
floor, and the superseded n=3 submission preserved unmodified.

**The scorer fails on exactly the scaffold class this set is about.** 0/6 repeat
proteins recover their interface; 2/2 non-repeat do (Fisher p = 0.036). And
`cases.yaml` obscured it — affitins are listed as repeat proteins but are
Sac7d-derived single domains. Both successes are the two cases that are not
repeat proteins.

**More compute made confidence worse, not better.** Diagnosis showed the binder's
own paratope is placed roughly right (Jaccard 0.36–0.54) while the target
epitope is missed outright (0.00 on three of six). Retrying at 5 seeds × 10
recycles, selecting by the model's own ipTM with the deposit sealed out, moved
C_1SVX from ipTM 0.664 → **0.927** with DockQ 0.0105 → **0.0106**. Knowledge
limit, not sampling limit.

## Part 2 — the GvpA design run

Target: a generated 5-rib GvpA surface patch. Reference: the native *Anabaena*
gas vesicle shell, 8GBS.

**It is not de novo design, and is not reported as such.** RFdiffusion3 was
declined at the approval prompt three times, so no backbone was generated. What
ran is inverse folding on the GvpC backbone placement — fold and pose given,
sequence designed.

That returned **+0.4 percentage points above a composition-shuffled baseline** —
no positional signal. The cause was visible in advance: GvpC in 8GBS is a
backbone-only poly-UNK trace with three GvpA residues inside 5 Å. There was no
interface to design against, and DockQ against it would have been meaningless.

**The result worth showing is the curvature check.** A patch cut from an 85 nm
cylinder has to bow like one, and curvature is a property the generator never
optimised against:

| | generated patch | Ana GV, measured |
|---|---|---|
| implied diameter | **35.6 nm** | **85 ± 4 nm** |
| rib sagitta over 121 Å | 10.35 Å | 4.33 Å expected |

**2.4× over-curved for its own species.** All five ribs agree exactly, which is
what makes the number trustworthy — a cylinder curves the same way in every rib.
Reference values are line-cited to Dutka et al. 2023, PMC10185304 L22 and L83 —
the source study for 8GBS itself.

My first attempt at this returned 12.6 nm and was wrong. It swept cylinder axes
and fitted circles, reporting a "well constrained" 137° arc; rendering the patch
in PyMOL showed a flat sheet, and a 137° arc would look like a letter C. The fit
had locked onto internal scatter. Both the error and the correction are in the
trajectory.

## Reproducing

```bash
python -m gonogo preflight    # free: what can run, what it costs, what blocks it
python -m gonogo inputs       # free: resolve 23 complexes, re-derive the case set's claims
python -m gonogo verify       # free: the bench's own verifier
python gvpa_curvature.py      # free: the curvature measurement
```

`inputs` re-derives every factual claim in `cases.yaml` from primary sources and
prints pass/fail per claim. That phase caught nearly every error in the case
set, at zero cost.

## What this pipeline currently licenses

Specificity claims, and triage of binding calls. **Not** epitope-directed design
and **not** affinity ranking — G1 and G3 both fail, and those are the two a
repeat-protein binder campaign actually needs. The GvpA run then found the
target structure itself to be 2.4× over-curved, which is a problem upstream of
any binder.
