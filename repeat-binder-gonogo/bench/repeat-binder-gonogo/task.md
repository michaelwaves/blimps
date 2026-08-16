---
schema_version: "1.3"
metadata:
  author_name: "reAgent"
  difficulty: hard
  category: capability
  tags: [protein-design, binder, validation, proto, repeat-proteins]
agent:
  timeout_sec: 5400
verifier:
  timeout_sec: 300
environment:
  cpus: 2
  memory_mb: 8192
---
# repeat-binder-gonogo

## prompt

Validate a binder-scoring pipeline against repeat-protein complexes whose
outcomes are already known, and report whether it is fit to be trusted on novel
designs. This is a control run: you are not designing anything.

`cases.yaml` holds the validation set (positives, hard negatives, and gates) and
`REFERENCES.md` holds the line-cited provenance for every measurement. Read both
before starting.

### What to do

1. Run gate **G5** (decoy rejection) first — it is CPU-cheap. Build the decoys
   yourself: sequence-shuffled binders, and non-cognate binder x target pairings
   drawn from cases already in the set.
2. Run gate **G2** (specificity discrimination) next: cases B1 and B4 are the
   same DARPin against two 83%-identical cathepsin B orthologs, one a binder and
   one not. Four complexes total. **If G2 fails, stop and report** — do not spend
   GPU budget on the remaining tiers.
3. Only if G2 passes, run **G1** (interface recovery) over Tier A and Tier C,
   then **G3** (affinity ranking) over the Tier A ladder, then **G4**
   (calibration) over everything scored.

Tier A's per-variant K_D table is marked `status: TODO` in `cases.yaml`. Either
extract it from the source cited there, or report G3 as not-run. Do not
substitute estimated affinities.

### Rules

- Every score must carry a confidence value from the model that produced it
  (ipTM, PAE, ensemble spread, likelihood). A bare point estimate is a failed
  submission, not a low-scoring one.
- Every prediction must record which tool produced it and where it ran.
- Log failures and dead ends. A gate that could not be run is reported as
  `not_run` with a reason — never silently omitted, never guessed.

### Output

Write `results/predictions.json`:

```json
{
  "predictions": [
    {"case_id": "B1", "score": 12.4, "confidence": 0.88,
     "confidence_metric": "iptm", "tool": "boltz2-prediction", "ran_on": "modal"}
  ],
  "gates": [
    {"gate_id": "G2", "status": "pass", "uncertainty_band": 1.2,
     "evidence": "score(B1)=12.4 > score(B4)=6.1, margin 6.3 > band 1.2"}
  ],
  "notes": "anything that did not work, and why"
}
```

**`score` is higher-is-better** — a stronger predicted interaction gets a larger
number. If your model emits an energy (lower is better), submit its negation and
say so in `notes`. One convention across all cases; the verifier compares them
directly.

`status` is one of `pass`, `fail`, `not_run`. A gate reported `pass` or `fail`
must carry `evidence`; `G2` must additionally carry `uncertainty_band`, the
margin its verdict is judged against. Decoy predictions use case ids prefixed
`decoy_`.
