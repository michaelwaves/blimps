# Lattice-aware GvpA binder pipeline

This directory is the completed continuation of the handoff. It adds true de
novo backbone generation and preserves the earlier control run's limits on what
the resulting scores are allowed to mean.

## Architecture

```text
7R1C-derived five-chain GvpA patch
        |
        +--> constrained arm: 28 cross-species GvpC-derived hotspots
        |         RFdiffusion3 -> 8 binder backbones
        |
        +--> free arm: no hotspots, center-of-mass origin
                  RFdiffusion3 -> 8 binder backbones
                                  |
                         geometry/referee gate
                         - footprint vs hotspots
                         - backbone clashes < 2 A
                                  |
                       ProteinMPNN (chain A only)
                       4 sequences/backbone -> 64 candidates
                                  |
                         inexpensive selection
                         - reject clashing backbones
                         - one sequence/backbone
                         - rank by MPNN perplexity
                                  |
                         Boltz-2 co-fold triage
                         2 backbones/arm at seed 0
                         winner/arm confirmed at seeds 1-2
                                  |
                         result + uncertainty + claim limit
```

The control run sits *above* this flow as an evidence gate:

- G5 pass: scores contain enough signal for decoy rejection.
- G2 pass: specificity comparisons are licensed.
- G1 fail: Boltz-2 does not recover repeat-protein interfaces reliably, so its
  co-fold scores are triage evidence, not proof of the epitope or pose.
- G3 fail: the score is not an affinity ranking.
- G4 pass: confidence can help decide where to spend follow-up effort.

## Experiment and result

Both arms used identical budgets: eight RFdiffusion3 backbones, four
ProteinMPNN sequences per backbone, and the same downstream selection policy.

| result | constrained | free |
|---|---:|---:|
| RFdiffusion3 backbones | 8 | 8 |
| ProteinMPNN candidates | 32 | 32 |
| hotspot residues contacted per backbone | **6–10** | **0** |
| closest hotspot distance | **0.7–2.5 A** | **14.5–22.7 A** |
| backbones with a <2 A backbone clash | **4/8** | **0/8** |
| winner mean binder-target ipTM, 3 seeds | **0.2605** | **0.0951** |
| winner ipTM seed range | **0.1998** | **0.0539** |

The hotspot constraint deterministically changes where RFdiffusion3 places the
binder and increases the downstream Boltz-2 triage score. It also creates a
50% backbone-clash failure rate, and the surviving winner remains low-confidence
and highly seed-sensitive. Therefore the pipeline has produced two novel
sequences and a useful localization experiment, **not a validated binder**.

## Outputs

- `pipeline_results.json` — complete experiment policy, screen, winners, and
  all Boltz-2 replicates.
- `candidates.json` / `candidates.csv` — all 64 sequences joined to their
  backbone geometry and ProteinMPNN metrics.
- `best_binders.fasta` — constrained and free winners.
- `rfd3_inputs/*.cif` — the 16 generated complexes with format-correct names.
- `boltz/` — cached raw outputs and parsed scores; reruns reuse completed calls.

## Reproduce

```bash
.venv/bin/python gvpa_rfd3.py --arm constrained --n-designs 8
.venv/bin/python gvpa_rfd3.py --arm free --n-designs 8
.venv/bin/python gvpa_full_pipeline.py mpnn --n-per-backbone 4
.venv/bin/python gvpa_full_pipeline.py compile
.venv/bin/python gvpa_full_pipeline.py boltz --top-per-arm 2
```

The Boltz stage is intentionally cost-gated: it screens two distinct,
clash-free backbones per arm and only replicates the winner from each arm.
