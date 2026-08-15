# blimps/blimp-binder

De novo binder design against gas vesicle shell protein (GvpA), on the Foundry
image (`rfd3` / `mpnn` / `rf3`, one A100).

## Targets

`environment/pdbs/` holds coordinate-only biological assemblies (residues 2-66
per GvpA copy), rebuilt from RCSB by `authoring/build_targets.py`:

- `7r1c.pdb` — GvpA pentamer, *Priestia megaterium*
  ([7R1C](https://www.rcsb.org/structure/7R1C)), chains A-E. Thick,
  pressure-resistant shell; no natural binder known.
- `8gbs.pdb` — GvpA tetramer, *Dolichospermum flos-aquae*
  ([8GBS](https://www.rcsb.org/structure/8GBS)), chains A-D, plus one GvpC
  repeat as chain G. GvpC is deposited as poly-UNK, so it marks the occupied
  surface but carries no sequence.

## Scoring

`tests/test_outputs.py` reads each binder sequence from the CA records of
`/logs/outputs/binder_<target>.pdb`, refolds it from sequence with fixed-seed
RF3 alongside the GvpA copies (GvpC omitted — no sequence), and requires the
best binder-to-GvpA chain-pair ipTM >= 0.60 for both targets, plus a 30-150
residue length window. Scores land in `/logs/verifier/scores.json`.

Deliberately minimal: one seed, one diffusion sample, and a sequence-only fold —
RF3 reassembles the GvpA copies itself rather than being held to the deposited
lattice, so the score is a proxy for binding the real vesicle surface.
Calibrate `MIN_IPTM` after the first real agent runs.

`solution/` is a smoke-test stub (a generic helix, expected to fail the ipTM
gate), not an oracle.

## Running

```bash
T=tasks/blimp-binder
uv run harbor task start-env -p $T -e modal
uv run harbor run -p $T -a claude-code -m claude-opus-5 -e modal
```

Smoke test of the full pipeline: `uv run harbor run -p $T -a oracle -e modal`
(expect reward 0.0; the junk-helix control scored ipTM 0.16 / 0.31).
