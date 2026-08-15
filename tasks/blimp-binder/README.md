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
RF3 alongside the GvpA copies (GvpC omitted — no sequence), and requires a
binder-GvpA interface PAE <= 10 A against the best-packed GvpA copy, plus a
30-150 residue length window. Scores land in `/logs/verifier/scores.json`.

Why PAE and not ipTM: RF3's `summary_confidences` exposes a single global
`iptm` and no per-chain-pair version (keys are `chain_pair_pae`,
`chain_pair_pae_min`, `chain_pair_pde`, `chain_pair_pde_min`, `chain_ptm`,
`iptm`, `ptm`, `overall_pae/pde/plddt`, `has_clash`, `ranking_score`). Global
ipTM here is dominated by GvpA-GvpA contacts: folded from sequence alone RF3
does not reassemble the shell lattice, scoring ~0.20 regardless of the binder.
The chain-pair PAE is local to the binder and unaffected by that.

Deliberately minimal: one seed, one diffusion sample, sequence-only fold — the
score is a proxy for binding the real vesicle surface, not a claim about it.
`MAX_INTERFACE_PAE_ANGSTROM` is a guess until a real RFD3 binder has been
scored; `scores.json` also records per-chain PAE, global ipTM, pLDDT and clash
for recalibration.

`solution/` is a smoke-test stub (a generic helix, expected to fail the ipTM
gate), not an oracle.

## Running

```bash
T=tasks/blimp-binder
uv run harbor task start-env -p $T -e modal
uv run harbor run -p $T -a claude-code -m claude-opus-5 -e modal
```

Smoke test of the full pipeline: `uv run harbor run -p $T -a oracle -e modal`
(expect reward 0.0 — the junk-helix control fails the PAE gate).
