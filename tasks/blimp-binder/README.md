# blimps/blimp-binder

De novo binder design against gas vesicle shell protein (GvpA), on the Foundry
image (`rfd3` / `mpnn` / `rf3`, one A100).

## Targets

`environment/pdbs/` holds coordinate-only single chains (residues 2-66) pulled
from RCSB:

- `7r1c.pdb` — GvpA, *Priestia megaterium* ([7R1C](https://www.rcsb.org/structure/7R1C)),
  thick pressure-resistant shell, no known natural binder.
- `8gbs.pdb` — GvpA, *Dolichospermum flos-aquae* ([8GBS](https://www.rcsb.org/structure/8GBS));
  the deposited structure also contains GvpC (modelled as poly-UNK), which is
  dropped here, so the agent has to find its own site.

## Scoring

`tests/test_outputs.py` reads each binder sequence from the CA records of
`/logs/outputs/binder_<target>.pdb`, refolds it from sequence with fixed-seed RF3
in complex with the target chain, and requires interface ipTM >= 0.60 for both
targets (plus a 30-150 residue length window). ipTM values land in
`/logs/verifier/scores.json`.

Deliberately minimal: single seed, one diffusion sample, no shell-lattice
context — ipTM against a GvpA monomer is a proxy, not a claim about binding the
assembled vesicle. Calibrate `MIN_IPTM` after the first real runs.

## Running

```bash
T=tasks/blimp-binder
uv run harbor task start-env -p $T -e modal
uv run harbor run -p $T -a claude-code -m claude-opus-5 -e modal
```
