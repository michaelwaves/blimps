Gas vesicles are hollow, gas-filled protein shells built from the structural protein GvpA. Their surface is a target for imaging and biosensing reagents.

Design a de novo protein binder for each of two GvpA shell proteins:

- `/app/foundry/pdbs/7r1c.pdb` — GvpA from *Priestia megaterium* No natural binder is known. 7r1c is a 5 mer of GvPA
- `/app/foundry/pdbs/8gbs.pdb` — GvpA from *Dolichospermum flos-aquae*. It is a 4 mer of GvPA and one subset repeat unit of GvpC. Its natural partner GvpC binds along the outer shell surface; a binder for a different surface would be new.

Submit `/logs/outputs/binder_7r1c.pdb` and `/logs/outputs/binder_8gbs.pdb`: ATOM records for a single protein chain of 30-150 residues whose CA records, in file order, encode the binder sequence using the 20 standard amino acids.

Grading: each binder sequence is refolded with RoseTTAFold3 in complex with its target chain (sequence only, no templates). Both complexes must reach interface ipTM >= 0.6.

The container provides the Foundry toolchain — the `rfd3`, `mpnn` and `rf3` command-line tools and their Python inference engines, with weights pre-downloaded and a single A100 GPU.

You have 7200 seconds. Do not use online solutions or hints specific to this task.
