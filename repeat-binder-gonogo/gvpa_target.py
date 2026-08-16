"""Put the GvpC site onto the GvpA patch, and cut a designable target out of it.

    .venv/bin/python gvpa_target.py

Two inputs that do not share a coordinate frame:

  GvpA_surface_patch_5ribs_x11.pdb   55 GvpA subunits, 3575 residues, 266 A across
  8GBS.cif                           4 GvpA + 1 GvpC, the native Ana shell model

The patch is the design target; 8GBS says where the natural binder sits. To use
one against the other, 8GBS has to be superimposed onto the patch: align a GvpA
subunit, then carry GvpC along with the same transform.

Then the patch is cut down. A binder-design pipeline wants a few hundred
residues, not 3575, so this keeps the subunits under the GvpC footprint plus one
shell of neighbours and writes them out as the target.

A caveat that belongs on the result, not in a footnote: GvpC in 8GBS is a
backbone-only poly-UNK trace (N, CA, C, O, CB; no residue identities) and only
three GvpA residues fall within 5 A of it. It marks roughly where the natural
binder lies. It is not a resolved interface, and nothing downstream should treat
it as one.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import biotite.structure.io.pdbx as pdbx
from biotite.structure import get_residues, superimpose
from biotite.structure.io.pdb import PDBFile

REPO = Path(__file__).resolve().parent
OUT = REPO / "designs" / "gvpa"
CONTACT_A = 8.0     # GvpC is a coarse placement; 5 A finds almost nothing
NEIGHBOUR_A = 12.0  # one shell of GvpA subunits around the footprint


def ca(arr):
    return arr[arr.atom_name == "CA"]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    shell = pdbx.get_structure(pdbx.CIFFile.read(REPO / "8GBS.cif"), model=1, use_author_fields=True)
    shell = shell[~shell.hetero]
    patch = PDBFile.read(str(REPO / "GvpA_surface_patch_5ribs_x11.pdb")).get_structure(model=1)

    gvpc = shell[shell.chain_id == "C"]
    ref = shell[shell.chain_id == "A1"]
    ref_ca = ca(ref)

    # Which patch subunit does a GvpA monomer land on best?
    best = None
    for chain in sorted(set(patch.chain_id)):
        cand = ca(patch[patch.chain_id == chain])
        n = min(len(ref_ca), len(cand))
        if n < 40:
            continue
        fitted, transform = superimpose(cand.coord[:n], ref_ca.coord[:n])
        rmsd = float(np.sqrt(((fitted - cand.coord[:n]) ** 2).sum(axis=1).mean()))
        if best is None or rmsd < best["rmsd"]:
            best = {"chain": chain, "rmsd": rmsd, "transform": transform, "n": n}

    print(f"best GvpA superposition: patch chain {best['chain']}, "
          f"CA RMSD {best['rmsd']:.2f} A over {best['n']} atoms")
    if best["rmsd"] > 3.0:
        print("  WARNING: poor fit — the patch may not be the same GvpA conformer as 8GBS")

    # Carry GvpC into the patch frame with the transform that aligned its GvpA.
    gvpc_moved = gvpc.copy()
    gvpc_moved.coord = best["transform"].apply(gvpc.coord)
    print(f"GvpC placed in patch frame: centroid "
          f"{np.round(gvpc_moved.coord.mean(axis=0), 1)}")

    # Footprint: patch residues near the placed GvpC.
    d = np.linalg.norm(patch.coord[:, None, :] - gvpc_moved.coord[None, :, :], axis=-1)
    near = d.min(axis=1)
    hot = {(c, int(r)) for c, r in zip(patch.chain_id[near < CONTACT_A], patch.res_id[near < CONTACT_A])}
    hot_chains = sorted({c for c, _ in hot})
    print(f"GvpC footprint on the patch: {len(hot)} residues across chains {hot_chains}")

    # Target = the footprint subunits themselves. Adding a neighbour shell pushed this
    # to 1690 residues, which is past what a binder-design pipeline will take; the five
    # subunits GvpC actually lies across are the surface a binder has to read.
    keep = set(hot_chains)
    target = patch[np.isin(patch.chain_id, sorted(keep))].copy()

    # Contig strings and hotspot selectors are parsed as <CHAIN><RESID>; the patch uses
    # lowercase ids (l..p), which those parsers generally do not accept. Relabel to A..E
    # and carry the mapping so the footprint can be quoted in the same alphabet.
    relabel = {old: chr(ord("A") + i) for i, old in enumerate(sorted(keep))}
    target.chain_id = np.array([relabel[c] for c in target.chain_id])
    hot = {(relabel[c], r) for c, r in hot if c in relabel}
    hot_chains = sorted(relabel[c] for c in hot_chains)
    print(f"relabelled chains: {relabel}")
    n_res = len(set(zip(target.chain_id, target.res_id)))
    print(f"trimmed target: {len(keep)} subunits, {n_res} residues "
          f"(from 55 subunits, 3575 residues)")

    tgt_file = PDBFile()
    tgt_file.set_structure(target)
    tgt_file.write(str(OUT / "gvpa_target.pdb"))

    ref_file = PDBFile()
    ref_file.set_structure(gvpc_moved)
    ref_file.write(str(OUT / "gvpc_reference.pdb"))

    hotspots = sorted(f"{c}{r}" for c, r in hot)
    (OUT / "site.json").write_text(json.dumps({
        "superposition": {"patch_chain": best["chain"], "ca_rmsd_A": round(best["rmsd"], 3),
                          "n_atoms": best["n"]},
        "contact_cutoff_A": CONTACT_A,
        "gvpc_footprint_residues": hotspots,
        "footprint_chains": hot_chains,
        "target_chains": sorted(keep),
        "target_n_residues": n_res,
        "caveat": "GvpC in 8GBS is a backbone-only poly-UNK trace with 3 GvpA residues inside 5 A. "
                  "It localises the natural binder; it is not a resolved interface and must not be "
                  "scored with DockQ.",
    }, indent=2) + "\n")

    print(f"\nwrote {OUT}/gvpa_target.pdb, gvpc_reference.pdb, site.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
