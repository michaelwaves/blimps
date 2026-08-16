"""Design a binder sequence for the GvpC site on the GvpA shell.

    .venv/bin/python gvpa_design.py [--n 32]

WHAT THIS IS, AND WHAT IT IS NOT

RFdiffusion3 was not deployed, so no new backbone is generated here. This is
inverse folding, not de novo design: the GvpC backbone placement from 8GBS is
reused as the scaffold, and ProteinMPNN designs a sequence onto it against the
GvpA surface. Calling that "de novo binder design" would overstate it — the
fold and the pose are given; only the sequence is designed.

What it does buy is a real, checkable question. GvpC's native sequence over the
same segment is known (P09413, residues 85-117), so the designs can be scored
for native recovery, and the co-folded complex can be checked against where GvpC
actually sits. Both are honest metrics with a ground truth behind them.

The standing caveat travels with every number: the GvpC placement in 8GBS is a
backbone-only poly-UNK trace with three GvpA residues inside 5 A. It says
roughly where the natural binder lies, and nothing about its side chains.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from biotite.structure import concatenate
from biotite.structure.io.pdb import PDBFile

REPO = Path(__file__).resolve().parent
WORK = REPO / "designs" / "gvpa"

# P09413 GVPC_DOLFA residues 85-117 — the segment 8GBS numbers its UNK trace as.
NATIVE_GVPC = "HKELQETSQQFLSATAQARIAQAEKQAQELLAF"
TARGET_CHAINS = ["A", "B", "C", "D", "E"]
BINDER_CHAIN = "F"


def build_complex() -> Path:
    """GvpA target (fixed) plus the placed GvpC backbone as the chain to design."""
    target = PDBFile.read(str(WORK / "gvpa_target.pdb")).get_structure(model=1)
    binder = PDBFile.read(str(WORK / "gvpc_reference.pdb")).get_structure(model=1).copy()
    binder.chain_id = np.array([BINDER_CHAIN] * len(binder))
    # ProteinMPNN reads residue identity off the backbone; UNK is not an amino acid.
    # Glycine is the neutral placeholder — no side chain to bias the design.
    binder.res_name = np.array(["GLY"] * len(binder))
    binder = binder[np.isin(binder.atom_name, ["N", "CA", "C", "O"])]

    combined = concatenate([target, binder])
    out = WORK / "design_complex.pdb"
    f = PDBFile()
    f.set_structure(combined)
    f.write(str(out))
    print(f"complex: {len(TARGET_CHAINS)} target chains fixed + chain {BINDER_CHAIN} "
          f"({len(set(binder.res_id))} residues) to design")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=32, help="sequences to sample")
    ap.add_argument("--temperature", type=float, default=0.2)
    args = ap.parse_args()

    from proto_tools.mcp.tools import run_tool

    complex_path = build_complex()
    print(f"designing chain {BINDER_CHAIN}, holding {','.join(TARGET_CHAINS)} fixed", flush=True)

    result = run_tool(
        "proteinmpnn-sample",
        # Only chain F is redesigned; every chain omitted here is held fixed, so the
        # GvpA surface keeps its native sequence and acts purely as context.
        inputs={"inputs": [{
            "structure": str(complex_path),
            "chains_to_redesign": {"chains": [BINDER_CHAIN]},
        }]},
        config={"num_sequences_per_structure": args.n, "temperature": args.temperature},
        device="modal",
    )
    if not result.get("ok", True):
        print(f"FAILED: {result.get('error')}")
        return 1

    payload = result["result"]
    designs = payload.get("designs") or payload.get("sequences") or payload
    (WORK / "proteinmpnn_raw.json").write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"\nwrote raw output to {WORK / 'proteinmpnn_raw.json'}")
    print(f"native GvpC 85-117: {NATIVE_GVPC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
