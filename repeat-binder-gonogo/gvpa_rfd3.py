"""De novo binder backbones against the GvpA surface, constrained and not.

    .venv/bin/python gvpa_rfd3.py --arm constrained   # hotspots = the GvpC site
    .venv/bin/python gvpa_rfd3.py --arm free          # no hotspots

Two arms, because the comparison is the experiment. G1 showed the scorer cannot
locate an epitope on a repeat protein (0/6, with the target epitope missed
outright on three). If that carries over here, the free arm should scatter
binders over the surface while the constrained arm sits on the supplied site.
Running only the constrained arm would hide that.

WHAT THE HOTSPOTS ARE, AND ARE NOT

They come from superimposing Anabaena GvpC (8GBS) onto this Megaterium GvpA
surface (7R1C-derived): 72% identity, 1.76 A CA RMSD. So they mark where the
natural binder sits on a *related* shell, not a measured Megaterium epitope.
Designs that satisfy them are designs against a hypothesis, and the writeup has
to say so.

Nothing here is scored against the GvpC placement with DockQ. That reference is
a backbone-only poly-UNK trace with three GvpA residues inside 5 A; it can say
roughly where a binder landed and nothing about how well it fits.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent
WORK = REPO / "designs" / "gvpa"

BINDER_LENGTH = "60-90"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=["constrained", "free"], required=True)
    ap.add_argument("--n-designs", type=int, default=8)
    args = ap.parse_args()

    from proto_tools.mcp.tools import run_tool

    site = json.loads((WORK / "site.json").read_text())
    hotspots = site["gvpc_footprint_residues"]
    target = WORK / "gvpa_target.pdb"

    spec = {
        "input_structure": str(target),
        # binder chain, then the five fixed target chains
        "contig": f"{BINDER_LENGTH},/0,A2-66,/0,B2-66,/0,C2-66,/0,D2-66,/0,E2-66",
        # Origin inference is a per-design parameter in the deployed schema,
        # not a run-level RFdiffusion3Config field.
        "infer_ori_strategy": "hotspots" if args.arm == "constrained" else "com",
    }
    if args.arm == "constrained":
        spec["select_hotspots"] = ",".join(hotspots)

    out_dir = WORK / f"rfd3_{args.arm}"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"arm={args.arm}  target=325 residues  designs={args.n_designs}")
    print(f"hotspots: {spec.get('select_hotspots', '(none — free arm)')}\n", flush=True)

    result = run_tool(
        "rfdiffusion3-design",
        inputs={"design_specs": [spec]},
        config={"n_batches": args.n_designs, "diffusion_batch_size": 1,
                "seed": 0},
        device="modal",
        output_dir=str(out_dir),
    )
    if not result.get("ok", True):
        print(f"FAILED: {result.get('error')}")
        return 1

    (WORK / f"rfd3_{args.arm}.json").write_text(
        json.dumps({"arm": args.arm, "hotspots": hotspots if args.arm == "constrained" else None,
                    "hotspot_provenance": "Anabaena GvpC (8GBS) superimposed onto Megaterium GvpA "
                                          "(7R1C-derived); 72% identity, 1.76 A CA RMSD. A hypothesis "
                                          "about the site, not a measured Megaterium epitope.",
                    "result": result.get("result")}, indent=2, default=str) + "\n")
    print(f"wrote {WORK / f'rfd3_{args.arm}.json'}  structures in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
