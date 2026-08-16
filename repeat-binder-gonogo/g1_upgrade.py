"""Is the epitope failure a sampling limit or a knowledge limit?

    .venv/bin/python g1_upgrade.py [CASE ...]

g1_diagnosis.py found that the Tier C failures are not register errors. The
binder's own paratope is roughly right (Jaccard 0.36-0.54) while the target
epitope is missed outright (0.00 on three of six). The model knows which face of
the DARPin binds; it does not know where on the target to put it.

Two explanations, and they call for different responses:

  sampling limit   the right epitope is inside the model's distribution but a
                   single seed with one diffusion sample rarely lands on it.
                   More samples plus confidence-based selection should recover it.
  knowledge limit  the model has no signal for this epitope at all. More
                   sampling buys nothing and the pipeline needs a different
                   input, not a bigger budget.

The test is a straight sampling increase: 5 seeds x 5 diffusion samples x 10
recycling steps, keeping the pose the MODEL ranks best.

WHAT THIS DELIBERATELY DOES NOT DO
Selection is by ipTM — the model's own confidence — never by DockQ. Choosing the
pose that best matches the deposit and then scoring it against the deposit would
report the quality of the search, not the quality of the model, and would make
G1 meaningless. The deposited structure enters exactly once, at scoring time,
after the pose is already fixed.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

from gonogo import caseset

REPO = Path(__file__).resolve().parent
WORK = REPO / "logs" / "g1_upgraded"
BASE = REPO / "logs" / "g1_structures"
OUT = REPO / "bench" / "repeat-binder-gonogo" / "results"

SEEDS = [0, 1, 2, 3, 4]
DIFFUSION_SAMPLES = 5
RECYCLING_STEPS = 10

# The six G1 failures, worst-diagnosed first.
DEFAULT_CASES = ["C_1SVX", "C_8RCI", "C_7YCO"]


def best_pose(spec, dest: Path) -> dict:
    """Fold repeatedly, keep the pose with the highest ipTM. No ground truth involved."""
    from proto_tools.mcp.tools import run_tool

    dest.mkdir(parents=True, exist_ok=True)
    best = None
    for seed in SEEDS:
        seed_dir = dest / f"seed{seed}"
        seed_dir.mkdir(exist_ok=True)
        result = run_tool(
            "boltz2-prediction",
            inputs={"complexes": [{"chains": [
                {"id": "B", "sequence": spec.binder.sequence, "entity_type": "protein"},
                {"id": "T", "sequence": spec.target.sequence, "entity_type": "protein"},
            ]}]},
            config={"seed": seed, "recycling_steps": RECYCLING_STEPS,
                    "use_msa": True, "diffusion_samples": DIFFUSION_SAMPLES},
            device="modal",
            output_dir=str(seed_dir),
        )
        if not result.get("ok", True):
            print(f"    seed {seed}: FAILED {str(result.get('error'))[:90]}", flush=True)
            continue
        for i, structure in enumerate(result["result"]["structures"]):
            iptm = float(structure["metrics"]["iptm"])
            print(f"    seed {seed} sample {i}: ipTM {iptm:.3f}", flush=True)
            if best is None or iptm > best["iptm"]:
                written = sorted(seed_dir.glob("*structure*"))
                if written:
                    keep = dest / "best.cif"
                    shutil.copyfile(written[min(i, len(written) - 1)], keep)
                    best = {"iptm": iptm, "seed": seed, "sample": i, "path": keep}
    return best


def main() -> int:
    cases = sys.argv[1:] or DEFAULT_CASES
    specs = {s.case_id: s for s in caseset.tier_c()[0]}
    baseline = {r["case_id"]: r for r in json.loads((OUT / "g1_dockq.json").read_text())["per_case"]}

    n_folds = len(cases) * len(SEEDS)
    print(f"{len(cases)} cases x {len(SEEDS)} seeds x {DIFFUSION_SAMPLES} samples, "
          f"{RECYCLING_STEPS} recycling steps")
    print(f"~{n_folds} trunk passes; selection by ipTM only\n", flush=True)

    rows = []
    for cid in cases:
        print(f"=== {cid} ===", flush=True)
        best = best_pose(specs[cid], WORK / cid)
        if best is None:
            rows.append({"case_id": cid, "error": "no pose produced"})
            continue
        row = {
            "case_id": cid,
            "baseline_DockQ": round(baseline.get(cid, {}).get("DockQ") or 0.0, 4),
            "baseline_iptm": None,
            "upgraded_iptm": round(best["iptm"], 4),
            "chosen_seed": best["seed"], "chosen_sample": best["sample"],
            "structure": str(best["path"].relative_to(REPO)),
        }
        rows.append(row)
        print(f"  best ipTM {best['iptm']:.3f} (seed {best['seed']}, sample {best['sample']})\n", flush=True)

    payload = {
        "question": "is the Tier C epitope failure a sampling limit or a knowledge limit?",
        "settings": {"seeds": SEEDS, "diffusion_samples": DIFFUSION_SAMPLES,
                     "recycling_steps": RECYCLING_STEPS, "selection": "highest ipTM"},
        "leakage_guard": "poses are ranked by the model's own ipTM. DockQ is computed only after "
                         "the pose is fixed, and never used to choose one.",
        "per_case": rows,
        "next": "score logs/g1_upgraded/*/best.cif against the deposits with g1_dockq-style DockQ",
    }
    (OUT / "g1_upgrade.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUT / 'g1_upgrade.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
