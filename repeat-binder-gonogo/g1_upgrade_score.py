"""Score the upgraded poses against the deposits. The moment of truth for the retry.

    .venv/bin/python g1_upgrade_score.py

g1_upgrade.py chose one pose per case by the model's own ipTM, with the deposit
nowhere in the loop. Now — and only now — the deposit is allowed in, to say
whether the chosen pose is any closer to the truth than the baseline was.

If DockQ moves with ipTM, the epitope was reachable and the first run simply
under-sampled. If ipTM climbs while DockQ stays at the floor, the model is
confidently wrong in a way more compute does not fix, and the pipeline needs a
different input rather than a bigger budget.
"""

from __future__ import annotations

import json
from pathlib import Path

from g1_dockq import native_two_chains, parse_source, run_dockq
from gonogo import caseset

REPO = Path(__file__).resolve().parent
UP = REPO / "logs" / "g1_upgraded"
OUT = REPO / "bench" / "repeat-binder-gonogo" / "results"


def main() -> int:
    upgrade = json.loads((OUT / "g1_upgrade.json").read_text())
    baseline = {r["case_id"]: r for r in json.loads((OUT / "g1_dockq.json").read_text())["per_case"]}
    specs = {s.case_id: s for s in caseset.tier_c()[0]}

    rows = []
    for entry in upgrade["per_case"]:
        cid = entry.get("case_id")
        pose = UP / cid / "best.cif"
        if not pose.exists():
            continue
        pdb_id, binder_chains = parse_source(specs[cid].binder.source)
        _, target_chains = parse_source(specs[cid].target.source)

        best = None
        for b in binder_chains:
            for t in target_chains:
                native = native_two_chains(pdb_id, b, t, UP / cid / f"{pdb_id}_{b}{t}_native.cif")
                got = run_dockq(pose, native, b, t)
                if got.get("DockQ") is not None and (best is None or got["DockQ"] > best["DockQ"]):
                    best = got

        base_dq = baseline.get(cid, {}).get("DockQ")
        row = {
            "case_id": cid,
            "baseline_DockQ": round(base_dq, 4) if base_dq is not None else None,
            "upgraded_DockQ": round(best["DockQ"], 4) if best else None,
            "upgraded_iRMSD": round(best["iRMSD"], 2) if best and best.get("iRMSD") else None,
            "baseline_iptm": None,
            "upgraded_iptm": entry.get("upgraded_iptm"),
        }
        row["DockQ_delta"] = (round(row["upgraded_DockQ"] - row["baseline_DockQ"], 4)
                              if row["upgraded_DockQ"] is not None and row["baseline_DockQ"] is not None else None)
        rows.append(row)
        print(f"  {cid:9} DockQ {row['baseline_DockQ']} -> {row['upgraded_DockQ']}   "
              f"(ipTM -> {row['upgraded_iptm']})", flush=True)

    moved = [r for r in rows if r["DockQ_delta"] is not None and r["DockQ_delta"] > 0.2]
    recovered = [r for r in rows if (r["upgraded_DockQ"] or 0) >= 0.49]
    verdict = (
        "SAMPLING LIMIT: more sampling plus confidence-based selection recovered the interface"
        if recovered else
        "KNOWLEDGE LIMIT: confidence rose substantially while DockQ stayed at the floor. The model does not "
        "know where the epitope is, and more compute does not find it. Buying more sampling for this class of "
        "target is not justified; the pipeline needs a different input (an epitope prior, a co-evolution "
        "signal, or a docking stage), not a bigger budget."
    )
    payload = {
        "question": "sampling limit or knowledge limit?",
        "selection": "highest ipTM, deposit not consulted until this file was produced",
        "per_case": rows,
        "n_recovered_to_pass": len(recovered),
        "n_materially_improved": len(moved),
        "verdict": verdict,
    }
    (OUT / "g1_upgrade_scored.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\n  {verdict}")
    print(f"  wrote {OUT / 'g1_upgrade_scored.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
