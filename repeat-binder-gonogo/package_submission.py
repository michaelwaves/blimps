"""Assemble the hackathon submission: one ordered trajectory plus its artifacts.

    .venv/bin/python package_submission.py

WHAT "CLEAN" MEANS HERE, AND WHAT IT DOES NOT

The trajectory rubric checks honest reporting as a pass/fail gate before it
scores anything, and D4 requires failed attempts and dead ends to stay in the
JSONL rather than be cleaned out. So this does not prune. Every abort, veto and
crash stays.

What it fixes is that the file is seventeen separate invocations concatenated:
`t` restarts on each one, and thirty lines are byte-identical repeats of the
per-attempt setup. Those repeats are evidence that setup re-ran, not noise, so
they are segmented rather than deleted — each record gains a run_id, a phase,
and a globally monotonic sequence number. A reader can then follow one ordered
narrative, or filter to a single attempt, without anything having been removed.

It also backfills the GvpA design work, which ran outside the harness and so
never reached a scribe. Those records are marked `logged: retrospective` so they
are not mistaken for live instrumentation.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent
OUT = REPO / "submission"

SOURCES = [
    ("control-run", REPO / "bench" / "repeat-binder-gonogo" / "results" / "trajectory.jsonl"),
    ("full-sweep-demo", REPO / "demo" / "full-sweep" / "trajectory.jsonl"),
    ("harness-selftest", REPO / "demo" / "selftest" / "trajectory.jsonl"),
]

# The GvpA design run happened through standalone scripts, outside the scribe.
# Reconstructed from the artifacts each step wrote, and labelled as such.
GVPA_STEPS = [
    ("librarian", "resolve design target", "GvpA_surface_patch_5ribs_x11.pdb: 55 subunits, 3575 residues, "
     "266 A across. 8GBS.cif: 4 GvpA + 1 GvpC, the native Ana shell model.", "ok",
     "designs/gvpa/site.json"),
    ("librarian", "superimpose 8GBS onto the patch", "the two inputs share no coordinate frame; aligned a GvpA "
     "subunit (CA RMSD 1.76 A over 65 atoms) and carried GvpC across with the same transform", "ok",
     "designs/gvpa/gvpc_reference.pdb"),
    ("referee", "assess the GvpC reference", "GvpC in 8GBS is a backbone-only poly-UNK trace: 5 atoms/residue, "
     "no residue identities, and only 3 GvpA residues within 5 A. It localises the natural binder but is NOT a "
     "resolved interface, so DockQ against it would be meaningless. Similarity must be measured as spatial "
     "footprint overlap instead.", "blocked", "designs/gvpa/site.json"),
    ("architect", "trim the target", "cut 3575 residues to the 5 subunits GvpC lies across (325 residues); a "
     "neighbour shell pushed it to 1690, past what a binder-design pipeline accepts. Relabelled chains l-p to "
     "A-E because contig and hotspot parsers reject lowercase ids.", "ok", "designs/gvpa/gvpa_target.pdb"),
    ("economist", "request rfdiffusion3 deploy", "de novo backbone generation needs rfdiffusion3, which is "
     "deployable but not deployed. Deployment declined three times at the approval prompt; workspace reports "
     "authenticated, environment present, deployable true, so this is not a credits or config failure.",
     "blocked", None),
    ("architect", "fall back to inverse folding", "without a backbone generator this is no longer de novo "
     "design: the GvpC backbone placement is reused as the scaffold and only the sequence is designed. "
     "Recorded as inverse folding, not as de novo binder design.", "ok", None),
    ("economist", "buy proteinmpnn-sample", "32 sequences for chain F against a fixed GvpA surface", "ok",
     "designs/gvpa/designed_sequences.json"),
    ("referee", "score native recovery", "mean recovery 11.8% against native GvpC 85-117 (P09413). Control: the "
     "same sequences shuffled recover 11.4%. Positional signal above composition is +0.4 pp -- effectively "
     "zero. Consistent with the 3-residue contact: there was no interface context to design against.", "veto",
     "designs/gvpa/designed_sequences.json"),
    ("referee", "curvature check, first attempt", "swept cylinder axes and fitted circles to the projection; "
     "returned radius 63 A over a 137 deg arc, reported as well constrained. WRONG -- PyMOL rendering showed "
     "the patch is a flat sheet, and a 137 deg arc would look like a letter C. The fit had locked onto "
     "internal scatter rather than global curvature.", "blocked", "designs/gvpa/patch_side.png"),
    ("referee", "curvature check, corrected", "measured sagitta of subunit centroids along each rib instead; "
     "R = L^2/8s, no axis search to overfit. All five ribs agree exactly: chord 121.3 A, sagitta 10.35 A, "
     "implied diameter 35.6 nm.", "ok", "designs/gvpa/curvature.json"),
    ("librarian", "literature comparison via paperclip", "Ana GVs measure 85 +/- 4 nm in diameter with a ~2.8 nm "
     "shell (Dutka et al. 2023, PMC10185304 L22 and L83) -- the source study for 8GBS. The generated patch is "
     "2.4x over-curved for its own species.", "ok", "designs/gvpa/curvature.json"),
]


def load(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]


def main() -> int:
    OUT.mkdir(exist_ok=True)
    merged: list[dict] = []
    seq = 0
    summary = []

    for phase, path in SOURCES:
        recs = load(path)
        # A drop in t marks the start of a new invocation.
        run_idx, prev_t = 0, -1.0
        for r in recs:
            t = float(r.get("t", 0.0))
            if t < prev_t:
                run_idx += 1
            prev_t = t
            seq += 1
            merged.append({
                "seq": seq,
                "phase": phase,
                "run_id": f"{phase}#{run_idx}",
                **r,
            })
        summary.append((phase, len(recs), run_idx + 1 if recs else 0))

    # Backfill the design run.
    for i, (role, step, rationale, outcome, evidence) in enumerate(GVPA_STEPS):
        seq += 1
        merged.append({
            "seq": seq, "phase": "gvpa-design", "run_id": "gvpa-design#0",
            "t": float(i), "role": role, "step": step, "rationale": rationale,
            "outcome": outcome, "cost": {"gpu_seconds": 0, "usd": 0.0},
            "evidence": evidence, "logged": "retrospective",
        })
    summary.append(("gvpa-design", len(GVPA_STEPS), 1))

    traj = OUT / "trajectory.jsonl"
    with traj.open("w") as fh:
        for r in merged:
            fh.write(json.dumps(r) + "\n")

    # ---- validate what was written, rather than assuming it -----------------
    written = [json.loads(l) for l in traj.read_text().splitlines() if l.strip()]
    seqs = [r["seq"] for r in written]
    outcomes = {}
    for r in written:
        outcomes[r.get("outcome")] = outcomes.get(r.get("outcome"), 0) + 1
    checks = {
        "every line parses as JSON": len(written) == len(merged),
        "seq is strictly increasing": seqs == sorted(seqs) and len(set(seqs)) == len(seqs),
        "every record has role/step/outcome": all(
            r.get("role") and r.get("step") and r.get("outcome") for r in written),
        "failures retained (blocked/veto present)": (outcomes.get("blocked", 0) + outcomes.get("veto", 0)) > 0,
    }

    print(f"{'phase':<20} {'records':>8} {'invocations':>12}")
    print("-" * 42)
    for phase, n, runs in summary:
        print(f"{phase:<20} {n:>8} {runs:>12}")
    print(f"{'TOTAL':<20} {len(written):>8}")
    print(f"\noutcomes: {outcomes}")
    print("\nvalidation:")
    for name, ok in checks.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")

    # ---- copy the artifacts every claim points at ---------------------------
    artifacts = OUT / "artifacts"
    artifacts.mkdir(exist_ok=True)
    for rel in ["bench/repeat-binder-gonogo/results/predictions.json",
                "bench/repeat-binder-gonogo/results/predictions.n3-range-band.json",
                "bench/repeat-binder-gonogo/results/g1_dockq.json",
                "bench/repeat-binder-gonogo/results/g1_diagnosis.json",
                "bench/repeat-binder-gonogo/results/g1_upgrade_scored.json",
                "bench/repeat-binder-gonogo/results/g2_seed_depth.json",
                "bench/repeat-binder-gonogo/results/inputs.json",
                "bench/repeat-binder-gonogo/cases.yaml",
                "bench/repeat-binder-gonogo/REFERENCES.md",
                "designs/gvpa/curvature.json",
                "designs/gvpa/site.json",
                "designs/gvpa/designed_sequences.json"]:
        src = REPO / rel
        if src.exists():
            shutil.copyfile(src, artifacts / src.name)

    manifest = {
        "submission": "repeat-binder-gonogo + GvpA design run",
        "trajectory": {"path": "trajectory.jsonl", "records": len(written),
                       "phases": {p: n for p, n, _ in summary},
                       "outcomes": outcomes,
                       "note": "seventeen invocations of the control run are segmented by run_id rather "
                               "than merged or pruned; t restarts within each. Nothing was removed."},
        "validation": checks,
        "artifacts": sorted(p.name for p in artifacts.iterdir()),
        "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                   for p in sorted(artifacts.iterdir())},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nwrote {traj} and {OUT / 'manifest.json'}")
    print(f"artifacts copied: {len(list(artifacts.iterdir()))}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
