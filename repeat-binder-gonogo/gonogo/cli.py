"""`gonogo` — the go/no-go control run, driven in phases.

    python -m gonogo preflight              what can run right now, and what it costs
    python -m gonogo inputs                 resolve every case to sequences, check the literature claims
    python -m gonogo run --scorer boltz2     the real run, in cost order, with the Referee's veto live
    python -m gonogo run --scorer synthetic --allow-synthetic --out demo/selftest
    python -m gonogo verify [--predictions PATH]   run the bench's own verifier

The phase split is the point. Everything before `run` is free and reveals
whether the expensive part is worth starting.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

from . import caseset, gates, registry
from .caseset import ComplexSpec
from .scorers import SCORERS, CachedScorer, Scorer, Unavailable
from .trajectory import Scribe, default_scribe

REPO = Path(__file__).resolve().parent.parent
TASK_DIR = REPO / "bench" / "repeat-binder-gonogo"
DEFAULT_OUT = TASK_DIR / "results"

TIER_B_POSITIVES = ("B1", "B2", "B3")


# --------------------------------------------------------------------------
# input resolution + the checks that make the case set trustworthy


def resolve_inputs(scribe: Scribe, with_tier_a: bool, tier_a_scaffold: str) -> tuple[list[ComplexSpec], list[caseset.Gap]]:
    specs: list[ComplexSpec] = []
    all_gaps: list[caseset.Gap] = []

    tier_b, gaps_b = caseset.tier_b()
    specs += tier_b
    all_gaps += gaps_b
    scribe.append("librarian", "resolve tier B", f"{len(tier_b)} complexes from PDB + UniProt, CPU only, no billing",
                  evidence=[s.case_id for s in tier_b])

    tier_c, gaps_c = caseset.tier_c()
    specs += tier_c
    all_gaps += gaps_c
    scribe.append("librarian", "resolve tier C", f"{len(tier_c)} deposited complexes across 5 scaffold families",
                  evidence=[s.case_id for s in tier_c])

    if with_tier_a:
        tier_a, gaps_a = caseset.tier_a(tier_a_scaffold)
        specs += tier_a
        all_gaps += gaps_a
        scribe.append(
            "librarian",
            "resolve tier A",
            f"{len(tier_a)} ladder variants built on {tier_a_scaffold}; K_D from Table S1 of PMC12202725 "
            f"(the TODO in cases.yaml, now closed)",
            evidence={"n": len(tier_a), "source": "gonogo/data/tier_a_kd.yaml"},
        )

    decoys = caseset.decoys(specs)
    specs += decoys
    scribe.append("architect", "build decoys", f"{len(decoys)} decoys: composition-preserving shuffles and "
                  f"non-cognate pairs drawn from unrelated systems", evidence=[s.case_id for s in decoys])

    for gap in all_gaps:
        scribe.append("referee", "input gap", f"{gap.what}: {gap.why}", outcome="gap", evidence={"blocks": gap.blocks})
    return specs, all_gaps


def check_case_set(scribe: Scribe) -> list[dict[str, Any]]:
    """Re-derive the case set's own claims before trusting any of them."""
    checks: list[dict[str, Any]] = []

    mouse, _ = registry.uniprot_region(caseset.MOUSE_CATB, "Chain", "Cathepsin B")
    human, _ = registry.uniprot_region(caseset.HUMAN_CATB, "Chain", "Cathepsin B")
    observed = registry.identity(mouse, human)
    ok = 0.80 <= observed <= 0.86
    checks.append({"check": "mouse/human cathepsin B identity", "claimed": "83%",
                   "observed": f"{observed:.2%}", "ok": ok})
    scribe.append("referee", "verify identity claim",
                  f"cases.yaml says 83% identical; global alignment of the two mature chains gives {observed:.2%}",
                  outcome="ok" if ok else "blocked")

    residues = {"mouse": (mouse[64], mouse[65]), "human": (human[64], human[65])}
    ok = residues["mouse"] == ("I", "Q") and residues["human"] == ("S", "M")
    checks.append({"check": "selectivity residues 65/66", "claimed": "mouse I65/Q66, human S65/M66",
                   "observed": f"mouse {''.join(residues['mouse'])}, human {''.join(residues['human'])}", "ok": ok})
    scribe.append("referee", "verify selectivity residues",
                  f"mature-chain positions 65/66 read mouse {''.join(residues['mouse'])} / "
                  f"human {''.join(residues['human'])}", outcome="ok" if ok else "blocked")

    deposited = caseset.deposited_chain("9S60", "target")
    ok = deposited.sequence == mouse
    checks.append({"check": "B1 target construct", "claimed": "deposited 9S60 target is the mature mouse enzyme",
                   "observed": "identical" if ok else "differs", "ok": ok})
    scribe.append("referee", "verify construct boundaries",
                  "B4 must use the same construct boundaries as B1 or the gate measures the construct, not the "
                  f"ortholog: deposited 9S60 target vs UniProt mature chain — {'identical' if ok else 'DIFFERS'}",
                  outcome="ok" if ok else "blocked")

    for pdb_id, scaffold in (("5AEI", "scaffold_5aei"), ("6SA8", "scaffold_6sa8")):
        specs, gaps = caseset.tier_a(scaffold, pockets=["arg"])
        ok = bool(specs) and not gaps
        checks.append({"check": f"tier A pocket mapping in {pdb_id}",
                       "claimed": "pocket positions from the paper read the Arg-binder pocket QWSQQEW",
                       "observed": "confirmed" if ok else (gaps[0].why if gaps else "no variants built"), "ok": ok})
        scribe.append("librarian", f"verify pocket mapping {pdb_id}",
                      "literature residue numbers are author numbering, which is not the sequence index; "
                      f"checked against the deposited coordinates — {'confirmed' if ok else 'MISMATCH'}",
                      outcome="ok" if ok else "blocked")
    return checks


# --------------------------------------------------------------------------
# phases


def cmd_preflight(args: argparse.Namespace) -> int:
    out = Path(args.out)
    scribe = default_scribe(out)
    print("\n== preflight: what can run, and what it costs ==\n")

    specs, gaps = resolve_inputs(scribe, with_tier_a=args.tier_a, tier_a_scaffold=args.scaffold)
    counts: dict[str, int] = {}
    for spec in specs:
        counts[spec.tier] = counts.get(spec.tier, 0) + 1
    print(f"\ncases resolved: {counts}")

    scorer = SCORERS[args.scorer]()
    try:
        scorer.preflight()
        available, reason = True, "ready"
    except Unavailable as exc:
        available, reason = False, str(exc)
    except Exception as exc:  # a missing backend must not look like a passing gate
        available, reason = False, f"{type(exc).__name__}: {exc}"

    scribe.append("economist", "scorer preflight", reason, outcome="ok" if available else "blocked",
                  evidence={"scorer": scorer.key})

    print(f"\nscorer {scorer.key}: {'AVAILABLE' if available else 'UNAVAILABLE'}")
    if not available:
        print(f"  reason: {reason}")

    print("\ncost, in the order the task requires (cheap and discriminating first).")
    print("Each row counts only the complexes that gate is the first to need; a later gate reuses what an")
    print("earlier one already paid for, which is why G3 can cost nothing once G1 has run over Tier A:")
    plan = [
        ("G5", "decoy rejection", [s for s in specs if s.tier == "decoy"] + [s for s in specs if s.case_id in TIER_B_POSITIVES]),
        ("G2", "specificity discrimination", [s for s in specs if s.tier == "B"]),
        ("G1", "interface recovery", [s for s in specs if s.tier in {"A", "C"}]),
        ("G3", "affinity ranking", [s for s in specs if s.tier == "A"]),
    ]
    total = 0.0
    already: set[str] = set()
    for gate_id, name, subset in plan:
        # each complex is bought once; a later gate reuses what an earlier one paid for
        subset = [spec for spec in subset if spec.case_id not in already]
        already.update(spec.case_id for spec in subset)
        estimate = scorer.cost_estimate(len(subset))
        usd = estimate.get("usd") or 0.0
        total += usd
        print(f"  {gate_id} {name:28s} {len(subset):3d} complexes  ~{estimate.get('gpu_seconds', 0)/60:6.1f} GPU-min  ~${usd}")
    print(f"  {'':36s}{'':3s}{'':22s}~${total:.2f} if every tier runs")
    print("\nnote: G2 is 6 complexes and vetoes everything after it. Nothing past G2 is worth buying until it passes.")

    if gaps:
        print("\nunresolved inputs:")
        for gap in gaps:
            print(f"  ? {gap.what}: {gap.why} (blocks {', '.join(gap.blocks)})")

    print(f"\ntrajectory: {scribe.path}")
    return 0 if available else 3


def cmd_inputs(args: argparse.Namespace) -> int:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    scribe = default_scribe(out)
    print("\n== inputs: resolve every case, then re-derive the claims made about it ==\n")

    specs, gaps = resolve_inputs(scribe, with_tier_a=args.tier_a, tier_a_scaffold=args.scaffold)
    print()
    checks = check_case_set(scribe)

    print("\ncase set checks:")
    for check in checks:
        print(f"  {'PASS' if check['ok'] else 'FAIL'}  {check['check']}: claimed {check['claimed']} -> observed {check['observed']}")

    payload = {
        "cases": [spec.as_record() for spec in specs],
        "gaps": [asdict(gap) for gap in gaps],
        "checks": checks,
    }
    (out / "inputs.json").write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {out / 'inputs.json'} ({len(specs)} complexes)")
    return 0 if all(check["ok"] for check in checks) else 4


def cmd_run(args: argparse.Namespace) -> int:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    scribe = default_scribe(out)
    print("\n== run: gates in cost order, Referee veto live ==\n")

    scorer_class = SCORERS[args.scorer]
    scorer: Scorer = scorer_class(seeds=args.seeds) if args.seeds else scorer_class()
    if not args.no_cache:
        scorer = CachedScorer(scorer, out / "score_cache")
    if args.scorer == "synthetic":
        if not args.allow_synthetic:
            print("refusing: the synthetic scorer produces meaningless numbers. Pass --allow-synthetic "
                  "to run it as a harness self-test.", file=sys.stderr)
            return 2
        if Path(args.out).resolve() == DEFAULT_OUT.resolve():
            print(f"refusing: synthetic output may not be written to the bench submission path "
                  f"({DEFAULT_OUT}). Use --out demo/selftest.", file=sys.stderr)
            return 2

    if args.ignore_veto and Path(args.out).resolve() == DEFAULT_OUT.resolve():
        print(f"refusing: --ignore-veto violates the task's cost policy, so its output may not be written to "
              f"the bench submission path ({DEFAULT_OUT}). Use --out demo/full-sweep.", file=sys.stderr)
        return 2

    specs, gaps = resolve_inputs(scribe, with_tier_a=args.tier_a, tier_a_scaffold=args.scaffold)
    by_id = {spec.case_id: spec for spec in specs}

    notes: list[str] = [
        f"score convention: higher is better; {scorer.score_definition}. No energies are submitted, so no sign flip.",
        "uncertainty band = spread of the score across independent seeds of the same model. It measures run-to-run "
        "instability, not model error; it is the only band a single model can honestly report.",
    ]
    for gap in gaps:
        notes.append(f"input gap: {gap.what} — {gap.why} (blocks {', '.join(gap.blocks)})")
    if args.scorer == "synthetic":
        notes.insert(0, "SYNTHETIC HARNESS SELF-TEST. Every number below is a hash of the input sequences and has "
                        "no physical meaning. This file exists to exercise the gate arithmetic and the verifier "
                        "wiring; it is not a scientific result.")

    try:
        scorer.preflight()
    except Unavailable as exc:
        scribe.append("referee", "abort before spending", str(exc), outcome="blocked")
        verdicts = [
            gates.Verdict("G5", "not_run", f"scorer unavailable: {exc}"),
            gates.Verdict("G2", "not_run", f"scorer unavailable: {exc}"),
            gates.g1_interface_recovery([], structures_available=False),
            gates.g3_affinity_ranking([], specs),
            gates.g4_calibration([], specs, decision_threshold=None),
        ]
        notes.append(f"no gate ran: {exc}")
        write_submission(out, [], verdicts, notes, scribe)
        print("\nnothing was scored. The submission records every gate as not_run with the reason above.")
        return 3

    # ---- G5, cheapest, first
    g5_specs = [spec for spec in specs if spec.tier == "decoy"] + [by_id[cid] for cid in TIER_B_POSITIVES if cid in by_id]
    scribe.append("economist", "buy G5", f"{len(g5_specs)} complexes, the cheapest discriminating gate",
                  cost=scorer.cost_estimate(len(g5_specs)))
    all_scores = list(scorer.score(g5_specs))
    g5 = g5_decoy = gates.g5_decoy_rejection(all_scores, TIER_B_POSITIVES)
    scribe.append("referee", "G5", g5.evidence, outcome="ok" if g5.status == "pass" else "veto")

    decoy_scores = [s.score for s in all_scores if s.case_id.startswith("decoy_")]
    threshold = max(decoy_scores) if decoy_scores else None

    # ---- G2, the gate the whole exercise exists for
    remaining_b = [spec for spec in specs if spec.tier == "B" and spec.case_id not in {s.case_id for s in all_scores}]
    scribe.append("economist", "buy G2", f"{len(remaining_b)} further complexes; failing this ends the run",
                  cost=scorer.cost_estimate(len(remaining_b)))
    all_scores += list(scorer.score(remaining_b))
    g2 = gates.g2_specificity(all_scores)
    scribe.append("referee", "G2", g2.evidence, outcome="ok" if g2.status == "pass" else "veto")

    verdicts = [g5, g2]
    if g2.status != "pass" and not args.ignore_veto:
        scribe.append("referee", "veto", "G2 did not pass; the remaining tiers are not bought", outcome="veto")
        notes.append("G1/G3/G4 were not run because G2 did not pass. The task's cost policy makes this the "
                     "required behaviour, not an omission.")
        verdicts += [
            gates.Verdict("G1", "not_run", "not run: G2 did not pass and the cost policy stops the run there"),
            gates.Verdict("G3", "not_run", "not run: G2 did not pass and the cost policy stops the run there"),
            gates.Verdict("G4", "not_run", "not run: G2 did not pass and the cost policy stops the run there"),
        ]
        write_submission(out, all_scores, verdicts, notes, scribe)
        return 1

    if g2.status != "pass":
        # The veto is the task's rule, not the harness's preference. Overriding it produces
        # a demonstration of the downstream gates, not a submission — say so loudly, in the
        # file itself, so the artefact cannot be mistaken for a compliant run later.
        scribe.append("referee", "veto overridden", "G2 did not pass; buying the remaining tiers anyway "
                      "by explicit request (--ignore-veto). This is a demo sweep, not a submission.",
                      outcome="veto")
        notes.append("NOT A VALID SUBMISSION. G2 failed and the remaining tiers were bought anyway by explicit "
                     "request (--ignore-veto), to exercise G1/G3/G4 for demonstration. task.md requires the run "
                     "to stop at a failed G2 — 'If G2 fails, stop and report'. Every gate verdict below G2 is "
                     "therefore reported against a scorer that already failed the gate qualifying it, and should "
                     "be read as a capability demo of the harness, not as evidence about the pipeline.")

    # ---- everything past here is only bought because G2 passed
    tier_ac = [spec for spec in specs if spec.tier in {"A", "C"}]
    if args.stop_after_g2:
        notes.append("run stopped after G2 by request (--stop-after-g2); G1/G3/G4 not attempted")
        verdicts += [
            gates.Verdict("G1", "not_run", "stopped after G2 by request"),
            gates.Verdict("G3", "not_run", "stopped after G2 by request"),
            gates.g4_calibration(all_scores, specs, threshold),
        ]
    else:
        justification = (
            f"justified by G2 passing with margin {g2.detail.get('margin')}"
            if g2.status == "pass"
            else f"NOT justified by a gate: G2 returned {g2.status} (margin {g2.detail.get('margin')}) and the "
                 f"veto was overridden by request"
        )
        scribe.append("economist", "buy tiers A and C",
                      f"{len(tier_ac)} complexes, {justification}", cost=scorer.cost_estimate(len(tier_ac)))
        all_scores += list(scorer.score(tier_ac))
        verdicts += [
            gates.g1_interface_recovery(all_scores, structures_available=args.scorer != "synthetic"),
            gates.g3_affinity_ranking(all_scores, specs),
            gates.g4_calibration(all_scores, specs, threshold),
        ]
    for verdict in verdicts[2:]:
        scribe.append("referee", verdict.gate_id, verdict.evidence,
                      outcome="ok" if verdict.status == "pass" else ("gap" if verdict.status == "not_run" else "veto"))

    if isinstance(scorer, CachedScorer):
        notes.append(f"checkpointing: {scorer.hits} complexes reused from {out / 'score_cache'}, "
                     f"{scorer.misses} newly computed")
    write_submission(out, all_scores, verdicts, notes, scribe)
    return 0


def write_submission(out: Path, scores, verdicts, notes: list[str], scribe: Scribe) -> Path:
    payload = {
        "predictions": [score.as_prediction() for score in scores],
        "gates": [verdict.as_record() for verdict in verdicts],
        "notes": " | ".join(notes),
    }
    out.mkdir(parents=True, exist_ok=True)
    path = out / "predictions.json"
    path.write_text(json.dumps(payload, indent=2))
    scribe.append("scribe", "write submission",
                  f"{len(payload['predictions'])} predictions, {len(payload['gates'])} gates",
                  evidence=str(path))
    print(f"\nwrote {path}")
    for verdict in verdicts:
        print(f"  {verdict.gate_id}: {verdict.status}")
    return path


def cmd_verify(args: argparse.Namespace) -> int:
    predictions = Path(args.predictions or (Path(args.out) / "predictions.json"))
    print(f"\n== verify: running the bench's own verifier against {predictions} ==\n")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(TASK_DIR / "verifier" / "test_outputs.py"), "-p", "no:randomly", "-q"],
        env={**dict(__import__("os").environ), "PREDICTIONS_PATH": str(predictions)},
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    print(result.stdout[-4000:])
    if result.stderr.strip():
        print(result.stderr[-2000:], file=sys.stderr)
    return result.returncode


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gonogo", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output directory (default: the bench results dir)")
    parser.add_argument("--scorer", default="boltz2", choices=sorted(SCORERS), help="which scorer to use")
    parser.add_argument("--tier-a", action="store_true", help="include the Tier A affinity ladder (63 variants/scaffold)")
    parser.add_argument("--scaffold", default="scaffold_5aei", choices=["scaffold_5aei", "scaffold_6sa8"])
    parser.add_argument(
        "--seeds", type=int, nargs="+", metavar="N",
        help="seeds to average over (default: 0 1 2). More seeds is a different cache key, so it "
             "re-scores every complex. Note that G2's band is a sample range, whose expectation grows "
             "with seed count — see docs/gonogo-feasibility.md §4 before reading a wider band as a "
             "worse model.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("preflight", help="what can run right now, and what it costs").set_defaults(func=cmd_preflight)
    sub.add_parser("inputs", help="resolve cases and re-derive the claims made about them").set_defaults(func=cmd_inputs)

    run = sub.add_parser("run", help="run the gates in cost order")
    run.add_argument("--allow-synthetic", action="store_true", help="permit the meaningless self-test scorer")
    run.add_argument("--stop-after-g2", action="store_true", help="buy nothing past the specificity gate")
    run.add_argument("--ignore-veto", action="store_true",
                     help="run G1/G3/G4 even when G2 fails. Violates the task's cost policy, so the output is "
                          "stamped NOT A VALID SUBMISSION and may not be written to the bench results dir.")
    run.add_argument("--no-cache", action="store_true", help="re-run every complex instead of reusing checkpoints")
    run.set_defaults(func=cmd_run)

    verify = sub.add_parser("verify", help="run the bench verifier on a submission")
    verify.add_argument("--predictions", help="path to predictions.json")
    verify.set_defaults(func=cmd_verify)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
