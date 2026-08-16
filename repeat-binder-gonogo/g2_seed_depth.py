"""How deep does G2's instability go? Extend B1 and B4 to more seeds.

    .venv/bin/python g2_seed_depth.py [--to 12]

The shipped run scored every complex at seeds 0,1,2. Those three replicates are
already paid for and cached, so this buys only the *new* seeds for the two cases
G2 actually turns on, and appends them to what is already there.

It writes a diagnostic, not a submission. `bench/.../results/predictions.json`
is not touched: its verdict was issued under a pre-declared band and re-issuing
it against a band chosen after seeing the numbers is the tuning this task
guards against everywhere else. The point here is to answer a different
question — is B4's spread real, or an artefact of n=3?
"""

from __future__ import annotations

import argparse
import json
import math
import statistics as st
from pathlib import Path

from gonogo.cli import resolve_inputs
from gonogo.trajectory import Scribe

REPO = Path(__file__).resolve().parent
RESULTS = REPO / "bench" / "repeat-binder-gonogo" / "results"
CASES = ("B1", "B4")


def cached_replicates(case_id: str) -> list[float]:
    """The seeds 0,1,2 the shipped run already bought."""
    for path in (RESULTS / "score_cache").glob(f"{case_id}.*.json"):
        return list(json.loads(path.read_text())["replicates"])
    raise SystemExit(f"no cached score for {case_id}; run `gonogo run` first")


def score_one_seed(spec, seed: int) -> float:
    from proto_tools.mcp.tools import run_tool

    result = run_tool(
        "boltz2-prediction",
        inputs={"complexes": [{"chains": [
            {"id": "B", "sequence": spec.binder.sequence, "entity_type": "protein"},
            {"id": "T", "sequence": spec.target.sequence, "entity_type": "protein"},
        ]}]},
        config={"seed": seed, "recycling_steps": 3, "use_msa": True, "diffusion_samples": 1},
        device="modal",
    )
    if not result.get("ok", True):
        raise RuntimeError(f"{spec.case_id} seed {seed}: {result.get('error')}")
    return 100.0 * float(result["result"]["structures"][0]["metrics"]["pair_chains_iptm"][0][1])


def welch(a: list[float], b: list[float]) -> tuple[float, float, float]:
    """Difference of means, its standard error, and Welch's degrees of freedom."""
    sea = st.stdev(a) / math.sqrt(len(a))
    seb = st.stdev(b) / math.sqrt(len(b))
    se = math.sqrt(sea**2 + seb**2)
    df = (sea**2 + seb**2) ** 2 / (sea**4 / (len(a) - 1) + seb**4 / (len(b) - 1))
    return st.mean(a) - st.mean(b), se, df


# t critical values at 95%, indexed by rounded df
T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31,
       9: 2.26, 10: 2.23, 12: 2.18, 15: 2.13, 20: 2.09, 30: 2.04}


def t_crit(df: float) -> float:
    return T95[min(T95, key=lambda k: abs(k - df))]


def verdicts(b1: list[float], b4: list[float]) -> dict:
    margin = st.mean(b1) - st.mean(b4)
    rng = max(max(b1) - min(b1), max(b4) - min(b4))
    diff, se, df = welch(b1, b4)
    ci_half = t_crit(df) * se
    return {
        "n": len(b1),
        "mean_B1": round(st.mean(b1), 3), "mean_B4": round(st.mean(b4), 3),
        "sd_B1": round(st.stdev(b1), 3), "sd_B4": round(st.stdev(b4), 3),
        "margin": round(margin, 3),
        "band_range": round(rng, 3), "verdict_range": "pass" if margin > rng else "fail",
        "band_ci_half": round(ci_half, 3), "verdict_ci": "pass" if margin > ci_half else "fail",
        "welch_t": round(diff / se, 3), "welch_df": round(df, 2),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--to", type=int, default=12, help="total seeds per case (default 12)")
    args = ap.parse_args()

    scribe = Scribe(RESULTS / "trajectory.jsonl")
    specs, _ = resolve_inputs(scribe, with_tier_a=False, tier_a_scaffold="scaffold_5aei")
    by_id = {s.case_id: s for s in specs}

    reps = {cid: cached_replicates(cid) for cid in CASES}
    have = len(reps["B1"])
    new_seeds = list(range(have, args.to))
    if not new_seeds:
        print(f"already at n={have}; nothing to buy")
        return 0

    runs = len(new_seeds) * len(CASES)
    print(f"reusing seeds 0..{have - 1} from cache; buying seeds {new_seeds} for {', '.join(CASES)}")
    print(f"{runs} model runs, ~{runs * 60 + 300} GPU-s, ~${(runs * 60 + 300) / 3600 * 4:.2f}\n", flush=True)
    scribe.append("economist", "buy G2 seed depth",
                  f"{runs} runs: seeds {new_seeds} on {', '.join(CASES)}, reusing {have} cached",
                  cost={"n_model_runs": runs, "gpu_seconds": runs * 60 + 300,
                        "usd": round((runs * 60 + 300) / 3600 * 4, 2)})

    for cid in CASES:
        for seed in new_seeds:
            value = score_one_seed(by_id[cid], seed)
            reps[cid].append(value)
            print(f"  {cid} seed {seed:>2}: {value:7.3f}", flush=True)

    print("\n== G2 as a function of how many seeds you buy ==\n")
    header = f"{'n':>3}  {'margin':>8}  {'range':>8} {'':>5}  {'CI/2':>8} {'':>5}  {'Welch t':>8}"
    print(header)
    print("-" * len(header))
    rows = []
    for n in range(3, args.to + 1):
        row = verdicts(reps["B1"][:n], reps["B4"][:n])
        rows.append(row)
        print(f"{n:>3}  {row['margin']:>8.2f}  {row['band_range']:>8.2f} {row['verdict_range']:>5}  "
              f"{row['band_ci_half']:>8.2f} {row['verdict_ci']:>5}  {row['welch_t']:>8.2f}")

    out = RESULTS / "g2_seed_depth.json"
    out.write_text(json.dumps({
        "purpose": "diagnostic only — does not amend the shipped G2 verdict",
        "shipped_verdict": {"band": "sample range at n=3", "value": 52.671, "status": "fail"},
        "replicates": {cid: [round(v, 4) for v in reps[cid]] for cid in CASES},
        "by_n": rows,
    }, indent=2) + "\n")
    final = rows[-1]
    scribe.append("referee", "G2 seed depth",
                  f"at n={final['n']}: margin {final['margin']}, range band {final['band_range']} "
                  f"({final['verdict_range']}), 95% CI half-width {final['band_ci_half']} ({final['verdict_ci']})",
                  outcome="ok")
    print(f"\nwrote {out}")
    print("the shipped predictions.json is unchanged — this is a diagnostic, not a re-scoring")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
