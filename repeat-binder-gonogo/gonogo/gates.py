"""Gate arithmetic. Each gate returns pass / fail / not_run, always with evidence.

A gate is never silently omitted. If it cannot be evaluated, that is a result
with a reason attached, and the reason names the thing that was missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import yaml

from .caseset import ComplexSpec
from .scorers import Score

DATA = Path(__file__).resolve().parent / "data"


@dataclass
class Verdict:
    gate_id: str
    status: str  # "pass" | "fail" | "not_run"
    evidence: str
    uncertainty_band: float | None = None
    detail: dict[str, Any] = field(default_factory=dict)

    def as_record(self) -> dict[str, Any]:
        record = {"gate_id": self.gate_id, "status": self.status, "evidence": self.evidence}
        if self.uncertainty_band is not None:
            record["uncertainty_band"] = round(self.uncertainty_band, 4)
        if self.detail:
            record["detail"] = self.detail
        return record


def _by_id(scores: Sequence[Score]) -> dict[str, Score]:
    return {score.case_id: score for score in scores}


# Student t, two-sided 95%, by rounded degrees of freedom. Enough resolution for a
# gate verdict; the band moves by less than a point across neighbouring df up here.
_T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26,
        10: 2.23, 11: 2.20, 12: 2.18, 14: 2.14, 16: 2.12, 18: 2.10, 20: 2.09, 25: 2.06,
        30: 2.04, 40: 2.02, 60: 2.00}


def _t_crit(df: float) -> float:
    return _T95[min(_T95, key=lambda k: abs(k - df))]


def _band_spec() -> dict[str, Any]:
    """The pre-registered G2 uncertainty band, read from cases.yaml.

    The estimator is part of the gate, not an implementation detail, so it lives
    in the case set beside the pass condition rather than in this file. An older
    case set without the key gets the original behaviour.
    """
    from .caseset import TASK_DIR

    cases = yaml.safe_load((TASK_DIR / "cases.yaml").read_text())
    g2 = next((g for g in cases.get("gates", []) if g.get("id") == "G2"), {})
    return g2.get("uncertainty_band") or {}


def _welch(a: Sequence[float], b: Sequence[float]) -> tuple[float, float]:
    """Standard error of the difference of means, and Welch-Satterthwaite df."""
    import statistics as st

    sea = st.stdev(a) / math.sqrt(len(a))
    seb = st.stdev(b) / math.sqrt(len(b))
    se = math.hypot(sea, seb)
    df = (sea**2 + seb**2) ** 2 / (sea**4 / (len(a) - 1) + seb**4 / (len(b) - 1))
    return se, df


# --------------------------------------------------------------------------


def g5_decoy_rejection(scores: Sequence[Score], positives: Sequence[str]) -> Verdict:
    table = _by_id(scores)
    decoys = {cid: s for cid, s in table.items() if cid.startswith("decoy_")}
    scored_positives = {cid: table[cid] for cid in positives if cid in table}

    if not decoys:
        return Verdict("G5", "not_run", "no decoys were scored")
    if not scored_positives:
        return Verdict("G5", "not_run", "no true positives were scored, so there is nothing to compare decoys against")

    weakest_id = min(scored_positives, key=lambda cid: scored_positives[cid].score)
    weakest = scored_positives[weakest_id].score
    offenders = {cid: round(s.score, 3) for cid, s in decoys.items() if s.score >= weakest}
    status = "fail" if offenders else "pass"
    evidence = (
        f"{len(decoys)} decoys scored against the weakest true positive "
        f"({weakest_id} = {weakest:.3f}); "
        + (f"{len(offenders)} scored at or above it: {offenders}" if offenders else "all scored below it")
    )
    return Verdict("G5", status, evidence, detail={"weakest_positive": weakest_id, "offenders": offenders})


def g2_specificity(scores: Sequence[Score]) -> Verdict:
    table = _by_id(scores)
    missing = [cid for cid in ("B1", "B4") if cid not in table]
    if missing:
        return Verdict("G2", "not_run", f"cases {missing} were not scored, so the specificity gate has no verdict")

    binder, non_binder = table["B1"], table["B4"]
    margin = binder.score - non_binder.score
    spec = _band_spec()

    # The retired estimator: the wider of the two per-case seed ranges. Kept so the
    # 2026-08-15 submission stays reproducible and so both numbers can be reported.
    legacy = max(
        b for b in (binder.band, non_binder.band) if not math.isnan(b)
    ) if not (math.isnan(binder.band) and math.isnan(non_binder.band)) else float("nan")

    if not spec:
        band, how = legacy, "widest seed spread of the two cases"
    else:
        n = min(len(binder.replicates), len(non_binder.replicates))
        floor = int(spec.get("min_seeds_per_case", 5))
        if n < floor:
            return Verdict(
                "G2", "not_run",
                f"{n} seeds per case is below the pre-registered minimum of {floor}: at this sample size "
                f"the t critical value dominates and the band reports the thinness of the sample rather "
                f"than the model's uncertainty. cases.yaml requires not_run here rather than a fail.",
                detail={"margin": round(margin, 4), "seeds_per_case": n, "min_seeds_per_case": floor},
            )
        se, df = _welch(binder.replicates, non_binder.replicates)
        band = _t_crit(df) * se
        how = f"95% CI half-width on the difference of means, Welch df={df:.1f}, pre-registered {spec.get('effective_from')}"

    if math.isnan(band):
        return Verdict(
            "G2",
            "not_run",
            "scored, but the scorer produced no uncertainty band (single seed); "
            "G2's verdict is a margin against a band and cannot be issued without one",
            detail={"score_B1": binder.score, "score_B4": non_binder.score, "margin": margin},
        )

    status = "pass" if margin > band else "fail"
    evidence = (
        f"score(B1, mouse cathepsin B, K_D 65.7 nM) = {binder.score:.3f}, "
        f"score(B4, human cathepsin B, no binding at 1000x) = {non_binder.score:.3f}, "
        f"margin {margin:.3f} {'>' if margin > band else '<='} band {band:.3f} ({how})"
    )
    detail = {
        "margin": round(margin, 4),
        "replicates_B1": [round(v, 3) for v in binder.replicates],
        "replicates_B4": [round(v, 3) for v in non_binder.replicates],
    }
    if spec and not math.isnan(legacy):
        # Both estimators, always, so a reader can see which one the verdict turned on.
        detail["band_retired_sample_range"] = round(legacy, 4)
        detail["verdict_under_retired_band"] = "pass" if margin > legacy else "fail"
    return Verdict("G2", status, evidence, uncertainty_band=band, detail=detail)


def g1_interface_recovery(scores: Sequence[Score], structures_available: bool) -> Verdict:
    """Interface recovery, as specified in DockQ.

    DockQ is not in the proto-tools catalogue and is not deployable through it,
    so the metric the gate is written in cannot be produced by this stack alone.
    That is reported rather than silently swapped for a different metric.
    """
    if not structures_available:
        return Verdict(
            "G1",
            "not_run",
            "no predicted structures were produced (the structure-prediction tool is not deployed), "
            "so no interface could be compared to a deposit",
        )
    return Verdict(
        "G1",
        "not_run",
        "structures exist but DockQ is not available: it is absent from the proto-tools catalogue "
        "(`search_tools('dockq')` returns no hit) and the gate's pass condition is written in DockQ >= 0.49. "
        "Closing it needs the standalone DockQ package, or the gate rewritten around a metric this stack has "
        "(USalign multimer TM-score, or interface RMSD computed locally)",
    )


def _spearman(xs: Sequence[float], ys: Sequence[float]) -> float:
    def rank(values: Sequence[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            shared = (i + j) / 2 + 1
            for k in range(i, j + 1):
                ranks[order[k]] = shared
            i = j + 1
        return ranks

    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def g3_affinity_ranking(scores: Sequence[Score], specs: Sequence[ComplexSpec]) -> Verdict:
    """Spearman rho against log K_D, per pocket, with the published bar attached."""
    truth = {spec.case_id: spec.truth for spec in specs if spec.tier == "A"}
    table = {cid: s for cid, s in _by_id(scores).items() if cid in truth}
    if not table:
        return Verdict(
            "G3",
            "not_run",
            "no Tier A variants were scored. The per-variant K_D table itself is available — it was extracted "
            "from Table S1 of PMC12202725 into gonogo/data/tier_a_kd.yaml (63 variants across 5 pockets) — so the "
            "blocker is model runs, not ground truth",
        )

    baselines = yaml.safe_load((DATA / "baselines.yaml").read_text())
    per_pocket: dict[str, Any] = {}
    for pocket in {truth[cid]["pocket"] for cid in table}:
        ids = [cid for cid in table if truth[cid]["pocket"] == pocket]
        if len(ids) < 4:
            per_pocket[pocket] = {"rho": None, "n": len(ids), "note": "too few variants scored to correlate"}
            continue
        predicted = [table[cid].score for cid in ids]
        # higher score should mean tighter binding, i.e. lower K_D
        observed = [-math.log10(truth[cid]["kd_nM"]) for cid in ids]
        rho = _spearman(predicted, observed)
        published = baselines["per_pocket"].get(pocket, {})
        per_pocket[pocket] = {"rho": round(rho, 3), "n": len(ids), "published_baselines": published}

    usable = [v["rho"] for v in per_pocket.values() if v.get("rho") is not None]
    if not usable:
        return Verdict("G3", "not_run", "no pocket had enough scored variants to correlate", detail=per_pocket)
    status = "pass" if all(rho >= 0.5 for rho in usable) else "fail"
    evidence = "; ".join(
        f"{pocket}: rho={value['rho']} over n={value['n']}" for pocket, value in per_pocket.items() if value.get("rho")
    )
    evidence += (
        ". Comparison to the published physics baselines is reported per pocket in `detail` and is not "
        "collapsed into a single number: the paper reports Pearson R for most pockets and Spearman rho for a few, "
        "and G3 is specified in rho"
    )
    return Verdict("G3", status, evidence, detail=per_pocket)


def _auc(positive: Sequence[float], negative: Sequence[float]) -> float:
    """Mann-Whitney AUC: P(confidence of a correct call > that of an incorrect one)."""
    if not positive or not negative:
        return float("nan")
    wins = sum((1.0 if p > n else 0.5 if p == n else 0.0) for p in positive for n in negative)
    return wins / (len(positive) * len(negative))


def g4_calibration(scores: Sequence[Score], specs: Sequence[ComplexSpec], decision_threshold: float | None) -> Verdict:
    """Is the confidence informative about whether the score was right?

    Correctness needs a decision threshold, and choosing one after seeing the
    answers would be cheating. The threshold used is the highest decoy score
    from G5 — fixed by the cheapest gate, before any of this is looked at.
    """
    if decision_threshold is None:
        return Verdict(
            "G4",
            "not_run",
            "no decision threshold: G5 did not run, and the threshold that separates 'called a binder' from "
            "'called a non-binder' is defined as the highest decoy score so that it cannot be tuned after the fact",
        )
    expectation = {spec.case_id: spec.expectation for spec in specs}
    correct_confidences, wrong_confidences = [], []
    for score in scores:
        label = expectation.get(score.case_id)
        if label not in {"binder", "non_binder", "decoy"}:
            continue
        called_binder = score.score > decision_threshold
        is_binder = label == "binder"
        (correct_confidences if called_binder == is_binder else wrong_confidences).append(score.confidence)

    if not wrong_confidences:
        return Verdict(
            "G4",
            "not_run",
            f"every one of {len(correct_confidences)} labelled predictions landed on the correct side of the "
            f"threshold ({decision_threshold:.3f}); AUC of confidence against correctness is undefined with one class. "
            "Reported as not_run rather than as a perfect score",
            detail={"n_correct": len(correct_confidences), "n_wrong": 0},
        )
    auc = _auc(correct_confidences, wrong_confidences)
    status = "pass" if auc >= 0.7 else "fail"
    return Verdict(
        "G4",
        status,
        f"AUC={auc:.3f} of confidence separating {len(correct_confidences)} correct from "
        f"{len(wrong_confidences)} incorrect calls, at decision threshold {decision_threshold:.3f} "
        f"(the highest decoy score, fixed by G5)",
        detail={"auc": round(auc, 4), "n_correct": len(correct_confidences), "n_wrong": len(wrong_confidences)},
    )
