"""Pytest-based verifier. Run by BenchFlow after the agent completes.

Checks the shape and honesty of results/predictions.json, then re-derives the
two gates that can be evaluated from the submitted numbers alone (G2, G5).
Gates that need external ground truth (G1, G3, G4) are checked for presence and
for a stated reason when not run — not for their scientific verdict.
"""

import json
import os
from pathlib import Path

import pytest
import yaml

TASK_DIR = Path(__file__).resolve().parents[1]
# Overridable so the verifier can be exercised outside the task container.
PREDICTIONS = Path(os.environ.get("PREDICTIONS_PATH", "/app/results/predictions.json"))

REQUIRED_PREDICTION_FIELDS = {"case_id", "score", "confidence", "confidence_metric", "tool", "ran_on"}
VALID_GATE_STATUS = {"pass", "fail", "not_run"}


@pytest.fixture(scope="module")
def cases():
    return yaml.safe_load((TASK_DIR / "cases.yaml").read_text())


@pytest.fixture(scope="module")
def submission():
    if not PREDICTIONS.exists():
        pytest.fail(f"no submission at {PREDICTIONS}")
    return json.loads(PREDICTIONS.read_text())


@pytest.fixture(scope="module")
def by_case(submission):
    return {p["case_id"]: p for p in submission.get("predictions", [])}


@pytest.fixture(scope="module")
def by_gate(submission):
    return {g["gate_id"]: g for g in submission.get("gates", [])}


def test_every_prediction_is_fully_attributed(submission):
    """No bare point estimates: each score carries confidence and provenance."""
    predictions = submission.get("predictions", [])
    assert predictions, "submission contains no predictions"
    for pred in predictions:
        missing = REQUIRED_PREDICTION_FIELDS - pred.keys()
        assert not missing, f"{pred.get('case_id', '?')} is missing {sorted(missing)}"
        assert isinstance(pred["confidence"], (int, float)), (
            f"{pred['case_id']} confidence is not numeric"
        )


def test_every_gate_is_accounted_for(cases, by_gate):
    """A gate that was not run must say so, with a reason. Silence fails."""
    for gate in cases["gates"]:
        gid = gate["id"]
        assert gid in by_gate, f"gate {gid} is absent from the submission"
        status = by_gate[gid].get("status")
        assert status in VALID_GATE_STATUS, f"gate {gid} has invalid status {status!r}"
        assert by_gate[gid].get("evidence"), f"gate {gid} reports {status} with no evidence"


def test_g5_decoys_score_below_every_true_positive(cases, by_case, by_gate):
    """Cheap gate, run first: shuffles and non-cognate pairs must lose."""
    if by_gate.get("G5", {}).get("status") == "not_run":
        pytest.skip("G5 reported not_run with a reason")

    decoys = {cid: p for cid, p in by_case.items() if cid.startswith("decoy_")}
    assert decoys, "G5 claims to have run but no decoy_* predictions were submitted"

    positive_ids = [p["id"] for p in cases["tier_b_specificity"]["positives"]]
    positives = {cid: by_case[cid] for cid in positive_ids if cid in by_case}
    assert positives, "no scored positives to compare decoys against"

    weakest = min(p["score"] for p in positives.values())
    offenders = {cid: p["score"] for cid, p in decoys.items() if p["score"] >= weakest}
    assert not offenders, f"decoys scored at or above the weakest true positive ({weakest}): {offenders}"


def test_g2_separates_mouse_from_human_cathepsin_b(by_case, by_gate):
    """The hard gate. Same DARPin, 83%-identical targets, opposite outcomes."""
    if by_gate.get("G2", {}).get("status") == "not_run":
        pytest.fail("G2 cannot be skipped — it is the gate this validation exists for")

    for cid in ("B1", "B4"):
        assert cid in by_case, f"case {cid} was not scored"

    binder, non_binder = by_case["B1"], by_case["B4"]
    margin = binder["score"] - non_binder["score"]
    assert margin > 0, (
        f"scored the non-binder (B4, human cathepsin B) at least as well as the "
        f"65.7 nM binder (B1, mouse): {binder['score']} vs {non_binder['score']}"
    )

    band = by_gate["G2"].get("uncertainty_band")
    assert band is not None, "G2 must state the uncertainty band the margin is judged against"
    assert margin > band, f"margin {margin} does not exceed the model's own uncertainty band {band}"


def test_reported_gate_verdicts_match_the_numbers(by_case, by_gate):
    """Guard against a pass being asserted in prose while the scores disagree."""
    g2 = by_gate.get("G2", {})
    if g2.get("status") == "pass" and {"B1", "B4"} <= by_case.keys():
        assert by_case["B1"]["score"] > by_case["B4"]["score"], (
            "G2 is reported as pass but the submitted scores say otherwise"
        )
