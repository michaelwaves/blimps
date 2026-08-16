from __future__ import annotations

from pathlib import Path

import gvpa_full_pipeline as pipeline


def candidate(arm: str, backbone: str, perplexity: float, clashes: int = 0) -> dict:
    return {
        "arm": arm,
        "backbone_id": backbone,
        "candidate_id": f"{backbone}_{perplexity}",
        "mpnn_perplexity": perplexity,
        "backbone_clashes_lt2A": clashes,
    }


def test_selection_uses_distinct_clash_free_backbones() -> None:
    rows = [
        candidate("constrained", "c00", 1.0, clashes=4),
        candidate("constrained", "c01", 1.2),
        candidate("constrained", "c01", 1.1),
        candidate("constrained", "c02", 1.3),
        candidate("free", "f00", 1.0),
        candidate("free", "f00", 0.9),
        candidate("free", "f01", 1.1),
    ]

    selected = pipeline._selected_for_boltz(rows, top_per_arm=2)

    assert [row["backbone_id"] for row in selected] == ["c01", "c02", "f00", "f01"]
    assert all(row["backbone_clashes_lt2A"] == 0 for row in selected)


def test_selection_fails_when_clash_free_budget_is_unavailable() -> None:
    rows = [
        candidate("constrained", "c00", 1.0, clashes=1),
        candidate("free", "f00", 1.0),
    ]

    try:
        pipeline._selected_for_boltz(rows, top_per_arm=1)
    except ValueError as exc:
        assert "constrained" in str(exc)
        assert "clash-free" in str(exc)
    else:
        raise AssertionError("selection accepted a clashing backbone")


def test_committed_backbones_are_complete_and_portable() -> None:
    rows = pipeline.backbone_rows()

    assert len(rows) == 16
    assert {row["arm"] for row in rows} == {"constrained", "free"}
    assert all(not Path(row["path"]).is_absolute() for row in rows)
    assert all(pipeline._repo_absolute(row["path"]).is_file() for row in rows)
