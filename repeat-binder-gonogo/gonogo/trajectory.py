"""The Scribe: an append-only record of what the run did, as it did it.

One line per step, written the moment the step happens. Dead ends stay in.
A run that is reconstructed afterwards loses the only thing that shows the
decisions responding to evidence.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROLES = ("librarian", "architect", "referee", "economist", "scribe")


@dataclass
class Scribe:
    path: Path
    echo: bool = True
    _t0: float = field(default_factory=time.time)

    def __post_init__(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def append(
        self,
        role: str,
        step: str,
        rationale: str,
        *,
        outcome: str = "ok",
        cost: dict[str, Any] | None = None,
        evidence: Any = None,
        **extra: Any,
    ) -> dict[str, Any]:
        if role not in ROLES:
            raise ValueError(f"unknown role {role!r}; expected one of {ROLES}")
        record = {
            "t": round(time.time() - self._t0, 3),
            "role": role,
            "step": step,
            "rationale": rationale,
            "outcome": outcome,
            "cost": cost or {"gpu_seconds": 0, "usd": 0.0},
            "evidence": evidence,
            **extra,
        }
        with self.path.open("a") as fh:
            fh.write(json.dumps(record, default=str) + "\n")
        if self.echo:
            marker = {"ok": "·", "blocked": "✗", "gap": "?", "veto": "⛔"}.get(outcome, "·")
            print(f"  {marker} [{role}] {step}: {rationale}", flush=True)
        return record

    def read(self) -> list[dict[str, Any]]:
        if not self.path.exists():
            return []
        return [json.loads(line) for line in self.path.read_text().splitlines() if line.strip()]


def default_scribe(run_dir: Path) -> Scribe:
    return Scribe(run_dir / "trajectory.jsonl", echo=os.environ.get("GONOGO_QUIET") != "1")
