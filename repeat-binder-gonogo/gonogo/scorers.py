"""Scorers: the only place a number enters the pipeline.

Contract, enforced by construction rather than by review:

* Every score is higher-is-better, in one convention across all cases.
* Every score arrives with a confidence value produced by the same model, the
  name of that confidence metric, the tool key, and where it ran.
* Every scorer states its own uncertainty band. A scorer that cannot produce
  one cannot be used for G2, because G2's verdict is a margin measured against
  that band.
* A scorer that cannot run raises `Unavailable` with a reason. Nothing
  downstream is allowed to substitute a guess.
"""

from __future__ import annotations

import hashlib
import statistics
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

from .caseset import ComplexSpec


class Unavailable(RuntimeError):
    """Raised when a scorer cannot run. Carries the reason into the submission."""


def _last_deploy_error(app: str) -> str:
    """Why the last deploy of `app` failed, if one was attempted here.

    "not deployed" is a symptom. The submission should carry the cause, which
    `proto-tools deploy` leaves in its log and nowhere else.
    """
    import re
    from pathlib import Path

    log = Path(__file__).resolve().parent.parent / "logs" / f"deploy.{app}.log"
    if not log.exists():
        return ""
    text = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", log.read_text(errors="replace"))
    errors = [
        line.strip(" │|")
        for line in text.splitlines()
        if ("Error" in line or "error" in line or "Please " in line) and len(line.strip(" │|")) > 12
    ]
    return f". Last deploy attempt failed with: {errors[-1]}" if errors else ""


@dataclass
class Score:
    case_id: str
    score: float
    confidence: float
    confidence_metric: str
    tool: str
    ran_on: str
    band: float
    replicates: list[float] = field(default_factory=list)
    extras: dict[str, Any] = field(default_factory=dict)

    def as_prediction(self) -> dict[str, Any]:
        record = {
            "case_id": self.case_id,
            "score": round(self.score, 4),
            "confidence": round(self.confidence, 4),
            "confidence_metric": self.confidence_metric,
            "tool": self.tool,
            "ran_on": self.ran_on,
            "uncertainty": {
                "band": round(self.band, 4),
                "source": "spread across independent seeds",
                "replicates": [round(v, 4) for v in self.replicates],
            },
        }
        record.update(self.extras)
        return record


class Scorer:
    key = "abstract"
    ran_on = "unknown"
    confidence_metric = "none"
    score_definition = "undefined"

    def preflight(self) -> None:
        """Raise Unavailable if this scorer cannot run right now."""

    def cost_estimate(self, n: int) -> dict[str, Any]:
        return {"n_predictions": n, "gpu_seconds": None, "usd": None}

    def score(self, specs: Sequence[ComplexSpec]) -> Iterable[Score]:
        raise NotImplementedError


# --------------------------------------------------------------------------


class Boltz2Scorer(Scorer):
    """Co-folding score from Boltz-2, run on the user's Modal deployment.

    score      100 x the binder/target pairwise interface ipTM (higher is better,
               already in the required convention, so no sign flip and no note).
    confidence ipTM for the complex, the model's own interface confidence.
    band       the spread of the score across independent seeds. This is the
               only honest uncertainty band available from a single model: it
               measures run-to-run instability, not model error, and G2's margin
               is judged against it.
    """

    key = "boltz2-prediction"
    ran_on = "modal"
    confidence_metric = "iptm"
    score_definition = "100 * pair_chains_iptm[binder][target]"

    def __init__(self, seeds: Sequence[int] = (0, 1, 2), recycling_steps: int = 3, use_msa: bool = True) -> None:
        self.seeds = list(seeds)
        self.recycling_steps = recycling_steps
        self.use_msa = use_msa

    def preflight(self) -> None:
        from proto_tools.mcp import tools as mcp_tools

        info = mcp_tools.workspace_info()
        if not info.get("authenticated"):
            raise Unavailable("no Modal credentials: `modal token new` has not been run for this workspace")
        if not info.get("environment_exists"):
            raise Unavailable(
                f"Modal environment {info.get('environment')!r} does not exist in workspace "
                f"{info.get('workspace')!r}; create it with "
                f"`proto-tools deploy --create-env --env {info.get('environment')}`"
            )
        catalogue = mcp_tools.list_tools(deployed_only=True)
        if isinstance(catalogue, dict):  # the MCP surface wraps the list, the function does not
            catalogue = catalogue.get("result", [])
        deployed = {tool["tool"] for tool in catalogue}
        if self.key not in deployed:
            raise Unavailable(
                f"{self.key} is in the catalogue but not deployed to this workspace; "
                f"deploy it before any GPU gate can run" + _last_deploy_error("boltz2")
            )

    def cost_estimate(self, n: int) -> dict[str, Any]:
        runs = n * len(self.seeds)
        # ~1 GPU-minute per medium complex per seed once the container is warm,
        # plus a one-off container start and weight load.
        gpu_seconds = runs * 60 + 300
        return {
            "n_predictions": n,
            "n_model_runs": runs,
            "gpu_seconds": gpu_seconds,
            "usd": round(gpu_seconds / 3600 * 4.0, 2),
            "assumption": "~60 GPU-s per complex per seed on an A100 at ~$4/GPU-hour, plus one 5-minute cold start",
        }

    def score(self, specs: Sequence[ComplexSpec]) -> Iterable[Score]:
        from proto_tools.mcp.tools import run_tool

        for spec in specs:
            values, confidences, paes = [], [], []
            for seed in self.seeds:
                result = run_tool(
                    self.key,
                    inputs={
                        "complexes": [
                            {
                                "chains": [
                                    {"id": "B", "sequence": spec.binder.sequence, "entity_type": "protein"},
                                    {"id": "T", "sequence": spec.target.sequence, "entity_type": "protein"},
                                ]
                            }
                        ]
                    },
                    config={
                        "seed": seed,
                        "recycling_steps": self.recycling_steps,
                        "use_msa": self.use_msa,
                        "diffusion_samples": 1,
                    },
                    device="modal",
                )
                if not result.get("ok", True):
                    raise Unavailable(f"{spec.case_id}: {self.key} failed — {result.get('error')}")
                metrics = result["result"]["structures"][0]["metrics"]
                pair = metrics["pair_chains_iptm"]
                values.append(100.0 * float(pair[0][1]))
                confidences.append(float(metrics["iptm"]))
                paes.append(float(metrics["avg_pae"]))

            band = (max(values) - min(values)) if len(values) > 1 else float("nan")
            yield Score(
                case_id=spec.case_id,
                score=statistics.mean(values),
                confidence=statistics.mean(confidences),
                confidence_metric=self.confidence_metric,
                tool=self.key,
                ran_on=self.ran_on,
                band=band,
                replicates=values,
                extras={"avg_pae": round(statistics.mean(paes), 3), "seeds": self.seeds},
            )


# --------------------------------------------------------------------------


class SyntheticScorer(Scorer):
    """NOT SCIENCE. A deterministic stand-in used to test the harness itself.

    It exists so the gate arithmetic, the submission schema, and the verifier
    wiring can be exercised without a GPU. Its numbers are a hash of the input
    sequences and mean nothing. Every prediction it produces is stamped
    `synthetic: true`, and the CLI refuses to write them into the bench
    submission path.
    """

    key = "synthetic-demo-scorer"
    ran_on = "local(synthetic)"
    confidence_metric = "synthetic_pseudo_iptm"
    score_definition = "deterministic hash of the input sequences — carries no physical meaning"

    # Five, not three: a seed costs nothing here (the value is a hash, not a fold), and
    # G2's pre-registered band declines to issue a verdict below five seeds per case.
    # A self-test that stopped short of that would exercise the abstention, not the gate.
    def __init__(self, seeds: Sequence[int] = (0, 1, 2, 3, 4), separate_known_outcomes: bool = True) -> None:
        self.seeds = list(seeds)
        self.separate_known_outcomes = separate_known_outcomes

    def _base(self, spec: ComplexSpec, seed: int) -> float:
        digest = hashlib.sha256(f"{spec.binder.sequence}|{spec.target.sequence}|{seed}".encode()).hexdigest()
        noise = int(digest[:8], 16) / 0xFFFFFFFF  # 0..1
        if not self.separate_known_outcomes:
            return 40.0 + 20.0 * noise
        anchor = {"binder": 70.0, "non_binder": 45.0, "decoy": 25.0}.get(spec.expectation, 50.0)
        return anchor + 6.0 * (noise - 0.5)

    def score(self, specs: Sequence[ComplexSpec]) -> Iterable[Score]:
        for spec in specs:
            values = [self._base(spec, seed) for seed in self.seeds]
            yield Score(
                case_id=spec.case_id,
                score=statistics.mean(values),
                confidence=min(0.99, statistics.mean(values) / 100.0),
                confidence_metric=self.confidence_metric,
                tool=self.key,
                ran_on=self.ran_on,
                band=max(values) - min(values),
                replicates=values,
                extras={"synthetic": True},
            )


# --------------------------------------------------------------------------


class CachedScorer(Scorer):
    """Persist each complex's score the moment it is produced.

    A gate sweep is dozens of independent GPU calls, and the submission is only
    written at the end. Without this, a run killed at 90% has bought nothing:
    the next attempt re-runs every complex from scratch. With it, a restart pays
    only for what is missing.

    The cache key covers the tool, the score definition, the seeds and both
    sequences, so changing any of them is a miss rather than a stale hit.
    """

    def __init__(self, inner: Scorer, cache_dir) -> None:
        from pathlib import Path

        self.inner = inner
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.hits = 0
        self.misses = 0

    # the wrapper must look exactly like what it wraps
    key = property(lambda self: self.inner.key)
    ran_on = property(lambda self: self.inner.ran_on)
    confidence_metric = property(lambda self: self.inner.confidence_metric)
    score_definition = property(lambda self: self.inner.score_definition)

    def preflight(self) -> None:
        self.inner.preflight()

    def cost_estimate(self, n: int) -> dict[str, Any]:
        return self.inner.cost_estimate(n)

    def _path(self, spec: ComplexSpec):
        material = "|".join(
            [
                self.inner.key,
                self.inner.score_definition,
                str(getattr(self.inner, "seeds", "")),
                spec.case_id,
                spec.binder.sequence,
                spec.target.sequence,
            ]
        )
        digest = hashlib.sha256(material.encode()).hexdigest()[:16]
        return self.cache_dir / f"{spec.case_id}.{digest}.json"

    def score(self, specs: Sequence[ComplexSpec]) -> Iterable[Score]:
        import json

        for spec in specs:
            path = self._path(spec)
            if path.exists():
                self.hits += 1
                yield Score(**json.loads(path.read_text()))
                continue
            for score in self.inner.score([spec]):
                self.misses += 1
                path.write_text(json.dumps(score.__dict__, default=str))
                yield score


SCORERS = {"boltz2": Boltz2Scorer, "synthetic": SyntheticScorer}
