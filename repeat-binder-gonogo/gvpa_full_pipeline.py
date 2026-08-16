"""Finish the lattice-aware GvpA binder pipeline.

Stages:
  mpnn     inverse-fold every RFdiffusion3 backbone (binder chain A only)
  compile  join backbone geometry and ProteinMPNN sequence metrics
  boltz    screen two candidates/arm, then replicate the winner/arm

The constrained and free arms are kept separate throughout.  Boltz-2 scores are
reported as triage evidence, not proof of epitope recovery: the control run's G1
gate failed on repeat proteins and therefore does not license that claim.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure


ROOT = Path(__file__).resolve().parent
WORK = ROOT / "designs" / "gvpa"
ARMS = ("constrained", "free")
TARGET_SEQUENCE = "SIQKSTNSSSLAEVIDRILDKGIVIDAFARVSVVGIEILTIEARVVIASVDTWLRYAEAVGLLRD"


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=str) + "\n")


def _repo_relative(path: Path) -> str:
    """Store paths relative to the repository so outputs survive a clone."""
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def _repo_absolute(path: str | Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def _saved_path(descriptor: dict, local_dir: Path) -> Path:
    """Resolve proto-tools descriptors copied from a different workstation."""
    saved = Path(descriptor["_saved_to"])
    if saved.exists():
        return saved
    local = local_dir / saved.name
    if local.exists():
        return local
    raise FileNotFoundError(f"neither {saved} nor portable fallback {local} exists")


def backbone_rows() -> list[dict]:
    rows: list[dict] = []
    cif_dir = WORK / "full_pipeline" / "rfd3_inputs"
    cif_dir.mkdir(parents=True, exist_ok=True)

    # The 16 materialized backbones are the portable, committed experiment
    # inputs. A fresh clone can start at ProteinMPNN without the raw Modal
    # response descriptors that created them.
    materialized = sorted(cif_dir.glob("[cf][0-9][0-9].cif"))
    if materialized:
        if len(materialized) != 16:
            raise ValueError(f"expected 16 materialized backbones, found {len(materialized)}")
        for path in materialized:
            arm = "constrained" if path.stem.startswith("c") else "free"
            rows.append({"arm": arm, "backbone": int(path.stem[1:]),
                         "backbone_id": path.stem, "path": _repo_relative(path),
                         "source_path": _repo_relative(path)})
        return rows

    for arm in ARMS:
        raw_dir = WORK / f"rfd3_{arm}"
        index = json.loads((raw_dir / "designed_structures_0_structures.json").read_text())
        for i, item in enumerate(index):
            path = _saved_path(item["structure"]["structure"], raw_dir)
            backbone_id = f"{arm[:1]}{i:02d}"
            cif_path = cif_dir / f"{backbone_id}.cif"
            shutil.copyfile(path, cif_path)
            rows.append({"arm": arm, "backbone": i, "backbone_id": backbone_id,
                         "path": _repo_relative(cif_path), "source_path": _repo_relative(cif_path)})
    return rows


def run_mpnn(n_per_backbone: int) -> None:
    from proto_tools.mcp.tools import run_tool

    rows = backbone_rows()
    out = WORK / "full_pipeline" / "mpnn"
    out.mkdir(parents=True, exist_ok=True)
    print(f"ProteinMPNN: {len(rows)} backbones x {n_per_backbone} sequences; redesign binder chain A only", flush=True)
    result = run_tool(
        "proteinmpnn-sample",
        inputs={"inputs": [
            {"structure": str(_repo_absolute(row["path"])), "chains_to_redesign": {"chains": ["A"]}}
            for row in rows
        ]},
        config={
            "num_sequences_per_structure": n_per_backbone,
            "batch_size": n_per_backbone,
            "temperature": 0.15,
            "excluded_amino_acids": ["C"],
            "seed": 0,
        },
        device="modal",
        output_dir=str(out),
    )
    if not result.get("ok", True):
        raise RuntimeError(result.get("error"))
    _write_json(out / "raw.json", result["result"])
    _write_json(out / "backbones.json", rows)
    print(f"wrote {out / 'raw.json'}")


def _load_saved_list(descriptor: object, local_dir: Path) -> list[dict]:
    if isinstance(descriptor, list):
        return descriptor
    if not isinstance(descriptor, dict) or "_saved_to" not in descriptor:
        raise ValueError(f"expected saved-list descriptor, got {descriptor!r}")
    return json.loads(_saved_path(descriptor, local_dir).read_text())


def _site_metrics(cif_path: str, hotspots: set[str]) -> dict:
    atoms = get_structure(CIFFile.read(_repo_absolute(cif_path)), model=1)
    heavy = atoms.element != "H"
    binder = atoms[(atoms.chain_id == "A") & heavy]
    target = atoms[(atoms.chain_id != "A") & heavy]

    # RFD3 inserts binder as A and shifts original target A:E to output B:F.
    output_hotspots = {(chr(ord(chain) + 1), int(resid)) for chain, resid in
                       ((x[0], x[1:]) for x in hotspots)}
    hot_mask = np.array([(c, int(r)) in output_hotspots for c, r in zip(target.chain_id, target.res_id)])
    hotspot_atoms = target[hot_mask]

    # The arrays are small (~500 x ~2300). Broadcasting avoids scipy, which is
    # intentionally unavailable in this DockQ-compatible environment.
    dist = np.linalg.norm(binder.coord[:, None, :] - target.coord[None, :, :], axis=2)
    contacted_atom = np.min(dist, axis=0) < 5.0
    contacted_output = {(str(c), int(r)) for c, r in zip(target.chain_id[contacted_atom], target.res_id[contacted_atom])}
    contacted_original = {f"{chr(ord(c) - 1)}{r}" for c, r in contacted_output}
    hotspot_hits = contacted_original & hotspots

    if len(hotspot_atoms):
        min_hotspot_distance = float(np.min(np.linalg.norm(
            binder.coord[:, None, :] - hotspot_atoms.coord[None, :, :], axis=2
        )))
    else:
        min_hotspot_distance = float("nan")
    backbone_names = {"N", "CA", "C", "O"}
    binder_bb = binder[np.isin(binder.atom_name, list(backbone_names))]
    target_bb = target[np.isin(target.atom_name, list(backbone_names))]
    bb_dist = np.linalg.norm(binder_bb.coord[:, None, :] - target_bb.coord[None, :, :], axis=2)
    return {
        "contact_residues": sorted(contacted_original),
        "n_contact_residues": len(contacted_original),
        "target_chains_contacted": len({x[0] for x in contacted_original}),
        "hotspot_hits": sorted(hotspot_hits),
        "n_hotspot_hits": len(hotspot_hits),
        "contact_precision": round(len(hotspot_hits) / max(1, len(contacted_original)), 4),
        "hotspot_recall": round(len(hotspot_hits) / max(1, len(hotspots)), 4),
        "min_hotspot_distance_A": round(min_hotspot_distance, 3),
        "min_target_distance_A": round(float(np.min(dist)), 3),
        "min_backbone_distance_A": round(float(np.min(bb_dist)), 3),
        "backbone_clashes_lt2A": int(np.sum(bb_dist < 2.0)),
    }


def compile_candidates() -> list[dict]:
    out = WORK / "full_pipeline"
    mpnn = json.loads((out / "mpnn" / "raw.json").read_text())
    # Re-resolve from the committed CIFs instead of trusting path strings in a
    # copied proto-tools descriptor, which may name the machine that ran it.
    backbones = backbone_rows()
    hotspots = set(json.loads((WORK / "site.json").read_text())["gvpc_footprint_residues"])
    mpnn_dir = out / "mpnn"
    design_sets = _load_saved_list(mpnn["design_sets"], mpnn_dir)
    if len(design_sets) != len(backbones):
        raise ValueError(f"{len(design_sets)} design sets for {len(backbones)} backbones")

    rows: list[dict] = []
    for source, design_set in zip(backbones, design_sets):
        site = _site_metrics(source["path"], hotspots)
        for seq_index, complex_ in enumerate(_load_saved_list(design_set["complexes"], mpnn_dir)):
            designed = [chain for chain, flag in zip(complex_["chains"], complex_["designed"]) if flag]
            if len(designed) != 1 or designed[0]["id"] != "A":
                raise ValueError(f"expected exactly binder chain A redesigned: {designed}")
            metrics = complex_["metrics"]
            rows.append({
                "candidate_id": f"{source['backbone_id']}_s{seq_index}",
                "arm": source["arm"],
                "backbone": source["backbone"],
                "backbone_id": source["backbone_id"],
                "backbone_path": source["path"],
                "sequence": designed[0]["sequence"],
                "length": len(designed[0]["sequence"]),
                "mpnn_perplexity": round(float(metrics["perplexity"]), 4),
                "mpnn_sequence_recovery": round(float(metrics["sequence_recovery"]), 4),
                **site,
            })

    rows.sort(key=lambda r: (r["arm"], r["mpnn_perplexity"], r["candidate_id"]))
    _write_json(out / "candidates.json", rows)
    with (out / "candidates.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k not in {"contact_residues", "hotspot_hits"}])
        writer.writeheader()
        writer.writerows([{k: v for k, v in row.items() if k in writer.fieldnames} for row in rows])
    print(f"compiled {len(rows)} candidates -> {out / 'candidates.json'}")
    return rows


def _selected_for_boltz(rows: list[dict], top_per_arm: int) -> list[dict]:
    selected: list[dict] = []
    for arm in ARMS:
        # First keep the best sequence on every distinct backbone, then rank
        # backbones by MPNN perplexity. This prevents four nearly identical
        # samples from one backbone consuming the validation budget.
        best_by_backbone: dict[str, dict] = {}
        for row in (r for r in rows if r["arm"] == arm):
            old = best_by_backbone.get(row["backbone_id"])
            if old is None or row["mpnn_perplexity"] < old["mpnn_perplexity"]:
                best_by_backbone[row["backbone_id"]] = row
        clash_free = [r for r in best_by_backbone.values() if r["backbone_clashes_lt2A"] == 0]
        if len(clash_free) < top_per_arm:
            raise ValueError(f"{arm}: requested {top_per_arm} backbones but only "
                             f"{len(clash_free)} are clash-free")
        selected.extend(sorted(clash_free, key=lambda r: r["mpnn_perplexity"])[:top_per_arm])
    return selected


def _complex(sequence: str) -> dict:
    return {"chains": [
        {"id": "A", "sequence": sequence, "entity_type": "protein"},
        *[{"id": chain, "sequence": TARGET_SEQUENCE, "entity_type": "protein"} for chain in "BCDEF"],
    ]}


def _run_boltz_batch(rows: list[dict], seed: int, label: str) -> list[dict]:
    from proto_tools.mcp.tools import run_tool

    out = WORK / "full_pipeline" / "boltz" / label
    out.mkdir(parents=True, exist_ok=True)
    raw_path = out / "raw.json"
    if raw_path.exists():
        print(f"Boltz-2 {label}: reusing completed raw result", flush=True)
        payload = json.loads(raw_path.read_text())
    else:
        print(f"Boltz-2 {label}: {len(rows)} complexes, seed={seed}, no MSA", flush=True)
        result = run_tool(
            "boltz2-prediction",
            inputs={"complexes": [_complex(row["sequence"]) for row in rows]},
            config={"seed": seed, "recycling_steps": 3, "use_msa": False,
                    "diffusion_samples": 1, "include_pae_matrix": False},
            device="modal",
            output_dir=str(out),
        )
        if not result.get("ok", True):
            raise RuntimeError(result.get("error"))
        payload = result["result"]
        _write_json(raw_path, payload)
    structures = _load_saved_list(payload["structures"], out)
    if len(structures) != len(rows):
        raise ValueError(f"Boltz returned {len(structures)} structures for {len(rows)} complexes")
    scored: list[dict] = []
    for row, structure in zip(rows, structures):
        metrics = structure["metrics"]
        pair = metrics["pair_chains_iptm"]
        scored.append({
            "candidate_id": row["candidate_id"],
            "arm": row["arm"],
            "seed": seed,
            "binder_target_iptm": round(max(float(pair[0][i]) for i in range(1, 6)), 4),
            "complex_iptm": round(float(metrics["iptm"]), 4),
            "confidence_score": round(float(metrics["confidence_score"]), 4),
            "avg_pae_A": round(float(metrics["avg_pae"]), 3),
        })
    _write_json(out / "scores.json", scored)
    return scored


def run_boltz(top_per_arm: int) -> dict:
    out = WORK / "full_pipeline"
    rows = json.loads((out / "candidates.json").read_text())
    selected = _selected_for_boltz(rows, top_per_arm)
    screen = _run_boltz_batch(selected, seed=0, label="screen_seed0")

    winners: list[dict] = []
    for arm in ARMS:
        best_score = max((s for s in screen if s["arm"] == arm), key=lambda s: s["binder_target_iptm"])
        winners.append(next(r for r in selected if r["candidate_id"] == best_score["candidate_id"]))

    replicated = list(screen)
    for seed in (1, 2):
        replicated.extend(_run_boltz_batch(winners, seed=seed, label=f"confirm_seed{seed}"))

    summaries: list[dict] = []
    for winner in winners:
        scores = [s for s in replicated if s["candidate_id"] == winner["candidate_id"]]
        values = [s["binder_target_iptm"] for s in scores]
        summaries.append({
            **winner,
            "boltz_replicates": scores,
            "boltz_binder_target_iptm_mean": round(statistics.mean(values), 4),
            "boltz_binder_target_iptm_range": round(max(values) - min(values), 4),
            "claim_limit": "triage only; G1 failed for repeat proteins, so this does not establish epitope recovery",
        })

    payload = {
        "experiment": "matched RFdiffusion3 constrained-vs-free binder design on a 5-chain GvpA patch",
        "arms": {"constrained": "28 cross-species GvpC-derived hotspots", "free": "center-of-mass origin, no hotspots"},
        "counts": {"backbones_per_arm": 8, "sequences_per_backbone": len(rows) // 16,
                   "candidates": len(rows), "boltz_screened": len(selected), "boltz_confirmed": len(winners)},
        "boltz_policy": "screen two distinct backbones per arm at seed 0; confirm the top score per arm at seeds 1 and 2",
        "boltz_use_msa": False,
        "screen": screen,
        "winners": summaries,
    }
    _write_json(out / "pipeline_results.json", payload)
    fasta = "".join(f">{row['candidate_id']} arm={row['arm']} iptm_mean={row['boltz_binder_target_iptm_mean']}\n{row['sequence']}\n" for row in summaries)
    (out / "best_binders.fasta").write_text(fasta)
    print(f"wrote {out / 'pipeline_results.json'} and {out / 'best_binders.fasta'}")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="stage", required=True)
    mpnn = sub.add_parser("mpnn")
    mpnn.add_argument("--n-per-backbone", type=int, default=4)
    sub.add_parser("compile")
    boltz = sub.add_parser("boltz")
    boltz.add_argument("--top-per-arm", type=int, default=2)
    args = parser.parse_args()

    if args.stage == "mpnn":
        run_mpnn(args.n_per_backbone)
    elif args.stage == "compile":
        compile_candidates()
    elif args.stage == "boltz":
        run_boltz(args.top_per_arm)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
