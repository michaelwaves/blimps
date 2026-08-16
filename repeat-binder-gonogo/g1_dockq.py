"""Close G1 locally: DockQ of the predicted interface against the deposit.

    .venv/bin/python g1_dockq.py

G1's pass condition is written in DockQ, which is absent from the proto-tools
catalogue — so `gates.g1_interface_recovery` reports `not_run` rather than
silently swapping in a different metric. That is the honest default, but it is
not the only option: DockQ is a standalone CPU package, and the gate can be
closed properly by installing it and comparing each predicted complex to its
deposit.

That is what this does. It re-folds each Tier C complex once (seed 0) with the
structure kept rather than discarded — the scorer throws structures away because
it only ever needed the ipTM — and runs DockQ against the deposited coordinates.

One fold per complex, ~8 GPU-minutes total. The scores are not re-used for any
other gate: this answers G1 only.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import statistics as st
import tempfile
from pathlib import Path

from gonogo import caseset, registry

REPO = Path(__file__).resolve().parent
OUT = REPO / "bench" / "repeat-binder-gonogo" / "results"
WORK = REPO / "logs" / "g1_structures"

# G1, from cases.yaml: DockQ >= 0.49 on >= 70% of positives; median iRMSD <= 2.0 A
DOCKQ_PASS = 0.49
FRACTION_REQUIRED = 0.70
IRMSD_MEDIAN_MAX = 2.0


def parse_source(source: str) -> tuple[str, list[str]]:
    """'pdb:7YCO:B' -> ('7YCO', ['B']); 'pdb:4CJ2:C/D' -> ('4CJ2', ['C', 'D']).

    Every copy is kept. A 2:2 assembly like 4CJ2 has two biological pairs and two
    pairs that never touch, and picking the first of each list lands on a
    non-contacting pair half the time — which reads as a failed prediction rather
    than as the chain-selection bug it is.
    """
    _, pdb_id, chains = source.split(":")
    return pdb_id, chains.split("/")


def native_two_chains(pdb_id: str, binder_chain: str, target_chain: str, dest: Path) -> Path:
    """The deposit, cut down to the two chains the prediction actually contains."""
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.pdbx as pdbx

    src = rcsb.fetch(pdb_id, "cif", registry.CACHE / "cif")
    cif = pdbx.CIFFile.read(src)
    # DockQ's mmCIF parser reads B_iso_or_equiv and occupancy off every atom;
    # biotite drops both unless they are asked for by name.
    array = pdbx.get_structure(cif, model=1, use_author_fields=True,
                               extra_fields=["b_factor", "occupancy"])
    array = array[~array.hetero]
    array = array[(array.chain_id == binder_chain) | (array.chain_id == target_chain)]
    if array.array_length() == 0:
        raise LookupError(f"{pdb_id}: no atoms for chains {binder_chain}/{target_chain}")
    out = pdbx.CIFFile()
    pdbx.set_structure(out, array)
    out.write(dest)
    return dest


def predict_structure(spec, dest_dir: Path) -> Path:
    """One fold, seed 0, structure kept. Same config the scorer used."""
    from proto_tools.mcp.tools import run_tool

    dest_dir.mkdir(parents=True, exist_ok=True)
    # A fold already on disk is a fold already paid for.
    existing = sorted(dest_dir.glob("*structure*"))
    if existing:
        return as_cif(existing[0])
    result = run_tool(
        "boltz2-prediction",
        inputs={"complexes": [{"chains": [
            {"id": "B", "sequence": spec.binder.sequence, "entity_type": "protein"},
            {"id": "T", "sequence": spec.target.sequence, "entity_type": "protein"},
        ]}]},
        config={"seed": 0, "recycling_steps": 3, "use_msa": True, "diffusion_samples": 1},
        device="modal",
        output_dir=str(dest_dir),
    )
    if not result.get("ok", True):
        raise RuntimeError(f"{spec.case_id}: {result.get('error')}")
    written = sorted(dest_dir.glob("*structure*"))
    if not written:
        raise RuntimeError(f"{spec.case_id}: run_tool wrote no structure into {dest_dir}")
    return as_cif(written[0])


def as_cif(path: Path) -> Path:
    """proto-tools writes mmCIF under a .txt name; DockQ picks its parser by suffix."""
    if path.suffix == ".cif":
        return path
    target = path.with_suffix(".cif")
    if not target.exists():
        shutil.copyfile(path, target)
    return target


def run_dockq(model: Path, native: Path, binder_chain: str, target_chain: str) -> dict:
    """DockQ 2.x CLI. Model chains are always B/T; the native keeps its author ids."""
    report = model.parent / "dockq.json"
    cmd = [str(REPO / ".venv" / "bin" / "DockQ"), str(model), str(native)]
    # --mapping is one character per chain, so it cannot express a multi-character
    # author id like 7QNP's 'AAA'. Let DockQ search for the pairing in that case.
    if len(binder_chain) == 1 and len(target_chain) == 1:
        cmd += ["--mapping", f"BT:{binder_chain}{target_chain}"]
    cmd += ["--json", str(report)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not report.exists():
        return {"error": (proc.stderr or proc.stdout).strip()[:400]}
    try:
        payload = json.loads(report.read_text())
    except json.JSONDecodeError:
        return {"error": f"unparseable DockQ report at {report}"}
    best = next(iter(payload.get("best_result", {}).values()), {})
    return {
        "DockQ": best.get("DockQ"),
        "iRMSD": best.get("iRMSD"),
        "LRMSD": best.get("LRMSD"),
        "fnat": best.get("fnat"),
    }


def main() -> int:
    if shutil.which("DockQ") is None and not (REPO / ".venv" / "bin" / "DockQ").exists():
        print("DockQ is not installed. `uv pip install DockQ` first.")
        return 2

    specs, _ = caseset.tier_c()
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []

    for spec in specs:
        pdb_id, binder_chains = parse_source(spec.binder.source)
        _, target_chains = parse_source(spec.target.source)
        pairs = [(b, t) for b in binder_chains for t in target_chains]
        print(f"\n=== {spec.case_id} ({pdb_id}: {len(pairs)} chain pairing(s) to try) ===", flush=True)

        case_dir = WORK / spec.case_id
        try:
            model = predict_structure(spec, case_dir)
            attempts = []
            for binder_chain, target_chain in pairs:
                native = native_two_chains(pdb_id, binder_chain, target_chain,
                                           case_dir / f"{pdb_id}_{binder_chain}{target_chain}_native.cif")
                got = run_dockq(model, native, binder_chain, target_chain)
                got["pairing"] = f"B->{binder_chain}, T->{target_chain}"
                attempts.append(got)
            # Best over the copies: the biological pair is the one that has an interface.
            usable = [a for a in attempts if a.get("DockQ") is not None]
            scores = max(usable, key=lambda a: a["DockQ"]) if usable else attempts[0]
            if len(pairs) > 1:
                scores["pairings_tried"] = [
                    {"pairing": a["pairing"], "DockQ": a.get("DockQ"), "error": a.get("error")} for a in attempts
                ]
        except Exception as exc:  # a failure here is a result, not a crash
            scores = {"error": f"{type(exc).__name__}: {exc}"}

        scores["case_id"] = spec.case_id
        scores["pdb"] = pdb_id
        rows.append(scores)
        if "error" in scores:
            print(f"  FAILED: {scores['error']}", flush=True)
        else:
            print(f"  DockQ={scores['DockQ']:.3f}  iRMSD={scores['iRMSD']:.2f}A  fnat={scores['fnat']:.3f}", flush=True)

    scored = [r for r in rows if r.get("DockQ") is not None]
    verdict = {"gate_id": "G1"}
    if not scored:
        verdict |= {"status": "not_run",
                    "evidence": "DockQ is installed but no complex produced a comparable pair; "
                                f"{len(rows)} attempted, all errored"}
    else:
        acceptable = [r for r in scored if r["DockQ"] >= DOCKQ_PASS]
        fraction = len(acceptable) / len(scored)
        median_irmsd = st.median(r["iRMSD"] for r in scored)
        ok = fraction >= FRACTION_REQUIRED and median_irmsd <= IRMSD_MEDIAN_MAX
        verdict |= {
            "status": "pass" if ok else "fail",
            "evidence": (
                f"DockQ >= {DOCKQ_PASS} on {len(acceptable)}/{len(scored)} Tier C positives "
                f"({fraction:.0%}, needs {FRACTION_REQUIRED:.0%}); median iRMSD {median_irmsd:.2f} A "
                f"(needs <= {IRMSD_MEDIAN_MAX})"
            ),
            "detail": {"fraction_acceptable": round(fraction, 3),
                       "median_irmsd": round(median_irmsd, 3),
                       "n_scored": len(scored), "n_attempted": len(rows)},
        }

    out = OUT / "g1_dockq.json"
    out.write_text(json.dumps({
        "note": "G1 closed locally with the standalone DockQ package (absent from proto-tools). "
                "One fold per complex at seed 0; structures kept rather than discarded.",
        "tool": "DockQ 2.1.3 (local, CPU) over boltz2-prediction structures from modal",
        "per_case": rows,
        "verdict": verdict,
    }, indent=2) + "\n")

    print(f"\n== G1: {verdict['status']} ==\n{verdict['evidence']}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
