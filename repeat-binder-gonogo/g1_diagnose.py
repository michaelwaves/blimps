"""Why did the interface predictions fail? Wrong epitope, or wrong register?

    .venv/bin/python g1_diagnose.py

G1 reports that 6 of 8 Tier C complexes miss their interface, all of them repeat
proteins. "iRMSD 17 A" says the answer is wrong; it does not say how. Two very
different failures produce the same number:

  wrong epitope   the binder is docked onto a different face of the target
                  entirely. The model has no idea where the site is.
  wrong register  the binder is on the right face, but shifted along it — the
                  wrong repeat contacts the right patch. Plausible on a repeat
                  protein, whose surface is periodic and near-degenerate.

The discriminator is which TARGET residues get contacted. Epitope error moves
that set; register error largely preserves it while moving which BINDER residues
sit opposite.

This is analysis of a finished run, not an input to a new one. Nothing here is
fed back into a prediction — doing that would leak the deposited interface into
a metric that is scored against the deposited interface.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent
WORK = REPO / "logs" / "g1_structures"
OUT = REPO / "bench" / "repeat-binder-gonogo" / "results"
CONTACT_A = 5.0   # heavy-atom contact cutoff


def load(path: Path):
    import biotite.structure.io.pdbx as pdbx
    a = pdbx.get_structure(pdbx.CIFFile.read(path), model=1, use_author_fields=True)
    return a[~a.hetero]


def chain_sequence(arr, chain_id: str):
    """One-letter sequence and the residue ids behind it, in order."""
    from biotite.sequence import ProteinSequence
    from biotite.structure import get_residues

    sub = arr[arr.chain_id == chain_id]
    ids, names = get_residues(sub)
    letters = []
    for name in names:
        try:
            letters.append(ProteinSequence.convert_letter_3to1(name))
        except Exception:
            letters.append("X")
    return "".join(letters), list(ids)


def contact_residues(arr, chain_a: str, chain_b: str) -> tuple[set, set]:
    """Residue ids on each side within CONTACT_A of the other chain."""
    A = arr[arr.chain_id == chain_a]
    B = arr[arr.chain_id == chain_b]
    if A.array_length() == 0 or B.array_length() == 0:
        return set(), set()
    d = np.linalg.norm(A.coord[:, None, :] - B.coord[None, :, :], axis=-1)
    close = d < CONTACT_A
    return set(A.res_id[close.any(axis=1)].tolist()), set(B.res_id[close.any(axis=0)].tolist())


def index_map(model_arr, model_chain, native_arr, native_chain) -> dict[int, int]:
    """Map native residue id -> model residue id by aligning the two chains.

    A predicted chain is numbered 1..N off the input sequence; a deposit carries
    author numbering with its own offset and gaps, and the input had purification
    tags stripped. Comparing residue ids across the two without aligning first
    compares nothing, which is exactly the bug this replaced.
    """
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal

    m_seq, m_ids = chain_sequence(model_arr, model_chain)
    n_seq, n_ids = chain_sequence(native_arr, native_chain)
    if not m_seq or not n_seq:
        return {}
    matrix = SubstitutionMatrix.std_protein_matrix()
    aln = align_optimal(ProteinSequence(m_seq.replace("X", "A")),
                        ProteinSequence(n_seq.replace("X", "A")),
                        matrix, gap_penalty=(-10, -1), terminal_penalty=False)[0]
    mapping = {}
    for m_i, n_i in aln.trace:
        if m_i >= 0 and n_i >= 0:
            mapping[n_ids[n_i]] = m_ids[m_i]
    return mapping


def jaccard(x: set, y: set) -> float:
    return len(x & y) / len(x | y) if (x or y) else 0.0


def main() -> int:
    from gonogo import caseset

    dq = {r["case_id"]: r for r in json.loads((OUT / "g1_dockq.json").read_text())["per_case"]}
    # Authoritative chain ids: the filename cannot express a multi-character author id
    # like 7QNP's 'AAA', so read the pairing DockQ scored, falling back to the case set.
    specs = {s.case_id: s for s in caseset.tier_c()[0]}
    rows = []

    for case_dir in sorted(WORK.iterdir()):
        if not case_dir.is_dir():
            continue
        cid = case_dir.name
        model_path = case_dir / "structures_0_structure.cif"
        if not model_path.exists() or cid not in specs:
            continue
        pairing = dq.get(cid, {}).get("pairing")
        if pairing:                                    # "B->C, T->B"
            nb = pairing.split("B->")[1].split(",")[0].strip()
            nt = pairing.split("T->")[1].strip()
        else:
            nb = specs[cid].binder.source.split(":")[2].split("/")[0]
            nt = specs[cid].target.source.split(":")[2].split("/")[0]
        native_path = case_dir / f"{specs[cid].binder.source.split(':')[1]}_{nb}{nt}_native.cif"
        if not native_path.exists():
            cands = sorted(case_dir.glob("*_native.cif"))
            if not cands:
                continue
            native_path = cands[0]

        model, native = load(model_path), load(native_path)
        m_binder, m_target = contact_residues(model, "B", "T")
        n_binder, n_target = contact_residues(native, nb, nt)

        # Project the deposit's interface onto the model's numbering before comparing.
        t_map = index_map(model, "T", native, nt)
        b_map = index_map(model, "B", native, nb)
        n_target_in_model = {t_map[r] for r in n_target if r in t_map}
        n_binder_in_model = {b_map[r] for r in n_binder if r in b_map}

        row = {
            "case_id": cid,
            "DockQ": round(dq.get(cid, {}).get("DockQ") or 0.0, 4),
            "n_target_epitope": len(n_target_in_model),
            "n_predicted_target_contacts": len(m_target),
            "target_epitope_jaccard": round(jaccard(m_target, n_target_in_model), 3),
            "target_epitope_recovered": round(
                len(m_target & n_target_in_model) / len(n_target_in_model), 3) if n_target_in_model else None,
            "binder_paratope_jaccard": round(jaccard(m_binder, n_binder_in_model), 3),
        }
        n_target = n_target_in_model
        # Classify only where there is a real interface on both sides.
        if not n_target or not m_target:
            row["diagnosis"] = "no comparable interface"
        elif row["DockQ"] and row["DockQ"] >= 0.49:
            row["diagnosis"] = "recovered"
        elif row["target_epitope_jaccard"] >= 0.30:
            row["diagnosis"] = "right epitope, wrong register/orientation"
        elif row["target_epitope_jaccard"] > 0.0:
            row["diagnosis"] = "epitope partially overlapping — drifted off site"
        else:
            row["diagnosis"] = "wrong epitope entirely — different face of the target"
        rows.append(row)
        print(f"  {cid:9} DockQ={str(row['DockQ']):>6}  target-epitope J={row['target_epitope_jaccard']:.3f}  "
              f"paratope J={row['binder_paratope_jaccard']:.3f}  {row['diagnosis']}", flush=True)

    failed = [r for r in rows if r["diagnosis"] not in ("recovered", "no comparable interface")]
    tally = {}
    for r in failed:
        tally[r["diagnosis"]] = tally.get(r["diagnosis"], 0) + 1

    payload = {
        "question": "for the complexes G1 failed, is the binder on the wrong face or the wrong part of the right face?",
        "method": f"heavy-atom contacts within {CONTACT_A} A; Jaccard overlap of the contacted TARGET residue set "
                  "between prediction and deposit. Epitope error moves that set; register error preserves it.",
        "not_a_feedback_signal": "this reads a finished run. It is not used to constrain, select or re-score any "
                                 "prediction — the deposited interface must not enter a metric measured against "
                                 "the deposited interface.",
        "per_case": rows,
        "tally_over_failures": tally,
    }
    (OUT / "g1_diagnosis.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\n  failure modes: {tally}")
    print(f"  wrote {OUT / 'g1_diagnosis.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
