"""Score each submitted binder by RF3 interface PAE against its GvpA target.

Binder sequences are read from the CA records of the submissions and refolded
from sequence alone (no templates) together with the GvpA copies of the target
assembly. RF3 reports no per-chain-pair ipTM - only a global one, which is
dominated by the GvpA-GvpA contacts it cannot reassemble from sequence - so the
score is the symmetrized chain-pair PAE between the binder and the GvpA copy it
binds best. The GvpC repeat in 8GBS is poly-UNK in the deposition and so has no
sequence to fold. Scores are dumped to /logs/verifier/scores.json.
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

OUTPUTS_DIR = Path(os.environ.get("OUTPUTS_DIR", "/logs/outputs"))
SCORES_PATH = Path(os.environ.get("SCORES_PATH", "/logs/verifier/scores.json"))

GVPA_SEQUENCES = {
    "7r1c": "SIQKSTNSSSLAEVIDRILDKGIVIDAFARVSVVGIEILTIEARVVIASVDTWLRYAEAVGLLRD",
    "8gbs": "AVEKTNSSSSLAEVIDRILDKGIVIDAWVRVSLVGIELLAIEARIVIASVETYLKYAEAVGLTQS",
}
GVPA_CHAIN_IDS = {"7r1c": "ABCDE", "8gbs": "ABCD"}
BINDER_CHAIN_ID = "Z"
MIN_BINDER_LENGTH = 30
MAX_BINDER_LENGTH = 150
MAX_INTERFACE_PAE_ANGSTROM = 10.0  # CALIBRATE against a real designed binder

RF3_SETTINGS = dict(n_recycles=10, diffusion_batch_size=1, num_steps=50, seed=0)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def read_binder_sequence(target_name: str) -> str:
    """Extract the binder sequence from CA records, in file order."""
    path = OUTPUTS_DIR / f"binder_{target_name}.pdb"
    if not path.exists():
        raise ValueError(f"submission not found: {path}")
    residues = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            residue_name = line[17:20].strip()
            if residue_name not in THREE_TO_ONE:
                raise ValueError(f"non-standard residue {residue_name!r} in {path}")
            residues.append(THREE_TO_ONE[residue_name])
    if not residues:
        raise ValueError(f"no CA records found in {path}")
    return "".join(residues)


def folded_chain_ids(output) -> list[str]:
    """Chain IDs in tokenization order; RF3 suffixes them with an entity index."""
    ordered = []
    for chain_id in output.confidences["token_chain_ids"]:
        base = chain_id.split("_", 1)[0]
        if base not in ordered:
            ordered.append(base)
    return ordered


def binder_interface_paes(output, target_chain_ids: str) -> dict[str, float]:
    """Symmetrized chain-pair PAE between the binder and each target chain."""
    ordered = folded_chain_ids(output)
    matrix = np.asarray(output.summary_confidences["chain_pair_pae"], dtype=float)
    assert matrix.shape == (len(ordered), len(ordered)), (
        f"chain_pair_pae is {matrix.shape}, expected {len(ordered)} chains"
    )
    binder = ordered.index(BINDER_CHAIN_ID)
    return {
        chain_id: float(matrix[binder][index] + matrix[index][binder]) / 2
        for chain_id, index in ((c, ordered.index(c)) for c in target_chain_ids)
    }


@pytest.fixture(scope="session")
def interface_scores() -> dict:
    """Fold each binder with the GvpA copies of its target and score the interface."""
    from rf3.inference_engines.rf3 import RF3InferenceEngine
    from rf3.utils.inference import InferenceInput

    inputs = [
        InferenceInput.from_json_dict({
            "name": target_name,
            "components": [
                *({"chain_id": chain_id, "seq": GVPA_SEQUENCES[target_name]}
                  for chain_id in GVPA_CHAIN_IDS[target_name]),
                {"chain_id": BINDER_CHAIN_ID, "seq": read_binder_sequence(target_name)},
            ],
        })
        for target_name in GVPA_SEQUENCES
    ]
    engine = RF3InferenceEngine(ckpt_path="rf3", verbose=False, **RF3_SETTINGS)
    outputs = engine.run(inputs=inputs, out_dir=None)
    scores = {name: score_output(outs[0], GVPA_CHAIN_IDS[name])
              for name, outs in outputs.items()}
    SCORES_PATH.write_text(json.dumps(scores, indent=2))
    return scores


def score_output(output, target_chain_ids: str) -> dict:
    """Gate metric plus the context needed to recalibrate it."""
    interface_paes = binder_interface_paes(output, target_chain_ids)
    summary = output.summary_confidences
    return {
        "interface_pae": min(interface_paes.values()),
        "interface_pae_per_chain": interface_paes,
        "global_iptm": float(summary["iptm"]),
        "plddt": float(summary["overall_plddt"]),
        "has_clash": bool(summary["has_clash"]),
    }


@pytest.mark.parametrize("target_name", GVPA_SEQUENCES)
def test_binder_format(target_name: str):
    """The binder is a single chain of standard amino acids within the length window."""
    length = len(read_binder_sequence(target_name))
    assert MIN_BINDER_LENGTH <= length <= MAX_BINDER_LENGTH, (
        f"binder for {target_name} has {length} residues, expected "
        f"{MIN_BINDER_LENGTH}-{MAX_BINDER_LENGTH}"
    )


@pytest.mark.parametrize("target_name", GVPA_SEQUENCES)
def test_binder_interface_confidence(target_name: str, interface_scores: dict):
    """The binder packs against a GvpA copy with a confident interface."""
    interface_pae = interface_scores[target_name]["interface_pae"]
    assert interface_pae <= MAX_INTERFACE_PAE_ANGSTROM, (
        f"{target_name} binder interface PAE {interface_pae:.2f} A > "
        f"{MAX_INTERFACE_PAE_ANGSTROM} A"
    )
