"""Score each submitted binder by RF3 interface ipTM against its GvpA target.

Binder sequences are read from the CA records of the submissions and refolded
from sequence alone (no templates) together with the GvpA copies of the target
assembly. The reported score is the best binder-to-GvpA chain-pair ipTM; the
GvpC repeat in 8GBS is poly-UNK in the deposition and so has no sequence to
fold. Scores are dumped to /logs/verifier/scores.json.
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
MIN_IPTM = 0.60

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


def binder_interface_iptm(output, target_chain_ids: str) -> float:
    """Best chain-pair ipTM between the binder and any target chain."""
    summary = output.summary_confidences
    chain_pair_iptm = summary.get("chain_pair_iptm")
    if chain_pair_iptm is None:
        return float(summary["iptm"])
    ordered = folded_chain_ids(output)
    matrix = np.asarray(chain_pair_iptm, dtype=float)
    binder = ordered.index(BINDER_CHAIN_ID)
    return max(float(matrix[binder][ordered.index(chain_id)])
               for chain_id in target_chain_ids)


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
    scores = {
        name: {
            "iptm": binder_interface_iptm(outs[0], GVPA_CHAIN_IDS[name]),
            "global_iptm": float(outs[0].summary_confidences["iptm"]),
            "plddt": float(outs[0].summary_confidences["overall_plddt"]),
        }
        for name, outs in outputs.items()
    }
    SCORES_PATH.write_text(json.dumps(scores, indent=2))
    return scores


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
    """The refolded binder-target complex is a confident interface."""
    iptm = interface_scores[target_name]["iptm"]
    assert iptm >= MIN_IPTM, f"{target_name} binder ipTM {iptm:.3f} < {MIN_IPTM}"
