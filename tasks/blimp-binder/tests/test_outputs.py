"""Score each submitted binder by RF3 interface ipTM against its GvpA target.

Binder sequences are read from the CA records of the submissions and refolded
from sequence alone (no templates) in complex with the target chain. Scores are
dumped to /logs/verifier/scores.json.
"""

import json
import os
from pathlib import Path

import pytest

OUTPUTS_DIR = Path(os.environ.get("OUTPUTS_DIR", "/logs/outputs"))
SCORES_PATH = Path(os.environ.get("SCORES_PATH", "/logs/verifier/scores.json"))

TARGET_SEQUENCES = {
    "7r1c": "SIQKSTNSSSLAEVIDRILDKGIVIDAFARVSVVGIEILTIEARVVIASVDTWLRYAEAVGLLRD",
    "8gbs": "AVEKTNSSSSLAEVIDRILDKGIVIDAWVRVSLVGIELLAIEARIVIASVETYLKYAEAVGLTQS",
}
TARGET_CHAIN_ID = "A"
BINDER_CHAIN_ID = "B"
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


def binder_path(target_name: str) -> Path:
    return OUTPUTS_DIR / f"binder_{target_name}.pdb"


def read_binder_sequence(target_name: str) -> str:
    """Extract the binder sequence from CA records, in file order."""
    path = binder_path(target_name)
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


@pytest.fixture(scope="session")
def interface_scores() -> dict:
    """Fold each binder with its target and return per-target ipTM."""
    from rf3.inference_engines.rf3 import RF3InferenceEngine
    from rf3.utils.inference import InferenceInput

    inputs = [
        InferenceInput.from_json_dict({
            "name": target_name,
            "components": [
                {"chain_id": TARGET_CHAIN_ID, "seq": target_sequence},
                {"chain_id": BINDER_CHAIN_ID, "seq": read_binder_sequence(target_name)},
            ],
        })
        for target_name, target_sequence in TARGET_SEQUENCES.items()
    ]
    engine = RF3InferenceEngine(ckpt_path="rf3", verbose=False, **RF3_SETTINGS)
    outputs = engine.run(inputs=inputs, out_dir=None)
    scores = {
        name: float(outs[0].summary_confidences["iptm"])
        for name, outs in outputs.items()
    }
    SCORES_PATH.write_text(json.dumps(scores, indent=2))
    return scores


@pytest.mark.parametrize("target_name", TARGET_SEQUENCES)
def test_binder_format(target_name: str):
    """The binder is a single chain of standard amino acids within the length window."""
    length = len(read_binder_sequence(target_name))
    assert MIN_BINDER_LENGTH <= length <= MAX_BINDER_LENGTH, (
        f"binder for {target_name} has {length} residues, expected "
        f"{MIN_BINDER_LENGTH}-{MAX_BINDER_LENGTH}"
    )


@pytest.mark.parametrize("target_name", TARGET_SEQUENCES)
def test_binder_interface_confidence(target_name: str, interface_scores: dict):
    """The refolded binder-target complex is a confident interface."""
    iptm = interface_scores[target_name]
    assert iptm >= MIN_IPTM, f"{target_name} ipTM {iptm:.3f} < {MIN_IPTM}"
