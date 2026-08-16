"""Sequence and structure retrieval, cached on disk, every result provenanced.

Everything here is CPU-only and free: it runs in-process through proto-tools
(`pdb-fetch-fasta`, `uniprot-fetch`) or fetches coordinates straight from RCSB.
Nothing in this module bills a Modal account, which is why the whole input
resolution phase can run before any budget decision is made.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

CACHE = Path(__file__).resolve().parent / ".cache"


def _cache_path(kind: str, key: str) -> Path:
    digest = hashlib.sha256(key.encode()).hexdigest()[:16]
    return CACHE / kind / f"{digest}.json"


def _cached(kind: str, key: str, produce) -> dict[str, Any]:
    path = _cache_path(kind, key)
    if path.exists():
        return json.loads(path.read_text())
    value = produce()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return value


def _run_tool(tool_key: str, inputs: dict[str, Any], config: dict[str, Any] | None = None) -> dict[str, Any]:
    from proto_tools.mcp.tools import run_tool

    result = run_tool(tool_key, inputs=inputs, config=config, device="local")
    if not result.get("ok", True):
        raise RuntimeError(f"{tool_key} failed: {result.get('error')}")
    return result["result"]


# --------------------------------------------------------------------------
# PDB


def pdb_chains(pdb_id: str) -> list[dict[str, Any]]:
    """Deposited chain sequences, as RCSB reports them."""

    def produce() -> dict[str, Any]:
        out = _run_tool("pdb-fetch-fasta", {"pdb_id": pdb_id})
        return {"chains": out["chains"], "source_url": out.get("source_url")}

    return _cached("pdb_fasta", pdb_id, produce)["chains"]


def pdb_chain_by_keyword(pdb_id: str, *keywords: str) -> dict[str, Any]:
    """Pick the deposited chain whose header mentions all of `keywords`."""
    wanted = [k.lower() for k in keywords]
    for chain in pdb_chains(pdb_id):
        header = chain["header"].lower()
        if all(k in header for k in wanted):
            return chain
    raise LookupError(f"{pdb_id}: no chain header matching {keywords}")


def pdb_residue_map(pdb_id: str, chain_id: str) -> dict[str, Any]:
    """Author-numbered residues actually present in the deposited coordinates.

    Needed because the literature specifies binding-pocket positions in author
    numbering, which is not the index into the SEQRES sequence.
    """

    def produce() -> dict[str, Any]:
        import biotite.database.rcsb as rcsb
        import biotite.structure.io.pdbx as pdbx
        from biotite.sequence import ProteinSequence
        from biotite.structure import get_residues

        CACHE.mkdir(parents=True, exist_ok=True)
        path = rcsb.fetch(pdb_id, "cif", CACHE / "cif")
        array = pdbx.get_structure(pdbx.CIFFile.read(path), model=1, use_author_fields=True)
        array = array[~array.hetero]
        array = array[array.chain_id == chain_id]
        if array.array_length() == 0:
            raise LookupError(f"{pdb_id} has no polymer chain {chain_id}")
        ids, names = get_residues(array)
        letters = []
        for name in names:
            try:
                letters.append(ProteinSequence.convert_letter_3to1(name))
            except Exception:
                letters.append("X")
        return {
            "residue_ids": [int(i) for i in ids],
            "sequence": "".join(letters),
            "source": f"rcsb:{pdb_id}:{chain_id} (mmCIF, author numbering)",
        }

    return _cached("pdb_residues", f"{pdb_id}:{chain_id}", produce)


# --------------------------------------------------------------------------
# UniProt


def uniprot(accession: str) -> dict[str, Any]:
    def produce() -> dict[str, Any]:
        out = _run_tool(
            "uniprot-fetch",
            {"uniprot_id": accession},
            {"fields": ["accession", "sequence", "gene_names", "ft_chain", "ft_propep", "ft_signal", "organism_name"]},
        )
        features = [
            {
                "type": f["type"],
                "start": f["location"]["start"]["value"],
                "end": f["location"]["end"]["value"],
                "description": f.get("description", ""),
            }
            for f in out["raw_entry"].get("features", [])
            if f["type"] in {"Chain", "Propeptide", "Signal", "Peptide"}
        ]
        return {"accession": accession, "sequence": out["sequence"], "features": features}

    return _cached("uniprot", accession, produce)


def uniprot_region(accession: str, feature_type: str, description: str | None = None) -> tuple[str, str]:
    """Return (sequence, provenance) for a named UniProt feature region."""
    entry = uniprot(accession)
    for feature in entry["features"]:
        if feature["type"] != feature_type:
            continue
        if description and description.lower() not in feature["description"].lower():
            continue
        start, end = feature["start"], feature["end"]
        return (
            entry["sequence"][start - 1 : end],
            f"uniprot:{accession} {feature_type}[{start}-{end}]"
            + (f" '{feature['description']}'" if feature["description"] else ""),
        )
    raise LookupError(f"{accession}: no {feature_type} feature matching {description!r}")


# --------------------------------------------------------------------------
# small sequence utilities


TERMINAL_WINDOW = 30
_N_TAG = re.compile(r"^M?[A-Z]{0,8}?H{5,}[GSA]{0,8}")
_C_TAG = re.compile(r"(?:LEVLFQGP|ENLYFQ[GS]?)?[A-Z]{0,8}?H{5,}[A-Z]{0,4}$")


@dataclass(frozen=True)
class Trim:
    sequence: str
    removed: str
    note: str


def strip_purification_tags(sequence: str) -> Trim:
    """Remove His / protease-site purification tags at either terminus.

    Deposited constructs carry them, they are disordered, and they are not part
    of the molecule whose binding was measured. Only the first and last
    `TERMINAL_WINDOW` residues are considered, so a histidine-rich patch in the
    middle of a real fold is never mistaken for a tag. Every removal is recorded
    on the chain, because it changes the sequence that gets folded.
    """
    removed_parts: list[str] = []
    head = _N_TAG.match(sequence[:TERMINAL_WINDOW])
    if head and "H" * 5 in head.group():
        removed_parts.append(f"N-terminal {head.group()}")
        sequence = sequence[head.end() :]
    tail_window = sequence[-TERMINAL_WINDOW:]
    tail = _C_TAG.search(tail_window)
    if tail and "H" * 5 in tail.group():
        cut = len(sequence) - len(tail_window) + tail.start()
        removed_parts.append(f"C-terminal {sequence[cut:]}")
        sequence = sequence[:cut]
    if not removed_parts:
        return Trim(sequence, "", "no terminal purification tag found")
    return Trim(sequence, " + ".join(removed_parts), "removed purification tag: " + "; ".join(removed_parts))


def identity(a: str, b: str) -> float:
    """Global-alignment sequence identity, for checking the case set's claims."""
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal, get_sequence_identity

    matrix = SubstitutionMatrix.std_protein_matrix()
    alignment = align_optimal(ProteinSequence(a), ProteinSequence(b), matrix, gap_penalty=(-10, -1))[0]
    return float(get_sequence_identity(alignment))
