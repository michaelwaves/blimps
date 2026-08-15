"""Rebuild environment/pdbs from the RCSB biological assemblies.

7R1C -> five GvpA copies (chains A-E); 8GBS -> four GvpA copies (chains A-D)
plus the modelled GvpC repeat (chain G, poly-UNK in the deposition).
Coordinates only: headers, waters and ligands are dropped.

    python authoring/build_targets.py
"""

import urllib.request
from pathlib import Path

PDBS_DIR = Path(__file__).resolve().parent.parent / "environment" / "pdbs"
ASSEMBLY_URL = "https://files.rcsb.org/download/{entry}-assembly1.cif"
CHAIN_RENAMES = {
    "7r1c": {"N": "A", "N-2": "B", "N-3": "C", "N-4": "D", "N-5": "E"},
    "8gbs": {"A1": "A", "A2": "B", "A3": "C", "A4": "D", "C": "G"},
}


def main() -> None:
    for entry, renames in CHAIN_RENAMES.items():
        assembly = fetch_assembly(entry)
        (PDBS_DIR / f"{entry}.pdb").write_text(to_pdb(assembly, renames))


def fetch_assembly(entry: str) -> list[str]:
    with urllib.request.urlopen(ASSEMBLY_URL.format(entry=entry.upper())) as response:
        text = response.read().decode()
    return [line for line in text.splitlines() if line.startswith("ATOM")]


def to_pdb(atom_lines: list[str], renames: dict[str, str]) -> str:
    records = []
    for serial, line in enumerate(atom_lines, 1):
        fields = line.split()
        element, residue_name = fields[2], fields[-4]
        atom_name = fields[3] if len(fields[3]) == 4 else f" {fields[3]}"
        chain_id, residue_number = renames[fields[-3]], int(fields[-5])
        x, y, z = (float(value) for value in fields[10:13])
        records.append(
            f"ATOM  {serial:5d} {atom_name:<4}{residue_name:>4} {chain_id}"
            f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {element:>2}"
        )
    return "\n".join([*records, "TER", "END"]) + "\n"


if __name__ == "__main__":
    main()
