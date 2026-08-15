"""Smoke-test stand-in: emits a generic helical bundle as both binders.

Not a real solution - it exists to exercise the environment and verifier
end to end. It is expected to fail the ipTM gate.
"""

from pathlib import Path

BINDER_SEQUENCE = "EELLKRAEELAKKAEELLKRAEELAKKAEELLKRAEELAKKG"
ONE_TO_THREE = {"A": "ALA", "E": "GLU", "K": "LYS",
                "L": "LEU", "R": "ARG", "G": "GLY"}
OUTPUTS_DIR = Path("/logs/outputs")


def main() -> None:
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    for target_name in ("7r1c", "8gbs"):
        (OUTPUTS_DIR / f"binder_{target_name}.pdb").write_text(extended_chain_pdb())


def extended_chain_pdb() -> str:
    lines = [
        f"ATOM  {index:5d}  CA  {ONE_TO_THREE[residue]} B{index:4d}    "
        f"{index * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C"
        for index, residue in enumerate(BINDER_SEQUENCE, 1)
    ]
    return "\n".join([*lines, "TER", "END"]) + "\n"


if __name__ == "__main__":
    main()
