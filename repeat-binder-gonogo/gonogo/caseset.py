"""Turn `cases.yaml` into concrete, scorable two-chain complexes.

Every complex records where each chain came from and every edit made to it.
Chain selection is an explicit table, not a keyword guess: a silently mispicked
chain would produce a confident, wrong, fully-provenanced number, which is the
worst failure this pipeline can have.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from . import registry

REPO = Path(__file__).resolve().parent.parent
TASK_DIR = REPO / "bench" / "repeat-binder-gonogo"
DATA = Path(__file__).resolve().parent / "data"

# Deposited chains, keyed (pdb_id) -> (binder_chain_label, target_chain_label).
# Labels are the strings RCSB puts in the FASTA header, matched case-insensitively.
COMPLEX_CHAINS: dict[str, tuple[str, str]] = {
    "9S60": ("DARPin 4m3", "Cathepsin B"),
    "5MBL": ("DARPin 81", "Cathepsin B"),
    "5MBM": ("DARPin 8h6", "Cathepsin B"),
    "1MJ0": ("SANK E3_5", ""),
    "7YCO": ("Repebody (A6)", "Spike protein S1"),
    "4CJ0": ("E12 AFFITIN", "ENDOGLUCANASE D"),
    "4CJ2": ("AFFITIN H4", "LYSOZYME C"),
    "7QNP": ("Designed Armadillo Repeat Protein", "Lysozyme"),
    "6G4J": ("alphaREP", "kinase YabT"),
    "1SVX": ("off7", "Maltose-binding"),
    "8RCI": ("DARPin C10", "p53"),
    "9SPO": ("DARPin C10", "p53"),
}

MOUSE_CATB = "P10605"
HUMAN_CATB = "P07858"


@dataclass
class Chain:
    name: str
    sequence: str
    source: str
    edits: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class ComplexSpec:
    case_id: str
    tier: str
    binder: Chain
    target: Chain
    expectation: str  # "binder" | "non_binder" | "decoy" | "unknown"
    truth: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def as_record(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "tier": self.tier,
            "expectation": self.expectation,
            "binder": {"name": self.binder.name, "len": self.binder.length, "source": self.binder.source,
                       "edits": self.binder.edits},
            "target": {"name": self.target.name, "len": self.target.length, "source": self.target.source,
                       "edits": self.target.edits},
            "truth": self.truth,
            "notes": self.notes,
        }


@dataclass
class Gap:
    """Something the harness could not resolve. Reported, never papered over."""

    what: str
    why: str
    blocks: list[str]


def deposited_chain(pdb_id: str, which: str) -> Chain:
    binder_label, target_label = COMPLEX_CHAINS[pdb_id]
    label = binder_label if which == "binder" else target_label
    chain = registry.pdb_chain_by_keyword(pdb_id, label)
    trim = registry.strip_purification_tags(chain["sequence"])
    edits = [] if not trim.removed else [trim.note]
    return Chain(
        name=chain["header"].split("|")[2] if chain["header"].count("|") >= 2 else label,
        sequence=trim.sequence,
        source=f"pdb:{pdb_id}:{'/'.join(chain['chain_ids'])}",
        edits=edits,
    )


# --------------------------------------------------------------------------
# Tier B — the gate the validation exists for


def tier_b() -> tuple[list[ComplexSpec], list[Gap]]:
    cases = yaml.safe_load((TASK_DIR / "cases.yaml").read_text())
    tier = cases["tier_b_specificity"]
    positives = {p["id"]: p for p in tier["positives"]}
    negatives = {n["id"]: n for n in tier["negatives"]}

    specs: list[ComplexSpec] = []
    gaps: list[Gap] = []

    darpin_4m3 = deposited_chain("9S60", "binder")

    mouse_mature, mouse_src = registry.uniprot_region(MOUSE_CATB, "Chain", "Cathepsin B")
    human_mature, human_src = registry.uniprot_region(HUMAN_CATB, "Chain", "Cathepsin B")
    deposited_mouse = deposited_chain("9S60", "target")

    # The negative's construct must match the positive's construct, or the
    # comparison measures the construct rather than the ortholog.
    boundary_check = (
        "identical to the 9S60 deposit"
        if deposited_mouse.sequence == mouse_mature
        else f"DIFFERS from the 9S60 deposit ({len(deposited_mouse.sequence)} vs {len(mouse_mature)} aa)"
    )

    specs.append(
        ComplexSpec(
            case_id="B1",
            tier="B",
            binder=darpin_4m3,
            target=Chain("mouse cathepsin B (mature)", mouse_mature, mouse_src),
            expectation="binder",
            truth={"KD_nM": positives["B1"]["KD_nM"], "Ki_nM": positives["B1"]["Ki_nM"], "method": "SPR + kinetics"},
            notes=[f"mature-chain construct {boundary_check}"],
        )
    )
    specs.append(
        ComplexSpec(
            case_id="B4",
            tier="B",
            binder=darpin_4m3,
            target=Chain("human cathepsin B (mature)", human_mature, human_src),
            expectation="non_binder",
            truth={"outcome": negatives["B4"]["outcome"], "selectivity_residues": negatives["B4"]["selectivity_residues"]},
            notes=["same construct boundaries as B1 — the only variable is the ortholog"],
        )
    )
    for case_id, pdb_id in (("B2", "5MBL"), ("B3", "5MBM")):
        specs.append(
            ComplexSpec(
                case_id=case_id,
                tier="B",
                binder=deposited_chain(pdb_id, "binder"),
                target=deposited_chain(pdb_id, "target"),
                expectation="binder",
                truth={"pdb": pdb_id, "note": "cross-reactive DARPin; the deposited target is HUMAN cathepsin B"},
            )
        )

    zymogen, zym_src = registry.uniprot_region(MOUSE_CATB, "Propeptide", "Activation peptide")
    specs.append(
        ComplexSpec(
            case_id="B5",
            tier="B",
            binder=darpin_4m3,
            target=Chain("mouse pro-cathepsin B", zymogen + mouse_mature, f"{zym_src} + {mouse_src}"),
            expectation="non_binder",
            truth={"outcome": negatives["B5"]["outcome"]},
            notes=["propeptide fused to the mature chain; the deposited zymogen structure is not in the case set"],
        )
    )
    try:
        e3_5 = deposited_chain("1MJ0", "binder")
        specs.append(
            ComplexSpec(
                case_id="B6",
                tier="B",
                binder=e3_5,
                target=Chain("mouse cathepsin B (mature)", mouse_mature, mouse_src),
                expectation="non_binder",
                truth={"outcome": negatives["B6"]["outcome"]},
                notes=["E3_5 sequence taken from its own deposit (1MJ0); cases.yaml cites no structure for B6"],
            )
        )
    except LookupError as exc:
        gaps.append(Gap("B6 binder sequence", str(exc), ["G2 completeness", "G5 positives pool"]))

    return specs, gaps


# --------------------------------------------------------------------------
# Tier C — scaffold generalisation, one binder/target pair per deposit


def tier_c() -> tuple[list[ComplexSpec], list[Gap]]:
    cases = yaml.safe_load((TASK_DIR / "cases.yaml").read_text())
    specs, gaps = [], []
    for entry in cases["tier_c_generalisation"]:
        pdb_id = entry["pdb"]
        if pdb_id not in COMPLEX_CHAINS:
            gaps.append(Gap(f"tier C {pdb_id}", "no chain assignment in COMPLEX_CHAINS", ["G1"]))
            continue
        try:
            specs.append(
                ComplexSpec(
                    case_id=f"C_{pdb_id}",
                    tier="C",
                    binder=deposited_chain(pdb_id, "binder"),
                    target=deposited_chain(pdb_id, "target"),
                    expectation="binder",
                    truth={"pdb": pdb_id, "resolution_A": entry["resolution_A"], "complex": entry["complex"]},
                )
            )
        except LookupError as exc:
            gaps.append(Gap(f"tier C {pdb_id}", str(exc), ["G1"]))
    return specs, gaps


# --------------------------------------------------------------------------
# Tier A — the affinity ladder, built from the closed K_D table


def tier_a(scaffold: str = "scaffold_5aei", pockets: list[str] | None = None) -> tuple[list[ComplexSpec], list[Gap]]:
    table = yaml.safe_load((DATA / "tier_a_kd.yaml").read_text())
    recipe = table["construct_recipe"][scaffold]
    pdb_id = {"scaffold_5aei": "5AEI", "scaffold_6sa8": "6SA8"}[scaffold]

    residues = registry.pdb_residue_map(pdb_id, recipe["chain"])
    index_of = {rid: i for i, rid in enumerate(residues["residue_ids"])}
    sequence = list(residues["sequence"])

    lo, hi = recipe["residue_range"]
    positions = recipe["pocket_positions"]
    missing = [p for p in positions if p not in index_of]
    if missing:
        return [], [Gap(f"tier A {pdb_id}", f"pocket positions absent from coordinates: {missing}", ["G3"])]

    wild_type_pocket = "".join(sequence[index_of[p]] for p in positions)
    gaps: list[Gap] = []
    if wild_type_pocket != table["pockets"]["arg"]["pocket_residues"]:
        gaps.append(
            Gap(
                f"tier A {pdb_id} pocket mapping",
                f"positions {positions} read {wild_type_pocket}, expected the Arg-binder pocket "
                f"{table['pockets']['arg']['pocket_residues']}",
                ["G3"],
            )
        )
        return [], gaps

    specs: list[ComplexSpec] = []
    for pocket_name, pocket in table["pockets"].items():
        if pockets and pocket_name not in pockets:
            continue
        scaffold_seq = list(sequence)
        for position, residue in zip(positions, pocket["pocket_residues"]):
            scaffold_seq[index_of[position]] = residue
        window = [i for i, rid in enumerate(residues["residue_ids"]) if lo <= rid <= hi]
        binder_seq = "".join(scaffold_seq[i] for i in window)

        for variant in pocket["variants"]:
            peptide = list(table["construct_recipe"]["base_peptide"])
            for flank in recipe["flanking_alanine_positions"]:
                peptide[flank - 1] = "A"
            peptide[recipe["variable_ligand_position"] - 1] = variant["residue"]
            specs.append(
                ComplexSpec(
                    case_id=f"A_{pocket_name}_{variant['residue']}_{pdb_id}",
                    tier="A",
                    binder=Chain(
                        f"dArmRP {pocket_name}-binder pocket ({pocket['pocket_residues']})",
                        binder_seq,
                        f"pdb:{pdb_id}:{recipe['chain']}[{lo}-{hi}]",
                        edits=[f"pocket positions {positions} -> {pocket['pocket_residues']}"],
                    ),
                    target=Chain(
                        f"peptide {''.join(peptide)}",
                        "".join(peptide),
                        f"{table['construct_recipe']['base_peptide']} per Table S1 recipe",
                        edits=[
                            f"position {recipe['variable_ligand_position']} -> {variant['residue']}",
                            f"positions {recipe['flanking_alanine_positions']} -> A "
                            f"({recipe['flanking_alanine_status']})",
                        ],
                    ),
                    expectation="binder",
                    truth={"kd_nM": variant["kd_nM"], "error_nM": variant["error_nM"], "pocket": pocket_name,
                           "scaffold": pdb_id},
                    notes=["K_D from Table S1 of PMC12202725 (L547-L629)"],
                )
            )
    return specs, gaps


# --------------------------------------------------------------------------
# Decoys — G5


def decoys(pool: list[ComplexSpec], seed: int = 0) -> list[ComplexSpec]:
    """Sequence-shuffled binders and non-cognate pairings, deterministically.

    Non-cognate pairs are drawn only from binders whose real target is in a
    different system entirely. DARPins 81 and 8h6 are cross-reactive across
    cathepsin B orthologs, so pairing either of them with the mouse enzyme
    would be a decoy that is plausibly a real binder — it is excluded.
    """
    by_id = {spec.case_id: spec for spec in pool}
    rng = random.Random(seed)
    out: list[ComplexSpec] = []

    def shuffled(sequence: str) -> str:
        letters = list(sequence)
        rng.shuffle(letters)
        return "".join(letters)

    for case_id in ("B1", "B2"):
        if case_id not in by_id:
            continue
        spec = by_id[case_id]
        for replicate in range(2):
            out.append(
                ComplexSpec(
                    case_id=f"decoy_shuffle_{case_id}_r{replicate}",
                    tier="decoy",
                    binder=Chain(
                        f"{spec.binder.name} [shuffled r{replicate}]",
                        shuffled(spec.binder.sequence),
                        spec.binder.source,
                        edits=[f"sequence shuffled, python random seed {seed}, replicate {replicate}"],
                    ),
                    target=spec.target,
                    expectation="decoy",
                    notes=["composition-preserving shuffle: same residues, no fold"],
                )
            )

    cross = [
        ("C_1SVX", "B1", "anti-MBP DARPin off7 against mouse cathepsin B"),
        ("C_4CJ2", "B1", "anti-lysozyme affitin H4 against mouse cathepsin B"),
        ("C_7YCO", "B1", "anti-RBD repebody A6 against mouse cathepsin B"),
        ("B1", "C_1SVX", "DARPin 4m3 against maltose binding protein"),
        ("B1", "C_4CJ2", "DARPin 4m3 against hen egg white lysozyme"),
    ]
    for binder_case, target_case, why in cross:
        if binder_case not in by_id or target_case not in by_id:
            continue
        out.append(
            ComplexSpec(
                case_id=f"decoy_noncognate_{binder_case}_x_{target_case}",
                tier="decoy",
                binder=by_id[binder_case].binder,
                target=by_id[target_case].target,
                expectation="decoy",
                notes=[why, "both chains are real, folded proteins; only the pairing is wrong"],
            )
        )
    return out
