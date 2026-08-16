"""How wide a gas vesicle does the generated patch imply?

    .venv/bin/python gvpa_curvature.py

The 5-rib GvpA patch was produced by a generative model. Sequence and fold can
look right while the *curvature* is wrong, and curvature is the one property a
cylindrical shell cannot fake: a patch cut from an 85 nm gas vesicle has to bow
like an 85 nm cylinder.

METHOD, AND WHY THE FIRST ONE WAS WRONG

The first version swept cylinder axes and fitted circles to the projection. It
returned radius 63 A / 137 deg of arc, and rendering the patch in PyMOL showed
that to be nonsense — the structure is visibly a flat sheet, and the fit had
locked onto internal scatter rather than global curvature. A 137 deg arc would
look like a letter C; the render looks like a plank.

This version measures the sagitta instead, which is hard to fool. Take a rib,
take the subunit centroids along it, fit a straight line, and measure the
maximum perpendicular deviation. For a chord of length L on a circle of radius
R, the sagitta is s = L^2 / 8R, so R = L^2 / 8s. Flat sheet -> s ~ 0 -> R
diverges. No axis search, no projection, nothing to overfit.

Reference: Ana GVs are 85 +/- 4 nm in diameter with a ~2.8 nm shell
(Dutka et al. 2023, PMC10185304 L22, L83).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile

REPO = Path(__file__).resolve().parent
OUT = REPO / "designs" / "gvpa"

ANA_DIAMETER_NM = 85.0
ANA_SHELL_NM = 2.8


def sagitta(points: np.ndarray) -> tuple[float, float]:
    """Max perpendicular deviation from the straight line through the endpoints,
    and the chord length. Both in the same units as `points`."""
    a, b = points[0], points[-1]
    chord = b - a
    L = float(np.linalg.norm(chord))
    if L < 1e-6:
        return 0.0, 0.0
    u = chord / L
    rel = points - a
    perp = rel - np.outer(rel @ u, u)
    return float(np.linalg.norm(perp, axis=1).max()), L


def main() -> int:
    patch = PDBFile.read(str(REPO / "GvpA_surface_patch_5ribs_x11.pdb")).get_structure(model=1)
    ca = patch[patch.atom_name == "CA"]
    chains = sorted(set(ca.chain_id))
    cent = np.array([ca.coord[ca.chain_id == c].mean(axis=0) for c in chains], dtype=float)

    # Ribs: 55 subunits as 5 ribs of 11. Group by projection on the principal axis
    # that separates ribs — the one with 5 clear clusters.
    centred = cent - cent.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    # Axis 0 (spread 206 A) separates into five clumps of eleven — that is the
    # rib-stacking direction. Axis 1 (126 A) runs along a rib. Axis 2 (11 A) is
    # the thin direction, and the whole out-of-plane bow has to live inside it.
    # An earlier version had 0 and 1 swapped, which sliced across the ribs and
    # produced a spurious symmetric curvature profile.
    proj = centred @ vt[0]
    order = np.argsort(proj)
    ribs = [order[i * 11:(i + 1) * 11] for i in range(5)]

    print(f"{len(chains)} subunits grouped into {len(ribs)} ribs of 11\n")
    print(f"{'rib':>4} {'chord L (A)':>12} {'sagitta s (A)':>14} {'implied R (A)':>14} {'diameter (nm)':>14}")
    print("-" * 62)

    rows = []
    for i, idx in enumerate(ribs):
        pts = cent[idx][np.argsort(centred[idx] @ vt[1])]   # order along the rib
        s, L = sagitta(pts)
        R = (L * L) / (8 * s) if s > 1e-3 else float("inf")
        d_nm = 2 * R / 10
        rows.append({"rib": i, "chord_A": round(L, 1), "sagitta_A": round(s, 2),
                     "implied_radius_A": None if np.isinf(R) else round(R, 1),
                     "implied_diameter_nm": None if np.isinf(d_nm) else round(d_nm, 1)})
        print(f"{i:>4} {L:>12.1f} {s:>14.2f} "
              f"{'flat' if np.isinf(R) else f'{R:>14.1f}'} "
              f"{'-' if np.isinf(d_nm) else f'{d_nm:>14.1f}'}")

    # What sagitta would a real Ana GV fragment show over this chord?
    mean_L = float(np.mean([r["chord_A"] for r in rows]))
    R_ana = ANA_DIAMETER_NM * 10 / 2
    expected_s = mean_L**2 / (8 * R_ana)
    observed_s = float(np.mean([r["sagitta_A"] for r in rows]))

    print(f"\nover a {mean_L:.0f} A rib chord:")
    print(f"  sagitta expected for an 85 nm Ana GV : {expected_s:6.2f} A")
    print(f"  sagitta observed in the generated patch: {observed_s:6.2f} A")

    # Shell thickness is a second, independent check on the same structure.
    thickness = float(np.percentile(np.linalg.norm(ca.coord - ca.coord.mean(axis=0), axis=1), 95) -
                      np.percentile(np.linalg.norm(ca.coord - ca.coord.mean(axis=0), axis=1), 5))
    verdict = ("consistent with an 85 nm GV" if 0.5 * expected_s < observed_s < 2 * expected_s else
               "FLATTER than any real GV — the patch carries essentially no curvature"
               if observed_s < 0.5 * expected_s else
               "MORE curved than a real Ana GV")
    print(f"\n  verdict: {verdict}")

    payload = {
        "method": "sagitta of subunit centroids along each rib; R = L^2/8s. Replaces an "
                  "earlier cylinder-axis sweep that returned radius 63 A / 137 deg of arc — "
                  "PyMOL rendering showed the patch is a flat sheet, so that fit was "
                  "locking onto internal scatter, not global curvature.",
        "reference": {"source": "Dutka et al. 2023, PMC10185304 L22 and L83",
                      "ana_gv_diameter_nm": ANA_DIAMETER_NM, "ana_shell_thickness_nm": ANA_SHELL_NM},
        "per_rib": rows,
        "mean_chord_A": round(mean_L, 1),
        "sagitta_observed_A": round(observed_s, 3),
        "sagitta_expected_for_85nm_A": round(expected_s, 3),
        "verdict": verdict,
    }
    (OUT / "curvature.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nwrote {OUT / 'curvature.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
