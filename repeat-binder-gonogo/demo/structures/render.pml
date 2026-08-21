# Predicted-vs-deposited structural comparison for the go/no-go demo.
# Run headless: /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq demo/structures/render.pml
#
# Chain naming: predicted structures from boltz2 use B (binder) / T (target).
# Deposited PDB files use their own author chain IDs (looked up per case below).

bg_color white
set ray_opaque_background, 0
set cartoon_fancy_helices, 1
set specular, 0.15

python

import os
OUT = "demo/structures"

def render_pair(name, predicted, truth, truth_binder_chain, truth_target_chain,
                 pred_binder_color, truth_binder_color, target_color, caption):
    cmd.reinitialize()
    cmd.load(predicted, "pred")
    cmd.load(truth, "truth")
    # Deposited files can carry extra copies in the asymmetric unit (9S60 has a
    # second A2/B2 pair). Keep only the chains that form the real complex, or
    # the unaligned extra copy renders as a stray, uncolored blob.
    cmd.remove(f"truth and not (chain {truth_target_chain}+{truth_binder_chain})")

    # Align on the TARGET chain only: the thing we want to see is whether the
    # binder lands in the same place, so the target is the fixed reference frame.
    cmd.create("pred_target", "pred and chain T")
    cmd.create("truth_target", f"truth and chain {truth_target_chain}")
    result = cmd.align("pred_target", "truth_target")
    rmsd, n_atoms = result[0], result[1]

    # Apply the same transform found above to the whole predicted object.
    cmd.matrix_copy("pred_target", "pred")

    cmd.hide("everything")
    cmd.show("cartoon", "pred or truth")
    cmd.color(target_color, f"truth and chain {truth_target_chain}")
    cmd.color("gray70", "pred and chain T")
    cmd.color(pred_binder_color, "pred and chain B")
    cmd.color(truth_binder_color, f"truth and chain {truth_binder_chain}")
    cmd.set("cartoon_transparency", 0.35, "pred and chain T")

    cmd.bg_color("white")
    cmd.orient(f"truth and chain {truth_target_chain}")
    cmd.zoom("truth or pred", buffer=5)

    cmd.set("ray_shadows", 0)
    cmd.png(f"{OUT}/{name}_overlay.png", width=1400, height=1050, dpi=150, ray=0)

    cmd.turn("y", 90)
    cmd.png(f"{OUT}/{name}_overlay_side.png", width=1400, height=1050, dpi=150, ray=0)

    cmd.save(f"{OUT}/{name}.pse")
    with open(f"{OUT}/{name}_rmsd.txt", "w") as fh:
        fh.write(f"{caption}\ntarget-chain alignment RMSD: {rmsd:.2f} A over {n_atoms} atoms\n")
    print(f"[{name}] target-chain RMSD = {rmsd:.2f} A over {n_atoms} atoms")

# B1: true binder, has its own deposited structure (9S60). Direct comparison.
render_pair(
    "B1", f"{OUT}/predicted_B1.cif", f"{OUT}/truth_9S60.cif",
    truth_binder_chain="B", truth_target_chain="A",
    pred_binder_color="marine", truth_binder_color="skyblue", target_color="wheat",
    caption="B1: DARPin 4m3 predicted vs deposited (9S60, true 65.7 nM binder)",
)

# B4: no deposited structure exists (SPR/kinetics-only negative). Aligned onto
# B1's deposited target instead, since mouse/human cathepsin B are 83% identical
# and near-superimposable -- this shows where the model PUTS the binder when
# there is nothing real to recover.
render_pair(
    "B4", f"{OUT}/predicted_B4.cif", f"{OUT}/truth_9S60.cif",
    truth_binder_chain="B", truth_target_chain="A",
    pred_binder_color="firebrick", truth_binder_color="skyblue", target_color="wheat",
    caption="B4: DARPin 4m3 + human cathepsin B (no binding observed) vs the "
            "DEPOSITED MOUSE complex (83%% identical target, aligned as a reference frame only -- "
            "human cathepsin B was never crystallized with this DARPin)",
)

# C_4CJ2: clean positive control, well-characterized, high confidence.
render_pair(
    "C_4CJ2", f"{OUT}/predicted_C_4CJ2.cif", f"{OUT}/truth_4CJ2.cif",
    truth_binder_chain="C", truth_target_chain="B",
    pred_binder_color="forest", truth_binder_color="limegreen", target_color="lightorange",
    caption="C_4CJ2: affitin H4 predicted vs deposited (4CJ2, hen egg lysozyme) -- positive control",
)

python end
