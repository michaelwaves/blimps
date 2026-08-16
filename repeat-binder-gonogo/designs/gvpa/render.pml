# Headless renders of the GvpA patch, the GvpC placement, and the curvature result.
#   pymol -cq render.pml
# PyMOL from Homebrew (/opt/homebrew/bin/pymol). The pip wheel is unusable here:
# it links libpng from the packager's own mamba prefix and forces numpy>=2,
# which breaks DockQ.

set ray_opaque_background, 1
bg_color white
set ray_shadows, 0
set orthoscopic, 1
set antialias, 2
set cartoon_transparency, 0

# ---------------------------------------------------------------- full patch
load ../../GvpA_surface_patch_5ribs_x11.pdb, patch
hide everything
show cartoon, patch
spectrum chain, rainbow, patch
set cartoon_transparency, 0.1, patch
orient patch
turn x, 90
png patch_axial.png, width=1600, height=1200, dpi=200, ray=1
turn x, -90
png patch_side.png, width=1600, height=1200, dpi=200, ray=1

# ------------------------------------------------- target + GvpC + footprint
delete patch
load gvpa_target.pdb, gvpa
load gvpc_reference.pdb, gvpc
hide everything
show cartoon, gvpa
color grey80, gvpa
color palecyan,  gvpa and chain A
color lightblue, gvpa and chain B
color skyblue,   gvpa and chain C
color lightblue, gvpa and chain D
color palecyan,  gvpa and chain E
set cartoon_transparency, 0.3, gvpa

select footprint, gvpa and ((chain A and resi 55+59+64) or (chain B and resi 52+55+56+59+60+64+65+66) or (chain C and resi 55+56+57+59+60+64+65+66) or (chain D and resi 56+59+60+64+65+66) or (chain E and resi 59+60+66))
show sticks, footprint and sidechain
color orange, footprint
set stick_radius, 0.2

show cartoon, gvpc
cartoon tube
set cartoon_tube_radius, 1.2, gvpc
color firebrick, gvpc

# The three sub-5A contacts. Only three -- that sparseness is the finding,
# not a rendering artefact, so it is drawn rather than asserted.
distance close_contacts, gvpc, gvpa, 5.0, mode=2
color red, close_contacts
set dash_width, 4

deselect
orient gvpa
png target_gvpc.png, width=1600, height=1200, dpi=200, ray=1
turn y, 90
png target_gvpc_rot90.png, width=1600, height=1200, dpi=200, ray=1

print "rendered: patch_axial, patch_side, target_gvpc, target_gvpc_rot90"
