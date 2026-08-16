# GvpA shell patch, the GvpC placement, and the designed-binder site.
#
#   pymol view.pml
#
# PyMOL could not be installed in the project venv on this machine: the
# pymol-open-source PyPI wheel links against libpng from the packager's own
# mamba prefix, and installing it also forces numpy>=2, which breaks DockQ.
# Use a conda/homebrew PyMOL and run this script against the files beside it.

load gvpa_target.pdb, gvpa
load gvpc_reference.pdb, gvpc
load design_complex.pdb, complex

bg_color white
hide everything
set ray_shadows, 0
set cartoon_transparency, 0.15
set orthoscopic, 1

# --- GvpA shell surface: the five subunits GvpC lies across -------------
show cartoon, gvpa
color grey70, gvpa
set cartoon_transparency, 0.35, gvpa

# one subunit per colour, so the rib periodicity is visible
color palecyan,   gvpa and chain A
color lightblue,  gvpa and chain B
color skyblue,    gvpa and chain C
color lightblue,  gvpa and chain D
color palecyan,   gvpa and chain E

# --- the GvpC footprint: what the natural binder touches ----------------
# 28 residues within 8 A of the placed GvpC helix. NOTE the cutoff: at the
# conventional 5 A only THREE residues qualify, because the 8GBS GvpC is an
# integrative placement rather than a resolved interface.
select footprint, gvpa and (chain A and resi 55+59+64) \
                        or (chain B and resi 52+55+56+59+60+64+65+66) \
                        or (chain C and resi 55+56+57+59+60+64+65+66) \
                        or (chain D and resi 56+59+60+64+65+66) \
                        or (chain E and resi 59+60+66)
show sticks, footprint and not (name C+N+O)
color orange, footprint
set stick_radius, 0.18, footprint

# --- GvpC as modelled: backbone only, no side chains, no sequence -------
show cartoon, gvpc
color firebrick, gvpc
set cartoon_transparency, 0.0, gvpc
# drawn as a ribbon deliberately: there are no side chains to show.
cartoon tube, gvpc

# --- the gap that matters ------------------------------------------------
# Nearest approach is 3.75 A and only 3 residues fall inside 5 A. Distances
# make that legible rather than leaving it to a caption.
distance contacts, gvpc, gvpa, 5.0, mode=2
color red, contacts
show dashes, contacts

deselect
orient gvpc
zoom gvpc, 12

set_view_help = "GvpC (dark red tube) lies across five GvpA subunits. Orange \
sticks are the 8 A footprint. Red dashes are the only sub-5 A contacts -- there \
are three. That sparseness is why ProteinMPNN, given this backbone as a \
scaffold, returned low-complexity exposed-helix sequences with no positional \
signal above a composition-matched shuffle."

print "Loaded. GvpC = dark red tube (backbone only, poly-UNK)."
print "Orange = 8 A footprint on GvpA. Red dashes = the 3 contacts under 5 A."
print "See site.json for the footprint residue list and the caveat."
