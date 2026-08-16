set ray_opaque_background, on
set antialias, 2
set cartoon_fancy_helices, on
set depth_cue, 0
set orthoscopic, on
set ray_shadows, 0
set specular, 0.18
bg_color white
viewport 960, 600

# The RFdiffusion3 complex is recentered relative to the original GvpA patch.
# Build a single reference object so the alignment moves the GvpC trace and the
# patch together, then hide the reference patch after alignment.
load designs/gvpa/gvpa_target.pdb, reference_patch
load designs/gvpa/gvpc_reference.pdb, reference_gvpc
alter reference_gvpc, chain='G'
create reference_complex, reference_patch or reference_gvpc
delete reference_patch
delete reference_gvpc

load designs/gvpa/full_pipeline/rfd3_inputs/c02.cif, constrained
align reference_complex and chain A+B+C+D+E and name CA, constrained and chain B+C+D+E+F and name CA, cycles=0
hide everything
show surface, constrained and chain B+C+D+E+F
color gray80, constrained and chain B+C+D+E+F
set transparency, 0.26, constrained and chain B+C+D+E+F
show cartoon, constrained and chain A
color cyan, constrained and chain A
show cartoon, reference_complex and chain G
color orange, reference_complex and chain G
set cartoon_transparency, 0.30, reference_complex and chain G
select constrained_hotspots, constrained and ((chain B and resi 55+59+64) or (chain C and resi 52+55+56+59+60+64+65+66) or (chain D and resi 55+56+57+59+60+64+65+66) or (chain E and resi 56+59+60+64+65+66) or (chain F and resi 59+60+66))
show spheres, constrained_hotspots and name CA
set sphere_scale, 0.55, constrained_hotspots and name CA
color magenta, constrained_hotspots
orient constrained
zoom constrained, 5
mset 1 x72
util.mroll 1,72,1
set ray_trace_frames, 1
mpng /private/tmp/gvpa_rotate_constrained/frame_

delete constrained
delete constrained_hotspots
load designs/gvpa/full_pipeline/rfd3_inputs/f02.cif, free
align reference_complex and chain A+B+C+D+E and name CA, free and chain B+C+D+E+F and name CA, cycles=0
hide everything
show surface, free and chain B+C+D+E+F
color gray80, free and chain B+C+D+E+F
set transparency, 0.26, free and chain B+C+D+E+F
show cartoon, free and chain A
color cyan, free and chain A
show cartoon, reference_complex and chain G
color orange, reference_complex and chain G
set cartoon_transparency, 0.30, reference_complex and chain G
select free_hotspots, free and ((chain B and resi 55+59+64) or (chain C and resi 52+55+56+59+60+64+65+66) or (chain D and resi 55+56+57+59+60+64+65+66) or (chain E and resi 56+59+60+64+65+66) or (chain F and resi 59+60+66))
show spheres, free_hotspots and name CA
set sphere_scale, 0.55, free_hotspots and name CA
color magenta, free_hotspots
orient free
zoom free, 5
mset 1 x72
util.mroll 1,72,1
set ray_trace_frames, 1
mpng /private/tmp/gvpa_rotate_free/frame_
quit
