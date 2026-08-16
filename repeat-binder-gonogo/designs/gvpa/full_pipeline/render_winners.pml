set ray_opaque_background, off
set antialias, 2
set cartoon_fancy_helices, on
set depth_cue, 0
bg_color white

load designs/gvpa/full_pipeline/rfd3_inputs/c02.cif, constrained
hide everything, constrained
show surface, constrained and chain B+C+D+E+F
color gray80, constrained and chain B+C+D+E+F
set transparency, 0.32, constrained and chain B+C+D+E+F
show cartoon, constrained and chain A
color cyan, constrained and chain A
select constrained_hotspots, constrained and ((chain B and resi 55+59+64) or (chain C and resi 52+55+56+59+60+64+65+66) or (chain D and resi 55+56+57+59+60+64+65+66) or (chain E and resi 56+59+60+64+65+66) or (chain F and resi 59+60+66))
show spheres, constrained_hotspots and name CA
set sphere_scale, 0.55, constrained_hotspots and name CA
color red, constrained_hotspots
orient constrained
zoom constrained, 4
ray 1600, 1000
png designs/gvpa/full_pipeline/constrained_winner.png, dpi=180
disable constrained

load designs/gvpa/full_pipeline/rfd3_inputs/f02.cif, free
hide everything, free
show surface, free and chain B+C+D+E+F
color gray80, free and chain B+C+D+E+F
set transparency, 0.32, free and chain B+C+D+E+F
show cartoon, free and chain A
color orange, free and chain A
select free_hotspots, free and ((chain B and resi 55+59+64) or (chain C and resi 52+55+56+59+60+64+65+66) or (chain D and resi 55+56+57+59+60+64+65+66) or (chain E and resi 56+59+60+64+65+66) or (chain F and resi 59+60+66))
show spheres, free_hotspots and name CA
set sphere_scale, 0.55, free_hotspots and name CA
color red, free_hotspots
orient free
zoom free, 4
ray 1600, 1000
png designs/gvpa/full_pipeline/free_winner.png, dpi=180
quit
