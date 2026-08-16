# References for the go/no-go validation set

Every measurement in `cases.yaml` is line-cited here. PDB fields were verified
directly against the PDB entry table, not taken from these papers.

## What each reference supplies

**[1]** supplies all of Tier B. DARPin 4m3 binds mouse cathepsin B at
K_D = 65.7 nM (pH 6) and 108.3 nM (pH 7) by SPR, and inhibits with a competitive
K_i of 26.7 nM. Against human cathepsin B it does nothing — no SEC complex, no
SPR response, no inhibition "even at 1000-fold molar excess", despite 83%
sequence identity and 0.57–0.64 Å Cα RMSD between the two orthologs. The paper
attributes selectivity to I65/Q66 in mouse against S65/M66 in human. It also
hands us two ready-made negative controls: the mouse zymogen (no complex) and
DARPin E3_5, an anti-MBP binder that produces no thermal stabilisation of
cathepsin B in cell lysate.

**[2]** supplies Tier A and its baselines. The dArmRP–peptide system was chosen
there for the same reason it is useful here: single-residue peptide
substitutions span K_D values from 1 to 1000 nM across five characterised
binding pockets, which is a much finer discrimination task than binder /
non-binder. It also reports what three established physics-based methods
(flex ddG, BBK*, PocketOptimizer) achieve on exactly that data, so Tier A comes
with a published bar rather than an invented one — and it shows predictions
shifting with the choice of input crystal structure (6SA8 against 5AEI), which
is worth reproducing as a robustness check.

--------
REFERENCES
[1] Zarić, M. et al. "Structural and Proteomic Analysis of the Mouse Cathepsin
    B-DARPin 4m3 Complex Reveals Species-Specific Binding Determinants."
    *International Journal of Molecular Sciences* (2025). doi:10.3390/ijms262411910
    https://paperclip.gxl.ai/citations/papers/PMC12732864#L10,L19,L21,L22,L24,L25,L27
[2] Ayyildiz, M., Noske, J., Gisdon, F. J., Kynast, J. P. & Höcker, B.
    "Evaluation of Physics-Based Protein Design Methods for Predicting Single
    Residue Effects on Peptide Binding Specificities." *Journal of Computational
    Chemistry* (2025). doi:10.1002/jcc.70160
    https://paperclip.gxl.ai/citations/papers/PMC12202725#L10,L17,L19,L20,L24
