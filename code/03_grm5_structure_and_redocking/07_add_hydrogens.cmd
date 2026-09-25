# Add hydrogens to 6FFI receptor variants
# UCSF Chimera 1.18
# No energy minimization performed

# ============================================================
# Dry receptor
# ============================================================

close all
open 6FFI_receptor_dry.pdb

addh hbond true

write format pdb relative #0 #0 6FFI_receptor_dry_H.pdb


# ============================================================
# Receptor with structural water HOH 4133
# ============================================================

close all
open 6FFI_receptor_w4133.pdb

addh hbond true

write format pdb relative #0 #0 6FFI_receptor_w4133_H.pdb

close all