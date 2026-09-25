# Generate receptor variants from corrected 6FFI structure
# UCSF Chimera 1.18

# ============================================================
# Dry receptor
# ============================================================

close all
open 6FFI_complex_clean_C634_fixed.pdb

# Remove crystallographic reference ligand
delete :D8B

# Remove all crystallographic waters
delete solvent

write format pdb relative #0 #0 6FFI_receptor_dry.pdb


# ============================================================
# Receptor retaining structural water HOH 4133
# ============================================================

close all
open 6FFI_complex_clean_C634_fixed.pdb

# Remove crystallographic reference ligand
delete :D8B

# Remove all solvent except structural water HOH 4133
select solvent & ~:4133.A
delete sel

write format pdb relative #0 #0 6FFI_receptor_w4133.pdb

close all