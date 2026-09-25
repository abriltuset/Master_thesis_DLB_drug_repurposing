# Preparation of PDB 6FFI for mGluR5 docking
# UCSF Chimera 1.18

close all
open 6FFI_raw.pdb

# Save crystallographic M-MPEP in the ORIGINAL PDB coordinate frame
select :D8B
write selected format pdb relative #0 #0 MMPEP_crystal.pdb
~select

# Remove T4 lysozyme fusion protein (residues 679-838, chain A)
delete :679-838.A

# Remove crystallization additives
delete :OLA
delete :MES

# Keep D8B and crystallographic waters
write format pdb relative #0 #0 6FFI_complex_clean.pdb

close all