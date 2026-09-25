# Restore native Cys634 in PDB 6FFI
# UCSF Chimera 1.18
#
# YCM634 is a carbamidomethylated cysteine.
# Native N, CA, CB, SG, C and O coordinates are retained.
# Only the artificial carbamidomethyl substituent is removed.

close all
open 6FFI_complex_clean.pdb

# Remove atoms belonging to the artificial YCM modification
delete :634.A@CD,CE,OZ1,NZ2

# Rename modified residue YCM634 as native cysteine
setattr r type CYS :634.A

write format pdb relative #0 #0 6FFI_complex_clean_C634.pdb

close all