#!/bin/bash
set -e

# Define the mGluR5 allosteric docking box
# Box based on crystallographic M-MPEP coordinates from PDB 6FFI
# 5 A padding around the experimental ligand

mk_prepare_receptor.py \
    --read_pdb 6FFI_receptor_dry_H.pdb \
    -o 6FFI_dry \
    -p \
    -v \
    --box_enveloping 6ffi_J_D8B.sdf \
    --padding 5
