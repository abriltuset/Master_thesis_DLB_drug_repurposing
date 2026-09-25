#!/bin/bash
set -e

# Redocking validation of crystallographic M-MPEP in PDB 6FFI
# AutoDock Vina
# Receptor: dry 6FFI
# Fixed seed for reproducibility

vina \
    --receptor 6FFI_dry.pdbqt \
    --ligand MMPEP.pdbqt \
    --config 6FFI_dry.box.txt \
    --exhaustiveness 32 \
    --num_modes 20 \
    --seed 2026 \
    --out MMPEP_redocking_dry.pdbqt \
    2>&1 | tee MMPEP_redocking_dry.log

# Export docked poses to SDF for RMSD analysis
mk_export.py \
    MMPEP_redocking_dry.pdbqt \
    -s MMPEP_redocking_dry.sdf