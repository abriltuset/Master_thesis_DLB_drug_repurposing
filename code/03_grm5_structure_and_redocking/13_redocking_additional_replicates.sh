#!/bin/bash
set -e

# Additional redocking replicates for reproducibility.
# Seed 2026 is run separately in 10_redock_MMPEP_dry.sh.
for seed in 2027 2028 2029
do
    vina \
        --receptor 6FFI_dry.pdbqt \
        --ligand MMPEP.pdbqt \
        --config 6FFI_dry.box.txt \
        --exhaustiveness 32 \
        --num_modes 20 \
        --seed $seed \
        --out MMPEP_redocking_seed${seed}.pdbqt \
        2>&1 | tee MMPEP_redocking_seed${seed}.log

    mk_export.py \
        MMPEP_redocking_seed${seed}.pdbqt \
        -s MMPEP_redocking_seed${seed}.sdf
done
