#!/usr/bin/env python3

"""
22_recover_failed_ligands_3d.py

Fallback 3D preparation for compounds that failed the standard
Molscrub ETKDG workflow.

Strategy:
    1. Read molecules written by Molscrub as failed.
    2. Apply the same Molscrub protonation strategy at pH 7.4,
       without generating 3D coordinates.
    3. Add explicit hydrogens.
    4. Generate 3D coordinates with RDKit ETKDGv3 using
       random-coordinate initialization and a fixed RNG seed.
    5. Minimize geometry using MMFF94s.
    6. Merge recovered molecules with the main successful SDF.

Inputs:
    approved_drugs_chembl37_ph74_failed.sdf
    approved_drugs_chembl37_ph74_3d.sdf

Outputs:
    approved_drugs_chembl37_ph74_3d_recovered.sdf
    approved_drugs_chembl37_ph74_3d_final.sdf
"""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from molscrub import Scrub


FAILED_INPUT = Path(
    "approved_drugs_chembl37_ph74_failed.sdf"
)

MAIN_INPUT = Path(
    "approved_drugs_chembl37_ph74_3d.sdf"
)

RECOVERED_OUTPUT = Path(
    "approved_drugs_chembl37_ph74_3d_recovered.sdf"
)

FINAL_OUTPUT = Path(
    "approved_drugs_chembl37_ph74_3d_final.sdf"
)

PH = 7.4
RNG_SEED = 2026


# =============================================================================
# Molscrub protonation only
# =============================================================================

scrub = Scrub(
    ph_low=PH,
    ph_high=PH,
    skip_tautomers=True,
    skip_gen3d=True,
)


print("Reading failed molecules...")

supplier = Chem.SDMolSupplier(
    str(FAILED_INPUT),
    removeHs=False
)

failed_mols = [
    mol for mol in supplier
    if mol is not None
]

print(f"Failed molecules to recover: {len(failed_mols)}")


# =============================================================================
# Recover
# =============================================================================

recovered = []
still_failed = []


for mol in failed_mols:

    name = (
        mol.GetProp("_Name")
        if mol.HasProp("_Name")
        else "UNKNOWN"
    )

    print(f"\nRecovering {name}...")

    try:

        # Same protonation framework as main Molscrub run,
        # but without the failing 3D generation step.
        states = scrub(mol)

        if len(states) != 1:
            print(
                f"  WARNING: generated {len(states)} "
                f"protonation states."
            )

        success_for_molecule = False

        for state in states:

            state = Chem.AddHs(
                state,
                addCoords=False
            )

            params = AllChem.ETKDGv3()

            params.randomSeed = RNG_SEED
            params.useRandomCoords = True
            params.enforceChirality = True

            conf_id = AllChem.EmbedMolecule(
                state,
                params
            )

            if conf_id < 0:
                continue

            # Match the main Molscrub force field as closely as possible.
            if AllChem.MMFFHasAllMoleculeParams(state):

                AllChem.MMFFOptimizeMolecule(
                    state,
                    mmffVariant="MMFF94s",
                    maxIters=400
                )

            else:

                # Conservative fallback if MMFF parameters
                # are unavailable.
                AllChem.UFFOptimizeMolecule(
                    state,
                    maxIters=400
                )

            state.SetProp(
                "_Name",
                name
            )

            state.SetProp(
                "3D_preparation_method",
                "ETKDGv3_random_coords_fallback"
            )

            state.SetProp(
                "ETKDG_random_seed",
                str(RNG_SEED)
            )

            recovered.append(
                state
            )

            success_for_molecule = True

        if not success_for_molecule:
            still_failed.append(name)

    except Exception as exc:

        print(
            f"  FAILED: {exc}"
        )

        still_failed.append(name)


# =============================================================================
# Write recovered compounds
# =============================================================================

writer = Chem.SDWriter(
    str(RECOVERED_OUTPUT)
)

for mol in recovered:
    writer.write(mol)

writer.close()


# =============================================================================
# Merge main + recovered
# =============================================================================

main_supplier = Chem.SDMolSupplier(
    str(MAIN_INPUT),
    removeHs=False
)

main_mols = [
    mol for mol in main_supplier
    if mol is not None
]


writer = Chem.SDWriter(
    str(FINAL_OUTPUT)
)

for mol in main_mols:
    writer.write(mol)

for mol in recovered:
    writer.write(mol)

writer.close()


# =============================================================================
# Report
# =============================================================================

print("\n======================================")
print("FAILED-LIGAND 3D RECOVERY")
print("======================================")

print(
    f"Failed inputs:              "
    f"{len(failed_mols)}"
)

print(
    f"Recovered states:           "
    f"{len(recovered)}"
)

print(
    f"Still failed:               "
    f"{len(still_failed)}"
)

print(
    f"Main successful structures: "
    f"{len(main_mols)}"
)

print(
    f"Final SDF records:          "
    f"{len(main_mols) + len(recovered)}"
)


if still_failed:

    print("\nStill failed:")

    for name in still_failed:
        print(f"  {name}")


print("\nFiles generated:")
print(f"  {RECOVERED_OUTPUT}")
print(f"  {FINAL_OUTPUT}")
