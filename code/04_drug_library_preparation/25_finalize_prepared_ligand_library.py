#!/usr/bin/env python3

"""
23_finalize_prepared_ligand_library.py

Finalize the ligand library after reproducible 3D preparation.

Inputs:
    approved_drugs_chembl37_docking_compatible.csv
    approved_drugs_chembl37_ph74_3d_final.sdf

Outputs:
    approved_drugs_chembl37_prepared_final.csv
    approved_drugs_chembl37_3d_failed_excluded.csv

Compounds absent from the final 3D-prepared SDF are considered
reproducible 3D-generation failures and are excluded from downstream docking.
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem


INPUT_CSV = Path(
    "approved_drugs_chembl37_docking_compatible.csv"
)

FINAL_SDF = Path(
    "approved_drugs_chembl37_ph74_3d_final.sdf"
)

OUTPUT_CSV = Path(
    "approved_drugs_chembl37_prepared_final.csv"
)

OUTPUT_EXCLUDED = Path(
    "approved_drugs_chembl37_3d_failed_excluded.csv"
)


# =============================================================================
# Read original library
# =============================================================================

df = pd.read_csv(INPUT_CSV)

print(f"Input compatible drugs: {len(df)}")


# =============================================================================
# Read final prepared SDF IDs
# =============================================================================

final_supplier = Chem.SDMolSupplier(
    str(FINAL_SDF),
    removeHs=False
)

final_ids = []

for mol in final_supplier:

    if mol is None:
        continue

    if not mol.HasProp("_Name"):
        raise ValueError(
            "Final SDF molecule without _Name."
        )

    final_ids.append(
        mol.GetProp("_Name").strip()
    )


final_ids = set(final_ids)

input_ids = set(
    df["parent_chembl_id"]
    .astype(str)
)


# =============================================================================
# Validate
# =============================================================================

unexpected_ids = (
    final_ids - input_ids
)

if unexpected_ids:

    raise ValueError(
        "Unexpected compounds in final SDF: "
        + ", ".join(sorted(unexpected_ids))
    )


failed_ids = (
    input_ids - final_ids
)


# =============================================================================
# Final and excluded tables
# =============================================================================

final_df = (
    df[
        df["parent_chembl_id"]
        .astype(str)
        .isin(final_ids)
    ]
    .copy()
    .reset_index(drop=True)
)


excluded_df = (
    df[
        df["parent_chembl_id"]
        .astype(str)
        .isin(failed_ids)
    ]
    .copy()
    .reset_index(drop=True)
)


excluded_df[
    "exclusion_reason"
] = (
    "reproducible_3D_generation_failure"
)


final_df.to_csv(
    OUTPUT_CSV,
    index=False
)

excluded_df.to_csv(
    OUTPUT_EXCLUDED,
    index=False
)


# =============================================================================
# Report
# =============================================================================

print("\n======================================")
print("FINAL PREPARED LIGAND LIBRARY")
print("======================================")

print(
    f"Input docking-compatible:   "
    f"{len(df)}"
)

print(
    f"Successfully prepared:      "
    f"{len(final_df)}"
)

print(
    f"Excluded after 3D failure:  "
    f"{len(excluded_df)}"
)

print(
    f"Retention:                  "
    f"{100 * len(final_df) / len(df):.2f}%"
)


if len(excluded_df):

    print("\nExcluded IDs:")

    for chembl_id in excluded_df[
        "parent_chembl_id"
    ]:
        print(f"  {chembl_id}")


print("\nFiles generated:")
print(f"  {OUTPUT_CSV}")
print(f"  {OUTPUT_EXCLUDED}")