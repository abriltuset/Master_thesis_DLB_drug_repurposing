#!/usr/bin/env python3

"""
25_finalize_pdbqt_library.py

Finalize the ligand library that will enter AutoDock Vina.

Inputs:
    approved_drugs_chembl37_prepared_final.csv
    ligands_pdbqt/

Outputs:
    approved_drugs_chembl37_vina_final.csv
    approved_drugs_chembl37_vina_excluded.csv
    approved_drugs_chembl37_vina_final_metadata.json

Only ligands with a valid generated PDBQT file are retained.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

INPUT_CSV = Path(
    "approved_drugs_chembl37_prepared_final.csv"
)

PDBQT_DIR = Path(
    "ligands_pdbqt"
)

OUTPUT_FINAL = Path(
    "approved_drugs_chembl37_vina_final.csv"
)

OUTPUT_EXCLUDED = Path(
    "approved_drugs_chembl37_vina_excluded.csv"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_vina_final_metadata.json"
)

ID_COLUMN = "parent_chembl_id"


# =============================================================================
# Read library
# =============================================================================

print("Reading prepared ligand library...")

df = pd.read_csv(INPUT_CSV)

expected_ids = set(
    df[ID_COLUMN]
    .astype(str)
    .str.strip()
)

print(f"Prepared ligands: {len(expected_ids)}")


# =============================================================================
# Read generated PDBQT files
# =============================================================================

pdbqt_files = sorted(
    PDBQT_DIR.glob("*.pdbqt")
)

pdbqt_ids = {
    file.stem
    for file in pdbqt_files
}


# =============================================================================
# Validate PDBQT files
# =============================================================================

valid_ids = set()

invalid_ids = set()


for file in pdbqt_files:

    chembl_id = file.stem

    text = file.read_text(
        encoding="utf-8",
        errors="replace"
    )

    required_markers = [
        "ROOT",
        "ENDROOT",
        "TORSDOF"
    ]

    if (
        file.stat().st_size > 0
        and all(
            marker in text
            for marker in required_markers
        )
    ):
        valid_ids.add(
            chembl_id
        )

    else:
        invalid_ids.add(
            chembl_id
        )


# =============================================================================
# Identify exclusions
# =============================================================================

missing_ids = (
    expected_ids - pdbqt_ids
)

excluded_ids = (
    missing_ids | invalid_ids
)


# =============================================================================
# Final library
# =============================================================================

final_df = (
    df[
        df[ID_COLUMN]
        .astype(str)
        .isin(valid_ids)
    ]
    .copy()
    .reset_index(drop=True)
)


final_df[
    "pdbqt_file"
] = final_df[
    ID_COLUMN
].astype(str).apply(
    lambda x:
        str(
            PDBQT_DIR
            / f"{x}.pdbqt"
        )
)


# =============================================================================
# Excluded library
# =============================================================================

excluded_df = (
    df[
        df[ID_COLUMN]
        .astype(str)
        .isin(excluded_ids)
    ]
    .copy()
    .reset_index(drop=True)
)


def exclusion_reason(chembl_id):

    chembl_id = str(
        chembl_id
    )

    if chembl_id in invalid_ids:
        return "invalid_PDBQT"

    if chembl_id in missing_ids:

        if chembl_id == "CHEMBL113178":
            return (
                "Meeko_atom_typing_failure_"
                "selenium"
            )

        return "Meeko_PDBQT_not_generated"

    return "unknown"


if len(excluded_df):

    excluded_df[
        "exclusion_reason"
    ] = excluded_df[
        ID_COLUMN
    ].apply(
        exclusion_reason
    )


# =============================================================================
# Export
# =============================================================================

final_df.to_csv(
    OUTPUT_FINAL,
    index=False
)

excluded_df.to_csv(
    OUTPUT_EXCLUDED,
    index=False
)


# =============================================================================
# Metadata
# =============================================================================

metadata = {

    "date_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "n_prepared_3d":
        int(len(df)),

    "n_valid_pdbqt":
        int(len(final_df)),

    "n_excluded_pdbqt_stage":
        int(len(excluded_df)),

    "n_missing_pdbqt":
        int(len(missing_ids)),

    "n_invalid_pdbqt":
        int(len(invalid_ids)),

    "vina_library_file":
        str(OUTPUT_FINAL),

    "pdbqt_directory":
        str(PDBQT_DIR),
}


with open(
    OUTPUT_METADATA,
    "w",
    encoding="utf-8"
) as handle:

    json.dump(
        metadata,
        handle,
        indent=4
    )


# =============================================================================
# Report
# =============================================================================

print("\n======================================")
print("FINAL VINA LIGAND LIBRARY")
print("======================================")

print(
    f"Prepared 3D ligands:        "
    f"{len(df)}"
)

print(
    f"Valid PDBQT ligands:        "
    f"{len(final_df)}"
)

print(
    f"Excluded at PDBQT stage:    "
    f"{len(excluded_df)}"
)

print(
    f"Retention:                  "
    f"{100 * len(final_df) / len(df):.2f}%"
)


if len(excluded_df):

    print("\nExcluded:")

    for _, row in excluded_df.iterrows():

        print(
            f"  {row[ID_COLUMN]}: "
            f"{row['exclusion_reason']}"
        )


print("\nFiles generated:")

print(
    f"  {OUTPUT_FINAL}"
)

print(
    f"  {OUTPUT_EXCLUDED}"
)

print(
    f"  {OUTPUT_METADATA}"
)
