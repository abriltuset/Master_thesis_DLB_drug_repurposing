#!/usr/bin/env python3

"""
21_prepare_ligands_ph74_3d.py

Prepare BBB-positive, docking-compatible approved drugs for Meeko.

Steps:
    1. Export ChEMBL parent SMILES to a .smi file.
    2. Protonate molecules at pH 7.4 using Molscrub.
    3. Generate 3D coordinates and explicit hydrogens.
    4. Disable tautomer enumeration, retaining the input tautomeric representation while allowing protonation-state preparation at pH 7.4.
    5. Validate the resulting SDF.

Input:
    approved_drugs_chembl37_docking_compatible.csv

Outputs:
    approved_drugs_chembl37_docking_compatible.smi
    approved_drugs_chembl37_ph74_3d.sdf
    approved_drugs_chembl37_ph74_3d.log
    approved_drugs_chembl37_ph74_3d_metadata.json
"""

import json
import shutil
import subprocess
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd
from rdkit import Chem


# =============================================================================
# Configuration
# =============================================================================

INPUT_CSV = Path(
    "approved_drugs_chembl37_docking_compatible.csv"
)

OUTPUT_SMI = Path(
    "approved_drugs_chembl37_docking_compatible.smi"
)

OUTPUT_SDF = Path(
    "approved_drugs_chembl37_ph74_3d.sdf"
)

OUTPUT_LOG = Path(
    "approved_drugs_chembl37_ph74_3d.log"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_ph74_3d_metadata.json"
)

OUTPUT_FAILED = Path(
    "approved_drugs_chembl37_ph74_failed.sdf"
)

SMILES_COLUMN = "canonical_smiles"
ID_COLUMN = "parent_chembl_id"

PH = 7.4


# =============================================================================
# Check Molscrub
# =============================================================================

scrub_executable = shutil.which("scrub.py")

if scrub_executable is None:
    raise RuntimeError(
        "scrub.py was not found. "
        "Install Molscrub before running this script."
    )


try:
    molscrub_version = version("molscrub")
except Exception:
    molscrub_version = "unknown"


print(f"Molscrub version: {molscrub_version}")


# =============================================================================
# Read library
# =============================================================================

print("\nReading docking-compatible library...")

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input file not found: {INPUT_CSV}"
    )

df = pd.read_csv(INPUT_CSV)

required = {
    SMILES_COLUMN,
    ID_COLUMN,
}

missing = required - set(df.columns)

if missing:
    raise ValueError(
        "Missing columns: "
        + ", ".join(sorted(missing))
    )


if df[SMILES_COLUMN].isna().any():
    raise ValueError(
        "Missing SMILES detected."
    )


if df[ID_COLUMN].isna().any():
    raise ValueError(
        "Missing ChEMBL IDs detected."
    )


if df[ID_COLUMN].duplicated().any():
    raise ValueError(
        "Duplicated parent ChEMBL IDs detected."
    )


print(f"Input drugs: {len(df)}")


# =============================================================================
# Write SMILES input
# =============================================================================

print("\nWriting SMILES file...")

with open(
    OUTPUT_SMI,
    "w",
    encoding="utf-8"
) as handle:

    for _, row in df.iterrows():

        smiles = str(
            row[SMILES_COLUMN]
        ).strip()

        chembl_id = str(
            row[ID_COLUMN]
        ).strip()

        handle.write(
            f"{smiles}\t{chembl_id}\n"
        )


print(f"SMILES written: {len(df)}")


# =============================================================================
# Run Molscrub
# =============================================================================

command = [
    scrub_executable,
    str(OUTPUT_SMI),
    "-o",
    str(OUTPUT_SDF),
    "--ph",
    str(PH),
    "--skip_tautomers",
    "--write_failed_mols",
    str(OUTPUT_FAILED),
    "--cpu",
    "1",
    "--etkdg_rng_seed",
    "2026",
]

print("\nRunning Molscrub...")
print(
    "Protonation: pH 7.4"
)
print(
    "Tautomer enumeration: disabled"
)


process = subprocess.run(
    command,
    capture_output=True,
    text=True
)


# Save complete log

with open(
    OUTPUT_LOG,
    "w",
    encoding="utf-8"
) as handle:

    handle.write(
        "COMMAND\n"
    )

    handle.write(
        " ".join(command)
        + "\n\n"
    )

    handle.write(
        "STDOUT\n"
    )

    handle.write(
        process.stdout
    )

    handle.write(
        "\n\nSTDERR\n"
    )

    handle.write(
        process.stderr
    )


if process.returncode != 0:

    print(
        "\nMolscrub failed."
    )

    print(
        f"Check: {OUTPUT_LOG}"
    )

    raise RuntimeError(
        f"Molscrub returned exit code "
        f"{process.returncode}"
    )


print("Molscrub completed.")


# =============================================================================
# Validate SDF
# =============================================================================

print("\nValidating generated SDF...")

supplier = Chem.SDMolSupplier(
    str(OUTPUT_SDF),
    removeHs=False
)


n_records = 0
n_valid = 0
n_invalid = 0
n_without_3d = 0


for mol in supplier:

    n_records += 1

    if mol is None:

        n_invalid += 1
        continue

    n_valid += 1

    if mol.GetNumConformers() == 0:

        n_without_3d += 1

    else:

        conformer = mol.GetConformer()

        if not conformer.Is3D():
            n_without_3d += 1


# =============================================================================
# Metadata
# =============================================================================

metadata = {

    "date_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "input_file":
        str(INPUT_CSV),

    "n_input_drugs":
        int(len(df)),

    "molscrub_version":
        molscrub_version,

    "ph":
        PH,

    "tautomer_enumeration":
        False,

    "three_dimensional_generation":
        True,

    "explicit_hydrogens":
        True,

    "n_sdf_records":
        int(n_records),

    "n_valid_sdf_records":
        int(n_valid),

    "n_invalid_sdf_records":
        int(n_invalid),

    "n_without_3d_coordinates":
        int(n_without_3d),

    "outputs": {
        "smiles":
            str(OUTPUT_SMI),

        "sdf":
            str(OUTPUT_SDF),

        "log":
            str(OUTPUT_LOG)
    }
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
print("LIGAND 3D PREPARATION")
print("======================================")

print(
    f"Input approved drugs:       "
    f"{len(df)}"
)

print(
    f"SDF molecular states:       "
    f"{n_records}"
)

print(
    f"Valid SDF records:          "
    f"{n_valid}"
)

print(
    f"Invalid SDF records:        "
    f"{n_invalid}"
)

print(
    f"Without 3D coordinates:     "
    f"{n_without_3d}"
)

print(
    f"Protonation pH:             "
    f"{PH}"
)

print(
    "Tautomer enumeration:       disabled"
)

print("\nFiles generated:")

print(f"  {OUTPUT_SMI}")
print(f"  {OUTPUT_SDF}")
print(f"  {OUTPUT_LOG}")
print(f"  {OUTPUT_METADATA}")
