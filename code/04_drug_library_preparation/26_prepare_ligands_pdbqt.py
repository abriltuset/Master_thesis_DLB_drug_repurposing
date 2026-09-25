#!/usr/bin/env python3

"""
24_prepare_ligands_pdbqt.py

Prepare the final approved-drug ligand library in PDBQT format
using Meeko for AutoDock Vina screening.

Inputs:
    approved_drugs_chembl37_ph74_3d_final.sdf
    approved_drugs_chembl37_prepared_final.csv

Outputs:
    ligands_pdbqt/
        CHEMBLxxxx.pdbqt
        ...
    approved_drugs_chembl37_pdbqt_failures.csv
    approved_drugs_chembl37_pdbqt_metadata.json
    approved_drugs_chembl37_pdbqt_preparation.log

Preparation:
    - Meeko mk_prepare_ligand.py
    - Gasteiger partial charges
    - default Meeko ligand flexibility treatment
"""

import json
import shutil
import subprocess
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

INPUT_SDF = Path(
    "approved_drugs_chembl37_ph74_3d_final.sdf"
)

INPUT_CSV = Path(
    "approved_drugs_chembl37_prepared_final.csv"
)

OUTPUT_DIR = Path(
    "ligands_pdbqt"
)

OUTPUT_FAILURES = Path(
    "approved_drugs_chembl37_pdbqt_failures.csv"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_pdbqt_metadata.json"
)

OUTPUT_LOG = Path(
    "approved_drugs_chembl37_pdbqt_preparation.log"
)

ID_COLUMN = "parent_chembl_id"


# =============================================================================
# Check inputs
# =============================================================================

print("Reading final prepared ligand library...")

if not INPUT_SDF.exists():
    raise FileNotFoundError(
        f"Input SDF not found: {INPUT_SDF}"
    )

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input CSV not found: {INPUT_CSV}"
    )

df = pd.read_csv(INPUT_CSV)

if ID_COLUMN not in df.columns:
    raise ValueError(
        f"Column '{ID_COLUMN}' not found."
    )

expected_ids = set(
    df[ID_COLUMN]
    .astype(str)
    .str.strip()
)

print(f"Expected ligands: {len(expected_ids)}")


# =============================================================================
# Check Meeko
# =============================================================================

meeko_executable = shutil.which(
    "mk_prepare_ligand.py"
)

if meeko_executable is None:
    raise RuntimeError(
        "mk_prepare_ligand.py not found "
        "in the active environment."
    )

try:
    meeko_version = version("meeko")
except Exception:
    meeko_version = "unknown"


print(f"Meeko version: {meeko_version}")


# =============================================================================
# Prepare clean output directory
# =============================================================================

if OUTPUT_DIR.exists():

    print(
        f"\nRemoving previous output directory: "
        f"{OUTPUT_DIR}"
    )

    shutil.rmtree(
        OUTPUT_DIR
    )

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# Run Meeko
# =============================================================================

command = [
    meeko_executable,
    "-i",
    str(INPUT_SDF),
    "--multimol_outdir",
    str(OUTPUT_DIR),
    "--charge_model",
    "gasteiger",
]


print("\nRunning Meeko ligand preparation...")

print(
    "Charge model: Gasteiger"
)


process = subprocess.run(
    command,
    capture_output=True,
    text=True
)


# Save full log

with open(
    OUTPUT_LOG,
    "w",
    encoding="utf-8"
) as handle:

    handle.write("COMMAND\n")
    handle.write(
        " ".join(command)
        + "\n\n"
    )

    handle.write("STDOUT\n")
    handle.write(
        process.stdout
    )

    handle.write(
        "\n\nSTDERR\n"
    )
    handle.write(
        process.stderr
    )


# Exit code 4 = partial Meeko failure.
# We still continue so the exact missing ligands
# can be identified and documented.

if process.returncode not in (0, 4):

    raise RuntimeError(
        f"Meeko returned unexpected exit code "
        f"{process.returncode}. "
        f"Check {OUTPUT_LOG}"
    )


# =============================================================================
# Audit generated PDBQT files
# =============================================================================

print("\nAuditing generated PDBQT files...")


pdbqt_files = sorted(
    OUTPUT_DIR.glob("*.pdbqt")
)


generated_ids = set(
    file.stem
    for file in pdbqt_files
)


missing_ids = sorted(
    expected_ids - generated_ids
)

unexpected_ids = sorted(
    generated_ids - expected_ids
)


empty_files = [
    file.name
    for file in pdbqt_files
    if file.stat().st_size == 0
]


# =============================================================================
# Basic PDBQT validation
# =============================================================================

invalid_pdbqt = []


for file in pdbqt_files:

    text = file.read_text(
        encoding="utf-8",
        errors="replace"
    )

    required_markers = [
        "ROOT",
        "ENDROOT",
        "TORSDOF",
    ]

    if not all(
        marker in text
        for marker in required_markers
    ):
        invalid_pdbqt.append(
            file.name
        )


# =============================================================================
# Failure table
# =============================================================================

failure_rows = []


for chembl_id in missing_ids:

    original = df[
        df[ID_COLUMN]
        .astype(str)
        .eq(chembl_id)
    ]

    name = None

    if (
        len(original) == 1
        and "name" in original.columns
    ):
        name = original.iloc[0]["name"]

    failure_rows.append({
        "parent_chembl_id":
            chembl_id,

        "name":
            name,

        "failure_stage":
            "Meeko_PDBQT_preparation",

        "reason":
            "no_PDBQT_generated",
    })


failure_df = pd.DataFrame(
    failure_rows
)


failure_df.to_csv(
    OUTPUT_FAILURES,
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

    "input_sdf":
        str(INPUT_SDF),

    "input_csv":
        str(INPUT_CSV),

    "meeko_version":
        meeko_version,

    "charge_model":
        "gasteiger",

    "n_expected_ligands":
        int(len(expected_ids)),

    "n_pdbqt_files":
        int(len(pdbqt_files)),

    "n_missing_ligands":
        int(len(missing_ids)),

    "n_unexpected_files":
        int(len(unexpected_ids)),

    "n_empty_pdbqt":
        int(len(empty_files)),

    "n_invalid_pdbqt":
        int(len(invalid_pdbqt)),

    "meeko_exit_code":
        int(process.returncode),

    "output_directory":
        str(OUTPUT_DIR),

    "log_file":
        str(OUTPUT_LOG),
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
print("MEEKO PDBQT PREPARATION")
print("======================================")

print(
    f"Expected ligands:           "
    f"{len(expected_ids)}"
)

print(
    f"PDBQT files generated:      "
    f"{len(pdbqt_files)}"
)

print(
    f"Missing ligands:            "
    f"{len(missing_ids)}"
)

print(
    f"Unexpected files:           "
    f"{len(unexpected_ids)}"
)

print(
    f"Empty PDBQT files:          "
    f"{len(empty_files)}"
)

print(
    f"Invalid PDBQT files:        "
    f"{len(invalid_pdbqt)}"
)

print(
    f"Meeko exit code:            "
    f"{process.returncode}"
)


if missing_ids:

    print("\nMissing ligand IDs:")

    for chembl_id in missing_ids:
        print(
            f"  {chembl_id}"
        )


if invalid_pdbqt:

    print("\nInvalid PDBQT files:")

    for filename in invalid_pdbqt:
        print(
            f"  {filename}"
        )


print("\nFiles generated:")

print(
    f"  {OUTPUT_DIR}/"
)

print(
    f"  {OUTPUT_FAILURES}"
)

print(
    f"  {OUTPUT_METADATA}"
)

print(
    f"  {OUTPUT_LOG}"
)
