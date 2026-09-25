#!/usr/bin/env python3

"""
26_vina_virtual_screening.py

Virtual screening of the final approved-drug library against
the validated mGluR5 NAM binding site in PDB 6FFI.

Validated docking protocol:
    receptor       = 6FFI_dry.pdbqt
    center         = (-23.720, -4.974, 41.888)
    size           = (13.367, 18.945, 18.362) Angstrom
    exhaustiveness = 32
    num_modes      = 20
    energy_range   = 3 kcal/mol
    spacing        = 0.375 Angstrom

Primary screening:
    fixed random seed = 2026

The script supports:
    --limit N       Run only the first N ligands (pilot mode)
    --output-dir    Select output directory
    --resume        Skip ligands already successfully docked

Input:
    approved_drugs_chembl37_vina_final.csv

Outputs:
    <output-dir>/poses/
    <output-dir>/logs/
    <output-dir>/vina_screening_results.csv
    <output-dir>/vina_screening_metadata.json
"""

import argparse
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# =============================================================================
# Validated protocol
# =============================================================================

INPUT_CSV = Path(
    "approved_drugs_chembl37_vina_final.csv"
)

RECEPTOR = Path(
    "6FFI_dry.pdbqt"
)

CENTER_X = -23.720
CENTER_Y = -4.974
CENTER_Z = 41.888

SIZE_X = 13.367
SIZE_Y = 18.945
SIZE_Z = 18.362

EXHAUSTIVENESS = 32
NUM_MODES = 20
ENERGY_RANGE = 3
SPACING = 0.375

SEED = 2026


# =============================================================================
# Arguments
# =============================================================================

parser = argparse.ArgumentParser()

parser.add_argument(
    "--limit",
    type=int,
    default=None,
    help="Only dock the first N ligands."
)

parser.add_argument(
    "--output-dir",
    default="vina_screening",
    help="Output directory."
)

parser.add_argument(
    "--resume",
    action="store_true",
    help="Skip ligands with existing valid output."
)

args = parser.parse_args()


OUTPUT_DIR = Path(
    args.output_dir
)

POSE_DIR = (
    OUTPUT_DIR / "poses"
)

LOG_DIR = (
    OUTPUT_DIR / "logs"
)

RESULTS_CSV = (
    OUTPUT_DIR
    / "vina_screening_results.csv"
)

METADATA_JSON = (
    OUTPUT_DIR
    / "vina_screening_metadata.json"
)


# =============================================================================
# Helpers
# =============================================================================

def get_vina_version(vina_executable):

    result = subprocess.run(
        [vina_executable, "--version"],
        capture_output=True,
        text=True
    )

    text = (
        result.stdout.strip()
        or result.stderr.strip()
    )

    return text


def parse_best_affinity(pdbqt_file):

    """
    Read the best Vina score from:
        REMARK VINA RESULT:
    """

    if not pdbqt_file.exists():
        return None

    with open(
        pdbqt_file,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as handle:

        for line in handle:

            if line.startswith(
                "REMARK VINA RESULT:"
            ):

                fields = line.split()

                try:
                    return float(
                        fields[3]
                    )

                except (
                    IndexError,
                    ValueError
                ):
                    return None

    return None


# =============================================================================
# Check Vina
# =============================================================================

vina = shutil.which("vina")

if vina is None:
    raise RuntimeError(
        "vina executable not found "
        "in the active environment."
    )


vina_version = get_vina_version(
    vina
)


print(f"Vina: {vina_version}")


# =============================================================================
# Check inputs
# =============================================================================

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Library not found: {INPUT_CSV}"
    )


if not RECEPTOR.exists():
    raise FileNotFoundError(
        f"Receptor not found: {RECEPTOR}"
    )


df = pd.read_csv(
    INPUT_CSV
)


required_columns = {
    "parent_chembl_id",
    "pdbqt_file",
}


missing = (
    required_columns
    - set(df.columns)
)

if missing:
    raise ValueError(
        "Missing required columns: "
        + ", ".join(sorted(missing))
    )


if args.limit is not None:

    if args.limit <= 0:
        raise ValueError(
            "--limit must be > 0"
        )

    df = df.head(
        args.limit
    ).copy()


print(
    f"Ligands selected: {len(df)}"
)


# =============================================================================
# Output directories
# =============================================================================

POSE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

LOG_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# Docking loop
# =============================================================================

results = []


print("\nStarting Vina screening...\n")


for number, (_, row) in enumerate(
    df.iterrows(),
    start=1
):

    chembl_id = str(
        row["parent_chembl_id"]
    )

    ligand = Path(
        row["pdbqt_file"]
    )

    output_pose = (
        POSE_DIR
        / f"{chembl_id}_out.pdbqt"
    )

    output_log = (
        LOG_DIR
        / f"{chembl_id}.log"
    )


    # -------------------------------------------------------------------------
    # Resume
    # -------------------------------------------------------------------------

    if (
        args.resume
        and output_pose.exists()
    ):

        previous_score = (
            parse_best_affinity(
                output_pose
            )
        )

        if previous_score is not None:

            print(
                f"[{number}/{len(df)}] "
                f"{chembl_id}: "
                f"already completed "
                f"({previous_score:.3f})"
            )

            results.append({
                "parent_chembl_id":
                    chembl_id,

                "vina_affinity_kcal_mol":
                    previous_score,

                "status":
                    "completed",

                "pose_file":
                    str(output_pose),

                "log_file":
                    str(output_log),
            })

            continue


    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------

    if not ligand.exists():

        print(
            f"[{number}/{len(df)}] "
            f"{chembl_id}: "
            f"MISSING PDBQT"
        )

        results.append({
            "parent_chembl_id":
                chembl_id,

            "vina_affinity_kcal_mol":
                None,

            "status":
                "missing_ligand_pdbqt",

            "pose_file":
                None,

            "log_file":
                None,
        })

        continue


    # -------------------------------------------------------------------------
    # Vina command
    # -------------------------------------------------------------------------

    command = [
        vina,

        "--receptor",
        str(RECEPTOR),

        "--ligand",
        str(ligand),

        "--center_x",
        str(CENTER_X),

        "--center_y",
        str(CENTER_Y),

        "--center_z",
        str(CENTER_Z),

        "--size_x",
        str(SIZE_X),

        "--size_y",
        str(SIZE_Y),

        "--size_z",
        str(SIZE_Z),

        "--exhaustiveness",
        str(EXHAUSTIVENESS),

        "--num_modes",
        str(NUM_MODES),

        "--energy_range",
        str(ENERGY_RANGE),

        "--spacing",
        str(SPACING),

        "--seed",
        str(SEED),

        "--out",
        str(output_pose),
    ]


    print(
        f"[{number}/{len(df)}] "
        f"Docking {chembl_id}..."
    )


    process = subprocess.run(
        command,
        capture_output=True,
        text=True
    )


    # -------------------------------------------------------------------------
    # Save Vina log
    # -------------------------------------------------------------------------

    with open(
        output_log,
        "w",
        encoding="utf-8"
    ) as handle:

        handle.write("COMMAND\n")
        handle.write(
            " ".join(command)
            + "\n\n"
        )

        handle.write("VINA_STDOUT\n")
        handle.write(
            process.stdout
        )

        handle.write(
            "\n\nVINA_STDERR\n"
        )
        handle.write(
            process.stderr
        )


    # -------------------------------------------------------------------------
    # Vina failure
    # -------------------------------------------------------------------------

    if process.returncode != 0:

        print(
            f"    FAILED "
            f"(exit {process.returncode})"
        )

        results.append({
            "parent_chembl_id":
                chembl_id,

            "vina_affinity_kcal_mol":
                None,

            "status":
                f"vina_error_{process.returncode}",

            "pose_file":
                None,

            "log_file":
                str(output_log),
        })

        continue

    # -------------------------------------------------------------------------
    # Parse affinity
    # -------------------------------------------------------------------------

    best_affinity = (
        parse_best_affinity(
            output_pose
        )
    )


    if best_affinity is None:

        status = (
            "output_generated_"
            "score_not_parsed"
        )

        print(
            "    Output generated, "
            "but score not parsed"
        )

    else:

        status = "completed"

        print(
            f"    Best affinity: "
            f"{best_affinity:.3f} "
            f"kcal/mol"
        )


    results.append({
        "parent_chembl_id":
            chembl_id,

        "vina_affinity_kcal_mol":
            best_affinity,

        "status":
            status,

        "pose_file":
            str(output_pose),

        "log_file":
            str(output_log),
    })

    # -------------------------------------------------------------------------
    # Progressive checkpoint
    # -------------------------------------------------------------------------

    checkpoint = pd.DataFrame(
        results
    )

    checkpoint.to_csv(
        RESULTS_CSV,
        index=False
    )


# =============================================================================
# Final results
# =============================================================================

results_df = pd.DataFrame(
    results
)


# Add original drug information

final_df = df.merge(
    results_df,
    on="parent_chembl_id",
    how="left",
    validate="one_to_one"
)


# Rank successful compounds

successful_mask = (
    final_df["status"]
    == "completed"
)


final_df[
    "vina_rank"
] = pd.NA


successful = (
    final_df[
        successful_mask
    ]
    .sort_values(
        "vina_affinity_kcal_mol",
        ascending=True
    )
    .copy()
)


successful[
    "vina_rank"
] = range(
    1,
    len(successful) + 1
)


rank_map = successful.set_index(
    "parent_chembl_id"
)["vina_rank"]


final_df.loc[
    successful_mask,
    "vina_rank"
] = (
    final_df.loc[
        successful_mask,
        "parent_chembl_id"
    ]
    .map(rank_map)
)


final_df = final_df.sort_values(
    [
        "vina_rank",
        "parent_chembl_id"
    ],
    na_position="last"
)


final_df.to_csv(
    RESULTS_CSV,
    index=False
)


# =============================================================================
# Metadata
# =============================================================================

n_completed = int(
    (
        final_df["status"]
        == "completed"
    ).sum()
)


n_failed = int(
    len(final_df)
    - n_completed
)


metadata = {

    "date_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "vina_version":
        vina_version,

    "receptor":
        str(RECEPTOR),

    "box": {
        "center_x": CENTER_X,
        "center_y": CENTER_Y,
        "center_z": CENTER_Z,

        "size_x": SIZE_X,
        "size_y": SIZE_Y,
        "size_z": SIZE_Z,
    },

    "exhaustiveness":
        EXHAUSTIVENESS,

    "num_modes":
        NUM_MODES,

    "energy_range_kcal_mol":
        ENERGY_RANGE,

    "seed":
        SEED,

    "n_selected":
        int(len(df)),

    "n_completed":
        n_completed,

    "n_failed":
        n_failed,

    "pilot_limit":
        args.limit,

    "results_file":
        str(RESULTS_CSV),

    "pose_directory":
        str(POSE_DIR),

    "log_directory":
        str(LOG_DIR),
}


with open(
    METADATA_JSON,
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
print("VINA VIRTUAL SCREENING")
print("======================================")

print(
    f"Ligands selected:          "
    f"{len(df)}"
)

print(
    f"Completed:                 "
    f"{n_completed}"
)

print(
    f"Failed:                    "
    f"{n_failed}"
)


if n_completed:

    best = successful.iloc[0]

    print(
        f"Best ligand:               "
        f"{best['parent_chembl_id']}"
    )

    print(
        f"Best affinity:             "
        f"{best['vina_affinity_kcal_mol']:.3f} "
        f"kcal/mol"
    )


print("\nFiles generated:")

print(
    f"  {RESULTS_CSV}"
)

print(
    f"  {METADATA_JSON}"
)

print(
    f"  {POSE_DIR}/"
)

print(
    f"  {LOG_DIR}/"
)
