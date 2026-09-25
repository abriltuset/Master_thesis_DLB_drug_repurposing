#!/usr/bin/env python3

"""
16_predict_bbb_admet_ai.py

Predict blood-brain barrier penetration for the normalized ChEMBL
approved-drug library using ADMET-AI v1.4.0.

Input:
    approved_drugs_chembl37_parent_structured.csv

Outputs:
    approved_drugs_chembl37_bbb_predictions.csv
    admet_ai_all_predictions.csv
    approved_drugs_chembl37_bbb_metadata.json

Important:
    No BBB threshold is applied in this script.
    BBB_Martins probabilities are retained as continuous predictions.
"""

import json
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd
from admet_ai import ADMETModel


# =============================================================================
# Configuration
# =============================================================================

INPUT_CSV = Path(
    "approved_drugs_chembl37_parent_structured.csv"
)

OUTPUT_BBB = Path(
    "approved_drugs_chembl37_bbb_predictions.csv"
)

OUTPUT_ALL_ADMET = Path(
    "admet_ai_all_predictions.csv"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_bbb_metadata.json"
)

SMILES_COLUMN = "canonical_smiles"

BATCH_SIZE = 256


# =============================================================================
# Read library
# =============================================================================

print("Reading normalized ChEMBL library...")

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input file not found: {INPUT_CSV}"
    )

df = pd.read_csv(INPUT_CSV)

if SMILES_COLUMN not in df.columns:
    raise ValueError(
        f"Column '{SMILES_COLUMN}' not found."
    )


if df[SMILES_COLUMN].isna().any():
    raise ValueError(
        "Input contains missing canonical SMILES."
    )


print(f"Input compounds: {len(df)}")


# =============================================================================
# Unique SMILES
# =============================================================================

unique_smiles = (
    df[SMILES_COLUMN]
    .astype(str)
    .drop_duplicates()
    .tolist()
)

print(f"Unique SMILES:   {len(unique_smiles)}")


# =============================================================================
# Load ADMET-AI model
# =============================================================================

print("\nLoading ADMET-AI models...")

model = ADMETModel()

print("ADMET-AI models loaded.")


# =============================================================================
# Predictions
# =============================================================================

print("\nRunning ADMET predictions...")

prediction_batches = []


for start in range(
    0,
    len(unique_smiles),
    BATCH_SIZE
):

    batch = unique_smiles[
        start:start + BATCH_SIZE
    ]

    predictions = model.predict(
        smiles=batch
    )

    # ADMET-AI returns SMILES as DataFrame index
    predictions = predictions.copy()

    predictions.insert(
        0,
        "canonical_smiles",
        predictions.index.astype(str)
    )

    predictions = predictions.reset_index(
        drop=True
    )

    prediction_batches.append(
        predictions
    )

    completed = min(
        start + BATCH_SIZE,
        len(unique_smiles)
    )

    print(
        f"\rPredicted "
        f"{completed}/{len(unique_smiles)}",
        end=""
    )


print()


all_predictions = pd.concat(
    prediction_batches,
    ignore_index=True
)


# =============================================================================
# Validate BBB endpoint
# =============================================================================

BBB_COLUMN = "BBB_Martins"

if BBB_COLUMN not in all_predictions.columns:

    print(
        "\nAvailable ADMET-AI columns:"
    )

    for column in all_predictions.columns:
        print(f"  {column}")

    raise ValueError(
        f"Expected endpoint "
        f"'{BBB_COLUMN}' not found."
    )


# =============================================================================
# Save all ADMET predictions
# =============================================================================

all_predictions.to_csv(
    OUTPUT_ALL_ADMET,
    index=False
)


# =============================================================================
# Merge BBB predictions with ChEMBL library
# =============================================================================

bbb_predictions = all_predictions[
    [
        "canonical_smiles",
        BBB_COLUMN
    ]
].copy()


bbb_predictions = bbb_predictions.rename(
    columns={
        BBB_COLUMN:
            "bbb_martins_probability"
    }
)


result = df.merge(
    bbb_predictions,
    on="canonical_smiles",
    how="left",
    validate="many_to_one"
)


# =============================================================================
# Validate merge
# =============================================================================

n_missing_predictions = int(
    result[
        "bbb_martins_probability"
    ]
    .isna()
    .sum()
)


if n_missing_predictions > 0:

    print(
        f"\nWARNING: "
        f"{n_missing_predictions} compounds "
        f"have no BBB prediction."
    )


# =============================================================================
# Save BBB dataset
# =============================================================================

result.to_csv(
    OUTPUT_BBB,
    index=False
)


# =============================================================================
# Descriptive statistics
# =============================================================================

bbb = result[
    "bbb_martins_probability"
].dropna()


summary = {
    "n_predictions":
        int(len(bbb)),

    "mean":
        float(bbb.mean()),

    "median":
        float(bbb.median()),

    "min":
        float(bbb.min()),

    "max":
        float(bbb.max()),

    "q25":
        float(bbb.quantile(0.25)),

    "q75":
        float(bbb.quantile(0.75))
}


# Counts shown only descriptively.
# No filtering is applied.

descriptive_counts = {

    "probability_ge_0.50":
        int((bbb >= 0.50).sum()),

    "probability_ge_0.70":
        int((bbb >= 0.70).sum()),

    "probability_ge_0.80":
        int((bbb >= 0.80).sum()),

    "probability_ge_0.90":
        int((bbb >= 0.90).sum())
}


# =============================================================================
# Metadata
# =============================================================================

try:
    admet_version = version(
        "admet-ai"
    )
except Exception:
    admet_version = "unknown"


metadata = {

    "date_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "input_file":
        str(INPUT_CSV),

    "admet_ai_version":
        admet_version,

    "endpoint":
        "BBB_Martins",

    "n_input_compounds":
        int(len(df)),

    "n_unique_smiles":
        int(len(unique_smiles)),

    "n_bbb_predictions":
        int(len(bbb)),

    "n_missing_bbb_predictions":
        n_missing_predictions,

    "bbb_probability_summary":
        summary,

    "descriptive_counts":
        descriptive_counts,

    "bbb_threshold_applied":
        False,

    "outputs": {
        "bbb_predictions":
            str(OUTPUT_BBB),

        "all_admet_predictions":
            str(OUTPUT_ALL_ADMET)
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
print("ADMET-AI BBB PREDICTION")
print("======================================")

print(f"ADMET-AI version:       {admet_version}")
print(f"Input compounds:        {len(df)}")
print(f"Unique SMILES:          {len(unique_smiles)}")
print(f"BBB predictions:        {len(bbb)}")
print(f"Missing predictions:    {n_missing_predictions}")

print("\nBBB_Martins probability:")

print(
    f"Mean:                    "
    f"{summary['mean']:.3f}"
)

print(
    f"Median:                  "
    f"{summary['median']:.3f}"
)

print(
    f"Range:                   "
    f"{summary['min']:.3f} - "
    f"{summary['max']:.3f}"
)

print(
    f"Q25-Q75:                 "
    f"{summary['q25']:.3f} - "
    f"{summary['q75']:.3f}"
)

print("\nDescriptive counts only:")

print(
    f"P >= 0.50:               "
    f"{descriptive_counts['probability_ge_0.50']}"
)

print(
    f"P >= 0.70:               "
    f"{descriptive_counts['probability_ge_0.70']}"
)

print(
    f"P >= 0.80:               "
    f"{descriptive_counts['probability_ge_0.80']}"
)

print(
    f"P >= 0.90:               "
    f"{descriptive_counts['probability_ge_0.90']}"
)

print("\nNo BBB threshold was applied.")

print("\nFiles generated:")

print(
    f"  {OUTPUT_BBB}"
)

print(
    f"  {OUTPUT_ALL_ADMET}"
)

print(
    f"  {OUTPUT_METADATA}"
)
