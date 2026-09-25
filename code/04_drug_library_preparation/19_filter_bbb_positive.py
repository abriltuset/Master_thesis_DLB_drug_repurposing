#!/usr/bin/env python3

"""
17_filter_bbb_positive.py

Filter the approved-drug library using the BBB_Martins prediction.

Input:
    approved_drugs_chembl37_bbb_predictions.csv

Selection criterion:
    BBB_Martins probability >= 0.50

Additional annotation:
    high_confidence_bbb = True if probability >= 0.80

No additional physicochemical or docking-related filters are applied here.

Outputs:
    approved_drugs_chembl37_bbb_positive.csv
    approved_drugs_chembl37_bbb_filter_metadata.json
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

INPUT_CSV = Path(
    "approved_drugs_chembl37_bbb_predictions.csv"
)

OUTPUT_CSV = Path(
    "approved_drugs_chembl37_bbb_positive.csv"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_bbb_filter_metadata.json"
)

BBB_COLUMN = "bbb_martins_probability"

BBB_THRESHOLD = 0.50
HIGH_CONFIDENCE_THRESHOLD = 0.80


# =============================================================================
# Read input
# =============================================================================

print("Reading BBB prediction library...")

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input file not found: {INPUT_CSV}"
    )

df = pd.read_csv(INPUT_CSV)

if BBB_COLUMN not in df.columns:
    raise ValueError(
        f"Column '{BBB_COLUMN}' not found."
    )

print(f"Input compounds: {len(df)}")


# =============================================================================
# Validate probabilities
# =============================================================================

bbb = pd.to_numeric(
    df[BBB_COLUMN],
    errors="coerce"
)

n_missing = int(
    bbb.isna().sum()
)

if n_missing > 0:
    raise ValueError(
        f"{n_missing} compounds have missing "
        f"BBB probabilities."
    )


if ((bbb < 0) | (bbb > 1)).any():
    raise ValueError(
        "BBB probabilities outside [0, 1] detected."
    )


df[BBB_COLUMN] = bbb


# =============================================================================
# BBB filtering
# =============================================================================

df["bbb_positive"] = (
    df[BBB_COLUMN] >= BBB_THRESHOLD
)

df["high_confidence_bbb"] = (
    df[BBB_COLUMN] >= HIGH_CONFIDENCE_THRESHOLD
)


selected = (
    df[df["bbb_positive"]]
    .copy()
    .sort_values(
        BBB_COLUMN,
        ascending=False
    )
    .reset_index(drop=True)
)


# =============================================================================
# Summary
# =============================================================================

n_input = len(df)

n_selected = len(selected)

n_excluded = (
    n_input - n_selected
)

n_high_confidence = int(
    selected[
        "high_confidence_bbb"
    ].sum()
)


# =============================================================================
# Export
# =============================================================================

selected.to_csv(
    OUTPUT_CSV,
    index=False
)


metadata = {

    "date_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "input_file":
        str(INPUT_CSV),

    "bbb_endpoint":
        "BBB_Martins",

    "bbb_threshold":
        BBB_THRESHOLD,

    "high_confidence_threshold":
        HIGH_CONFIDENCE_THRESHOLD,

    "n_input_compounds":
        int(n_input),

    "n_bbb_positive":
        int(n_selected),

    "n_bbb_negative":
        int(n_excluded),

    "n_high_confidence_bbb":
        int(n_high_confidence),

    "additional_physicochemical_filters_applied":
        False,

    "output_file":
        str(OUTPUT_CSV)
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
print("BBB POSITIVE FILTER")
print("======================================")

print(
    f"Input compounds:          "
    f"{n_input}"
)

print(
    f"BBB+ (P >= {BBB_THRESHOLD:.2f}):       "
    f"{n_selected}"
)

print(
    f"BBB- (P < {BBB_THRESHOLD:.2f}):        "
    f"{n_excluded}"
)

print(
    f"High confidence "
    f"(P >= {HIGH_CONFIDENCE_THRESHOLD:.2f}): "
    f"{n_high_confidence}"
)

print(
    f"Retention:                "
    f"{100 * n_selected / n_input:.1f}%"
)

print("\nNo other filters were applied.")

print("\nFiles generated:")
print(f"  {OUTPUT_CSV}")
print(f"  {OUTPUT_METADATA}")
