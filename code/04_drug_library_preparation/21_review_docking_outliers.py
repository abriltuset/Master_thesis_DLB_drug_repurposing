#!/usr/bin/env python3

"""
19_review_docking_outliers.py

Summarize compounds requiring review before ligand preparation.

NO compounds are removed.

Input:
    approved_drugs_chembl37_bbb_positive_docking_audit.csv

Outputs:
    approved_drugs_chembl37_manual_review.csv
    approved_drugs_chembl37_unusual_element_counts.csv
"""

from pathlib import Path
from collections import Counter

import pandas as pd


INPUT = Path(
    "approved_drugs_chembl37_bbb_positive_docking_audit.csv"
)

OUTPUT_REVIEW = Path(
    "approved_drugs_chembl37_manual_review.csv"
)

OUTPUT_ELEMENTS = Path(
    "approved_drugs_chembl37_unusual_element_counts.csv"
)


print("Reading docking audit...")

df = pd.read_csv(INPUT)

print(f"Input compounds: {len(df)}")


# =============================================================================
# Additional descriptive flags
# =============================================================================

df["mw_lt_100_flag"] = (
    df["mol_wt_rdkit"] < 100
)

df["mw_gt_700_flag"] = (
    df["mol_wt_rdkit"] > 700
)

df["abs_charge_gt_1_flag"] = (
    df["abs_formal_charge"] > 1
)


# =============================================================================
# Generate review reasons
# =============================================================================

def get_reasons(row):

    reasons = []

    if not row["rdkit_parse_ok"]:
        reasons.append("RDKit_parse_failure")

    if row["multifragment_flag"]:
        reasons.append("multiple_fragments")

    if row["special_element_flag"]:
        reasons.append("B_or_Si")

    if row["other_element_flag"]:
        reasons.append("unusual_element")

    if row["torsion_gt32_flag"]:
        reasons.append("rotatable_bonds_gt32")

    if row["mw_lt_100_flag"]:
        reasons.append("MW_lt_100")

    if row["mw_gt_700_flag"]:
        reasons.append("MW_gt_700")

    if row["abs_charge_gt_1_flag"]:
        reasons.append("abs_charge_gt_1")

    return ";".join(reasons)


df["review_reasons"] = df.apply(
    get_reasons,
    axis=1
)

df["review_required"] = (
    df["review_reasons"] != ""
)


review = (
    df[df["review_required"]]
    .copy()
    .sort_values(
        [
            "other_element_flag",
            "multifragment_flag",
            "mol_wt_rdkit"
        ],
        ascending=[
            False,
            False,
            False
        ]
    )
    .reset_index(drop=True)
)


# =============================================================================
# Count unusual elements
# =============================================================================

element_counter = Counter()

for value in df["other_elements"].dropna():

    for element in str(value).split(";"):

        element = element.strip()

        if element:
            element_counter[element] += 1


element_df = pd.DataFrame(
    [
        {
            "element": element,
            "n_compounds": count
        }
        for element, count
        in element_counter.most_common()
    ]
)


# =============================================================================
# Export
# =============================================================================

review.to_csv(
    OUTPUT_REVIEW,
    index=False
)

element_df.to_csv(
    OUTPUT_ELEMENTS,
    index=False
)


# =============================================================================
# Report
# =============================================================================

print("\n======================================")
print("DOCKING OUTLIER REVIEW")
print("======================================")

print(f"Total compounds:             {len(df)}")
print(f"Require review:              {len(review)}")

print(
    f"Multiple fragments:          "
    f"{int(df['multifragment_flag'].sum())}"
)

print(
    f"B / Si:                     "
    f"{int(df['special_element_flag'].sum())}"
)

print(
    f"Other unusual elements:      "
    f"{int(df['other_element_flag'].sum())}"
)

print(
    f"MW < 100 Da:                 "
    f"{int(df['mw_lt_100_flag'].sum())}"
)

print(
    f"MW > 700 Da:                 "
    f"{int(df['mw_gt_700_flag'].sum())}"
)

print(
    f"|formal charge| > 1:         "
    f"{int(df['abs_charge_gt_1_flag'].sum())}"
)


print("\nUnusual elements:")

if len(element_df) == 0:
    print("  None")
else:
    for _, row in element_df.iterrows():
        print(
            f"  {row['element']}: "
            f"{row['n_compounds']}"
        )


print("\nFiles generated:")
print(f"  {OUTPUT_REVIEW}")
print(f"  {OUTPUT_ELEMENTS}")
