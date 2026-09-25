#!/usr/bin/env python3

"""
20_filter_docking_compatible_drugs.py

Apply conservative technical compatibility criteria before ligand preparation.

Input:
    approved_drugs_chembl37_bbb_positive_docking_audit.csv

Exclusion criteria:
    1. More than one disconnected molecular fragment.
    2. No carbon atom.
    3. Presence of elements outside the element set used
       in the standard Vina screening workflow.

No MW, charge, BBB-confidence or safety filtering is applied here.

Outputs:
    approved_drugs_chembl37_docking_compatible.csv
    approved_drugs_chembl37_docking_excluded.csv
"""

from pathlib import Path
import pandas as pd


INPUT = Path(
    "approved_drugs_chembl37_bbb_positive_docking_audit.csv"
)

OUTPUT_KEEP = Path(
    "approved_drugs_chembl37_docking_compatible.csv"
)

OUTPUT_EXCLUDED = Path(
    "approved_drugs_chembl37_docking_excluded.csv"
)


# Elements accepted for this standard small-molecule Vina workflow.
#
# Selenium is retained at this stage and its practical compatibility
# is evaluated during subsequent PDBQT generation with Meeko.
# Silicon is explicitly supported.
#
# Boron is not included because it is not represented as
# a standard accepted atom type in the current Vina implementation.

ALLOWED_ELEMENTS = {
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "Cl",
    "Br",
    "I",
    "Si",
    "Se",
}


print("Reading docking audit...")

df = pd.read_csv(INPUT)

print(f"Input compounds: {len(df)}")


# =============================================================================
# Helpers
# =============================================================================

def parse_elements(value):

    if pd.isna(value):
        return set()

    return {
        x.strip()
        for x in str(value).split(";")
        if x.strip()
    }


def exclusion_reason(row):

    reasons = []

    elements = parse_elements(
        row["elements"]
    )

    # Disconnected structure
    if int(row["n_fragments"]) > 1:
        reasons.append(
            "multiple_disconnected_fragments"
        )

    # Require an organic scaffold
    if "C" not in elements:
        reasons.append(
            "no_carbon"
        )

    # Unsupported/non-standard elements
    unsupported = sorted(
        elements - ALLOWED_ELEMENTS
    )

    if unsupported:
        reasons.append(
            "unsupported_elements:"
            + ",".join(unsupported)
        )

    return ";".join(reasons)


# =============================================================================
# Apply criteria
# =============================================================================

df["technical_exclusion_reason"] = (
    df.apply(
        exclusion_reason,
        axis=1
    )
)

df["docking_compatible"] = (
    df["technical_exclusion_reason"] == ""
)


keep = (
    df[df["docking_compatible"]]
    .copy()
    .reset_index(drop=True)
)

excluded = (
    df[~df["docking_compatible"]]
    .copy()
    .reset_index(drop=True)
)


# =============================================================================
# Exclusion counts
# =============================================================================

n_multifragment = int(
    (
        df["n_fragments"] > 1
    ).sum()
)


n_no_carbon = 0
n_unsupported = 0

unsupported_counter = {}


for _, row in df.iterrows():

    elements = parse_elements(
        row["elements"]
    )

    if "C" not in elements:
        n_no_carbon += 1

    unsupported = (
        elements - ALLOWED_ELEMENTS
    )

    if unsupported:

        n_unsupported += 1

        for element in unsupported:

            unsupported_counter[element] = (
                unsupported_counter.get(
                    element,
                    0
                ) + 1
            )


# =============================================================================
# Export
# =============================================================================

keep.to_csv(
    OUTPUT_KEEP,
    index=False
)

excluded.to_csv(
    OUTPUT_EXCLUDED,
    index=False
)


# =============================================================================
# Report
# =============================================================================

print("\n======================================")
print("DOCKING TECHNICAL FILTER")
print("======================================")

print(
    f"Input compounds:             "
    f"{len(df)}"
)

print(
    f"Docking compatible:          "
    f"{len(keep)}"
)

print(
    f"Excluded:                    "
    f"{len(excluded)}"
)

print(
    f"Multiple fragments:          "
    f"{n_multifragment}"
)

print(
    f"No carbon:                   "
    f"{n_no_carbon}"
)

print(
    f"Unsupported elements:        "
    f"{n_unsupported}"
)


print("\nUnsupported-element counts:")

for element, count in sorted(
    unsupported_counter.items()
):
    print(
        f"  {element}: {count}"
    )


print(
    "\nNote: categories overlap; "
    "exclusion counts therefore do not sum "
    "to the total excluded."
)

print("\nFiles generated:")

print(
    f"  {OUTPUT_KEEP}"
)

print(
    f"  {OUTPUT_EXCLUDED}"
)
