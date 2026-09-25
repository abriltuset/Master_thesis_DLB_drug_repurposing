#!/usr/bin/env python3
# Note: the outputs used in the thesis were generated from ChEMBL 37.
# Because the script queries the live API, later reruns may reflect newer releases.
"""
17_normalize_chembl_parent_drugs.py

Normalize the raw ChEMBL approved-drug library to parent compounds.

Input:
    approved_drugs_chembl_raw.csv

Outputs:
    approved_drugs_chembl37_parent_all.csv
    approved_drugs_chembl37_parent_structured.csv
    approved_drugs_chembl37_parent_metadata.json

Steps:
    1. Assign each approved ChEMBL record to its parent compound.
    2. Retrieve parent-compound information from the ChEMBL API.
    3. Collapse alternative forms/salts into one parent record.
    4. Retain traceability to all original approved ChEMBL IDs.
    5. Generate a second library restricted to compounds
       with a parent canonical SMILES.

No BBB/CNS filtering is performed here.
Withdrawn drugs are NOT removed at this stage.
"""

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


# =============================================================================
# Configuration
# =============================================================================

INPUT_CSV = Path("approved_drugs_chembl_raw.csv")

OUTPUT_ALL = Path(
    "approved_drugs_chembl37_parent_all.csv"
)

OUTPUT_STRUCTURED = Path(
    "approved_drugs_chembl37_parent_structured.csv"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_parent_metadata.json"
)

MOLECULE_URL = (
    "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
)

MOLECULE_DETAIL_URL = (
    "https://www.ebi.ac.uk/chembl/api/data/molecule"
)

BATCH_SIZE = 100
TIMEOUT = 60
MAX_RETRIES = 3


# =============================================================================
# Helpers
# =============================================================================

def clean_string(value):
    """Convert missing/empty values to None."""
    if pd.isna(value):
        return None

    value = str(value).strip()

    if not value or value.lower() == "nan":
        return None

    return value


def flag_to_int(value):
    """Convert ChEMBL boolean/integer flags to 0/1."""

    if pd.isna(value):
        return 0

    if isinstance(value, bool):
        return int(value)

    text = str(value).strip().lower()

    if text in {"1", "true", "yes"}:
        return 1

    return 0


def request_json(url, params=None):
    """GET JSON with simple retry logic."""

    last_exception = None

    for attempt in range(1, MAX_RETRIES + 1):

        try:
            response = requests.get(
                url,
                params=params,
                timeout=TIMEOUT
            )

            response.raise_for_status()

            return response.json()

        except Exception as exc:
            last_exception = exc

            print(
                f"\nRequest failed "
                f"(attempt {attempt}/{MAX_RETRIES}): {exc}"
            )

            if attempt < MAX_RETRIES:
                time.sleep(2)

    raise RuntimeError(
        f"ChEMBL request failed after "
        f"{MAX_RETRIES} attempts"
    ) from last_exception


def first_non_null(series):
    """Return the first non-empty value in a pandas Series."""

    for value in series:

        cleaned = clean_string(value)

        if cleaned is not None:
            return cleaned

    return None


# =============================================================================
# Read input
# =============================================================================

print("Reading raw ChEMBL library...")

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input file not found: {INPUT_CSV}"
    )

raw = pd.read_csv(INPUT_CSV)

required_columns = {
    "chembl_id",
    "name",
    "parent_chembl_id",
    "withdrawn_flag"
}

missing_columns = (
    required_columns - set(raw.columns)
)

if missing_columns:
    raise ValueError(
        "Missing required columns: "
        + ", ".join(sorted(missing_columns))
    )


print(f"Raw records: {len(raw)}")


# =============================================================================
# Resolve parent IDs
# =============================================================================

def resolve_parent(row):

    parent = clean_string(
        row["parent_chembl_id"]
    )

    if parent is not None:
        return parent

    return clean_string(
        row["chembl_id"]
    )


raw["resolved_parent_chembl_id"] = raw.apply(
    resolve_parent,
    axis=1
)


parent_ids = sorted(
    raw["resolved_parent_chembl_id"]
    .dropna()
    .unique()
)


print(
    f"Unique parent ChEMBL IDs: "
    f"{len(parent_ids)}"
)


# =============================================================================
# Retrieve parent molecules from ChEMBL
# =============================================================================

print("\nRetrieving parent compounds from ChEMBL...")

parent_records = {}


for start in range(
    0,
    len(parent_ids),
    BATCH_SIZE
):

    batch = parent_ids[
        start:start + BATCH_SIZE
    ]

    params = {
        "molecule_chembl_id__in":
            ",".join(batch),

        "limit":
            BATCH_SIZE
    }

    data = request_json(
        MOLECULE_URL,
        params=params
    )

    molecules = data.get(
        "molecules",
        []
    )

    for molecule in molecules:

        chembl_id = molecule.get(
            "molecule_chembl_id"
        )

        if chembl_id:
            parent_records[
                chembl_id
            ] = molecule

    downloaded = min(
        start + BATCH_SIZE,
        len(parent_ids)
    )

    print(
        f"\rRetrieved "
        f"{downloaded}/{len(parent_ids)} "
        f"parent IDs",
        end=""
    )


print()


# =============================================================================
# Fallback for missing parent records
# =============================================================================

missing_parent_ids = [
    parent_id
    for parent_id in parent_ids
    if parent_id not in parent_records
]


if missing_parent_ids:

    print(
        f"\n{len(missing_parent_ids)} parent IDs "
        f"were not returned in batch queries."
    )

    print(
        "Trying individual ChEMBL queries..."
    )

    for i, parent_id in enumerate(
        missing_parent_ids,
        start=1
    ):

        url = (
            f"{MOLECULE_DETAIL_URL}/"
            f"{parent_id}.json"
        )

        try:
            molecule = request_json(url)

            if molecule:
                parent_records[
                    parent_id
                ] = molecule

        except Exception as exc:

            print(
                f"\nWARNING: could not retrieve "
                f"{parent_id}: {exc}"
            )

        print(
            f"\rFallback queries: "
            f"{i}/{len(missing_parent_ids)}",
            end=""
        )

    print()


# =============================================================================
# Build one row per parent compound
# =============================================================================

print("\nCollapsing approved forms into parent compounds...")


rows = []


for parent_id, family in raw.groupby(
    "resolved_parent_chembl_id",
    sort=True
):

    molecule = parent_records.get(
        parent_id,
        {}
    )

    structures = (
        molecule.get("molecule_structures")
        or {}
    )

    properties = (
        molecule.get("molecule_properties")
        or {}
    )

    parent_name = clean_string(
        molecule.get("pref_name")
    )

    # If the parent has no preferred name,
    # retain one original approved-form name
    if parent_name is None:
        parent_name = first_non_null(
            family["name"]
        )

    source_ids = sorted(
        set(
            family["chembl_id"]
            .dropna()
            .astype(str)
        )
    )

    source_names = sorted(
        {
            str(x)
            for x in family["name"].dropna()
            if str(x).strip()
        }
    )

    raw_withdrawn = max(
        flag_to_int(x)
        for x in family["withdrawn_flag"]
    )

    parent_withdrawn = flag_to_int(
        molecule.get("withdrawn_flag")
    )

    family_withdrawn = max(
        raw_withdrawn,
        parent_withdrawn
    )

    canonical_smiles = clean_string(
        structures.get(
            "canonical_smiles"
        )
    )

    rows.append({

        # Parent identity
        "parent_chembl_id":
            parent_id,

        "name":
            parent_name,

        "molecule_type":
            molecule.get(
                "molecule_type"
            ),

        # Approval
        "max_phase":
            molecule.get(
                "max_phase"
            ),

        "first_approval":
            molecule.get(
                "first_approval"
            ),

        # Safety flag
        "withdrawn_flag":
            family_withdrawn,

        # Structure
        "canonical_smiles":
            canonical_smiles,

        "standard_inchi":
            structures.get(
                "standard_inchi"
            ),

        "standard_inchi_key":
            structures.get(
                "standard_inchi_key"
            ),

        # Physicochemical properties
        "full_mwt":
            properties.get(
                "full_mwt"
            ),

        "mw_freebase":
            properties.get(
                "mw_freebase"
            ),

        "alogp":
            properties.get(
                "alogp"
            ),

        "psa":
            properties.get(
                "psa"
            ),

        "hbd":
            properties.get(
                "hbd"
            ),

        "hba":
            properties.get(
                "hba"
            ),

        "num_ro5_violations":
            properties.get(
                "num_ro5_violations"
            ),

        # Traceability
        "n_approved_forms":
            len(source_ids),

        "approved_form_chembl_ids":
            ";".join(source_ids),

        "approved_form_names":
            ";".join(source_names),

        # Structure availability
        "has_parent_smiles":
            int(
                canonical_smiles
                is not None
            )
    })


parent_df = pd.DataFrame(rows)


# =============================================================================
# Sort and export
# =============================================================================

parent_df = parent_df.sort_values(
    by=[
        "name",
        "parent_chembl_id"
    ],
    na_position="last"
).reset_index(drop=True)


parent_df.to_csv(
    OUTPUT_ALL,
    index=False
)


structured_df = (
    parent_df[
        parent_df["has_parent_smiles"] == 1
    ]
    .copy()
    .reset_index(drop=True)
)


structured_df.to_csv(
    OUTPUT_STRUCTURED,
    index=False
)


# =============================================================================
# Summary metrics
# =============================================================================

n_raw = len(raw)

n_parent = len(parent_df)

n_structured = len(
    structured_df
)

n_without_structure = (
    n_parent - n_structured
)

n_withdrawn_parent = int(
    parent_df[
        "withdrawn_flag"
    ]
    .eq(1)
    .sum()
)

n_multiform = int(
    parent_df[
        "n_approved_forms"
    ]
    .gt(1)
    .sum()
)


metadata = {

    "date_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "input_file":
        str(INPUT_CSV),

    "raw_records":
        int(n_raw),

    "unique_parent_compounds":
        int(n_parent),

    "parents_with_structure":
        int(n_structured),

    "parents_without_structure":
        int(n_without_structure),

    "parent_families_with_multiple_approved_forms":
        int(n_multiform),

    "parents_with_withdrawn_flag":
        int(n_withdrawn_parent),

    "withdrawn_drugs_removed":
        False,

    "bbb_filter_applied":
        False,

    "outputs": {
        "all_parents":
            str(OUTPUT_ALL),

        "structured_parents":
            str(OUTPUT_STRUCTURED)
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
# Final report
# =============================================================================

print("\n======================================")
print("ChEMBL PARENT NORMALIZATION")
print("======================================")

print(
    f"Raw approved records:       "
    f"{n_raw}"
)

print(
    f"Unique parent compounds:    "
    f"{n_parent}"
)

print(
    f"Parents with SMILES:        "
    f"{n_structured}"
)

print(
    f"Parents without SMILES:     "
    f"{n_without_structure}"
)

print(
    f"Families with >1 form:      "
    f"{n_multiform}"
)

print(
    f"Withdrawn flag = 1:         "
    f"{n_withdrawn_parent}"
)

print("\nFiles generated:")

print(
    f"  {OUTPUT_ALL}"
)

print(
    f"  {OUTPUT_STRUCTURED}"
)

print(
    f"  {OUTPUT_METADATA}"
)
