#!/usr/bin/env python3

"""
16_download_chembl_approved_drugs.py

Download the raw library of approved small-molecule drugs from ChEMBL.

Selection criteria:
    - max_phase = 4
    - molecule_type = "Small molecule"

No additional filtering is performed at this stage.
Salts, alternative forms, withdrawn drugs and molecules without structures
are retained so that exclusions can be documented explicitly later.

Outputs:
    - approved_drugs_chembl_raw.csv
    - approved_drugs_chembl_metadata.json
"""

import json
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests


# =============================================================================
# Configuration
# =============================================================================

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
STATUS_URL = "https://www.ebi.ac.uk/chembl/api/data/status.json"

OUTPUT_CSV = Path("approved_drugs_chembl_raw.csv")
OUTPUT_METADATA = Path("approved_drugs_chembl_metadata.json")

LIMIT = 1000
TIMEOUT = 60


# =============================================================================
# Helper functions
# =============================================================================

def get_nested(dictionary, key):
    """Return a value from a dictionary safely."""
    if isinstance(dictionary, dict):
        return dictionary.get(key)
    return None


def request_json(url, params=None):
    """Perform a GET request and return JSON."""
    response = requests.get(
        url,
        params=params,
        timeout=TIMEOUT
    )
    response.raise_for_status()
    return response.json()


# =============================================================================
# ChEMBL version
# =============================================================================

print("Checking ChEMBL version...")

try:
    status = request_json(STATUS_URL)
    chembl_version = (
        status.get("chembl_db_version")
        or status.get("chembl_version")
        or status.get("version")
        or "unknown"
    )

    print(f"ChEMBL database version: {chembl_version}")

except Exception as exc:
    print(f"WARNING: ChEMBL version could not be retrieved: {exc}")
    status = {}
    chembl_version = "unknown"


# =============================================================================
# Download approved small molecules
# =============================================================================

print("\nDownloading approved small molecules...")

rows = []
offset = 0
total_count = None

while True:

    params = {
        "max_phase": 4,
        "molecule_type": "Small molecule",
        "limit": LIMIT,
        "offset": offset
    }

    data = request_json(BASE_URL, params=params)

    molecules = data.get("molecules", [])

    if total_count is None:
        total_count = data.get("page_meta", {}).get("total_count")
        print(f"Records reported by ChEMBL: {total_count}")

    if not molecules:
        break

    for mol in molecules:

        structures = mol.get("molecule_structures") or {}
        properties = mol.get("molecule_properties") or {}
        hierarchy = mol.get("molecule_hierarchy") or {}

        rows.append({

            # Identification
            "chembl_id":
                mol.get("molecule_chembl_id"),

            "name":
                mol.get("pref_name"),

            "molecule_type":
                mol.get("molecule_type"),

            # Approval information
            "max_phase":
                mol.get("max_phase"),

            "first_approval":
                mol.get("first_approval"),

            "withdrawn_flag":
                mol.get("withdrawn_flag"),

            # Parent/alternative-form information
            "parent_chembl_id":
                hierarchy.get("parent_chembl_id"),

            # Structure
            "canonical_smiles":
                structures.get("canonical_smiles"),

            "standard_inchi":
                structures.get("standard_inchi"),

            "standard_inchi_key":
                structures.get("standard_inchi_key"),

            # ChEMBL physicochemical properties
            "full_mwt":
                properties.get("full_mwt"),

            "mw_freebase":
                properties.get("mw_freebase"),

            "alogp":
                properties.get("alogp"),

            "psa":
                properties.get("psa"),

            "hbd":
                properties.get("hbd"),

            "hba":
                properties.get("hba"),

            "num_ro5_violations":
                properties.get("num_ro5_violations")
        })

    offset += len(molecules)

    print(
        f"\rDownloaded {offset}"
        + (f"/{total_count}" if total_count else "")
        + " molecules",
        end=""
    )

    if total_count is not None and offset >= total_count:
        break


print("\nDownload completed.")


# =============================================================================
# Save raw dataset
# =============================================================================

df = pd.DataFrame(rows)

df = df.sort_values(
    by=["chembl_id"],
    na_position="last"
).reset_index(drop=True)

df.to_csv(
    OUTPUT_CSV,
    index=False
)


# =============================================================================
# Summary
# =============================================================================

n_total = len(df)

n_with_smiles = df["canonical_smiles"].notna().sum()

n_without_smiles = df["canonical_smiles"].isna().sum()

n_with_parent = df["parent_chembl_id"].notna().sum()

n_withdrawn = (
    pd.to_numeric(
        df["withdrawn_flag"],
        errors="coerce"
    )
    .fillna(0)
    .astype(int)
    .eq(1)
    .sum()
)


metadata = {
    "download_date_utc":
        datetime.utcnow().isoformat() + "Z",

    "chembl_database_version":
        chembl_version,

    "selection":
        {
            "max_phase": 4,
            "molecule_type": "Small molecule"
        },

    "n_records":
        int(n_total),

    "n_with_canonical_smiles":
        int(n_with_smiles),

    "n_without_canonical_smiles":
        int(n_without_smiles),

    "n_with_parent_chembl_id":
        int(n_with_parent),

    "n_withdrawn_flag":
        int(n_withdrawn),

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
print("RAW ChEMBL APPROVED-DRUG LIBRARY")
print("======================================")

print(f"Total records:          {n_total}")
print(f"With canonical SMILES:  {n_with_smiles}")
print(f"Without SMILES:         {n_without_smiles}")
print(f"With parent ID:         {n_with_parent}")
print(f"Withdrawn flag = 1:     {n_withdrawn}")

print("\nFiles generated:")
print(f"  {OUTPUT_CSV}")
print(f"  {OUTPUT_METADATA}")
