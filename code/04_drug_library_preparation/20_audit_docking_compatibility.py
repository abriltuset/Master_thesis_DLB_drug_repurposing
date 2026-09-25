#!/usr/bin/env python3

"""
18_audit_docking_compatibility.py

Audit BBB-positive approved drugs before ligand preparation and docking.

Input:
    approved_drugs_chembl37_bbb_positive.csv

This script DOES NOT exclude compounds.

It evaluates:
    - SMILES validity
    - number of disconnected fragments
    - elemental composition
    - molecular weight
    - heavy atoms
    - rotatable bonds
    - formal charge
    - TPSA
    - cLogP
    - HBD / HBA

Flags are generated for compounds that may require manual review
before 3D/PDBQT preparation.

Output:
    approved_drugs_chembl37_bbb_positive_docking_audit.csv
    approved_drugs_chembl37_docking_audit_metadata.json
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from rdkit import Chem
from rdkit.Chem import (
    Crippen,
    Descriptors,
    Lipinski,
    rdMolDescriptors,
)


# =============================================================================
# Configuration
# =============================================================================

INPUT_CSV = Path(
    "approved_drugs_chembl37_bbb_positive.csv"
)

OUTPUT_CSV = Path(
    "approved_drugs_chembl37_bbb_positive_docking_audit.csv"
)

OUTPUT_METADATA = Path(
    "approved_drugs_chembl37_docking_audit_metadata.json"
)

SMILES_COLUMN = "canonical_smiles"


# Common elements handled routinely in conventional organic ligand docking.
# B and Si are kept separately because AutoDock/Meeko can require
# special handling depending on the workflow.
ROUTINE_ELEMENTS = {
    "C", "N", "O", "F", "P",
    "S", "Cl", "Br", "I"
}

SPECIAL_ELEMENTS = {
    "B", "Si"
}


# =============================================================================
# Load input
# =============================================================================

print("Reading BBB-positive drug library...")

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input file not found: {INPUT_CSV}"
    )

df = pd.read_csv(INPUT_CSV)

if SMILES_COLUMN not in df.columns:
    raise ValueError(
        f"Column '{SMILES_COLUMN}' not found."
    )

print(f"Input compounds: {len(df)}")


# =============================================================================
# Audit function
# =============================================================================

def audit_smiles(smiles):

    result = {
        "rdkit_parse_ok": False,
        "n_fragments": None,
        "multifragment_flag": False,
        "elements": None,
        "special_elements": None,
        "other_elements": None,
        "special_element_flag": False,
        "other_element_flag": False,
        "heavy_atoms": None,
        "rotatable_bonds": None,
        "torsion_gt32_flag": False,
        "formal_charge": None,
        "abs_formal_charge": None,
        "mol_wt_rdkit": None,
        "tpsa_rdkit": None,
        "clogp_rdkit": None,
        "hbd_rdkit": None,
        "hba_rdkit": None,
        "ring_count": None,
        "manual_review_flag": True,
    }

    if pd.isna(smiles):
        return result

    smiles = str(smiles).strip()

    if not smiles:
        return result

    # -------------------------------------------------------------------------
    # Parse molecule
    # -------------------------------------------------------------------------

    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None

    if mol is None:
        return result

    result["rdkit_parse_ok"] = True

    # -------------------------------------------------------------------------
    # Fragments
    # -------------------------------------------------------------------------

    fragments = Chem.GetMolFrags(
        mol,
        asMols=False,
        sanitizeFrags=False
    )

    n_fragments = len(fragments)

    result["n_fragments"] = n_fragments
    result["multifragment_flag"] = (
        n_fragments > 1
    )

    # -------------------------------------------------------------------------
    # Elements
    # -------------------------------------------------------------------------

    elements = sorted({
        atom.GetSymbol()
        for atom in mol.GetAtoms()
    })

    special = sorted(
        set(elements) & SPECIAL_ELEMENTS
    )

    other = sorted(
        set(elements)
        - ROUTINE_ELEMENTS
        - SPECIAL_ELEMENTS
    )

    result["elements"] = ";".join(elements)

    result["special_elements"] = (
        ";".join(special)
        if special
        else ""
    )

    result["other_elements"] = (
        ";".join(other)
        if other
        else ""
    )

    result["special_element_flag"] = (
        len(special) > 0
    )

    result["other_element_flag"] = (
        len(other) > 0
    )

    # -------------------------------------------------------------------------
    # Molecular descriptors
    # -------------------------------------------------------------------------

    heavy_atoms = mol.GetNumHeavyAtoms()

    rotatable = Lipinski.NumRotatableBonds(
        mol
    )

    formal_charge = sum(
        atom.GetFormalCharge()
        for atom in mol.GetAtoms()
    )

    result["heavy_atoms"] = int(
        heavy_atoms
    )

    result["rotatable_bonds"] = int(
        rotatable
    )

    # >32 is retained only as a technical warning.
    # No molecule is excluded here.
    result["torsion_gt32_flag"] = (
        rotatable > 32
    )

    result["formal_charge"] = int(
        formal_charge
    )

    result["abs_formal_charge"] = abs(
        int(formal_charge)
    )

    result["mol_wt_rdkit"] = float(
        Descriptors.MolWt(mol)
    )

    result["tpsa_rdkit"] = float(
        rdMolDescriptors.CalcTPSA(mol)
    )

    result["clogp_rdkit"] = float(
        Crippen.MolLogP(mol)
    )

    result["hbd_rdkit"] = int(
        Lipinski.NumHDonors(mol)
    )

    result["hba_rdkit"] = int(
        Lipinski.NumHAcceptors(mol)
    )

    result["ring_count"] = int(
        rdMolDescriptors.CalcNumRings(mol)
    )

    # -------------------------------------------------------------------------
    # Manual review flag
    # -------------------------------------------------------------------------

    result["manual_review_flag"] = any([
        result["multifragment_flag"],
        result["special_element_flag"],
        result["other_element_flag"],
        result["torsion_gt32_flag"],
    ])

    return result


# =============================================================================
# Run audit
# =============================================================================

print("\nAuditing structures...")

audit_rows = []

for i, smiles in enumerate(
    df[SMILES_COLUMN],
    start=1
):

    audit_rows.append(
        audit_smiles(smiles)
    )

    if i % 100 == 0 or i == len(df):

        print(
            f"\rAudited {i}/{len(df)}",
            end=""
        )

print()


audit_df = pd.DataFrame(
    audit_rows
)

result = pd.concat(
    [
        df.reset_index(drop=True),
        audit_df
    ],
    axis=1
)


# =============================================================================
# Summary
# =============================================================================

n_total = len(result)

n_parse_failed = int(
    (~result["rdkit_parse_ok"]).sum()
)

n_multifragment = int(
    result["multifragment_flag"].sum()
)

n_special_elements = int(
    result["special_element_flag"].sum()
)

n_other_elements = int(
    result["other_element_flag"].sum()
)

n_torsion_gt32 = int(
    result["torsion_gt32_flag"].sum()
)

n_manual_review = int(
    result["manual_review_flag"].sum()
)


# Additional descriptive counts.
# These are NOT exclusion criteria.

valid = result[
    result["rdkit_parse_ok"]
].copy()

descriptor_summary = {

    "molecular_weight": {
        "median":
            float(
                valid["mol_wt_rdkit"].median()
            ),
        "min":
            float(
                valid["mol_wt_rdkit"].min()
            ),
        "max":
            float(
                valid["mol_wt_rdkit"].max()
            ),
    },

    "rotatable_bonds": {
        "median":
            float(
                valid["rotatable_bonds"].median()
            ),
        "max":
            int(
                valid["rotatable_bonds"].max()
            ),
    },

    "formal_charge": {
        "min":
            int(
                valid["formal_charge"].min()
            ),
        "max":
            int(
                valid["formal_charge"].max()
            ),
    }
}


# =============================================================================
# Save
# =============================================================================

result.to_csv(
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

    "n_input_compounds":
        int(n_total),

    "n_rdkit_parse_failed":
        int(n_parse_failed),

    "n_multifragment":
        int(n_multifragment),

    "n_special_elements_B_or_Si":
        int(n_special_elements),

    "n_other_elements":
        int(n_other_elements),

    "n_rotatable_bonds_gt_32":
        int(n_torsion_gt32),

    "n_manual_review":
        int(n_manual_review),

    "descriptor_summary":
        descriptor_summary,

    "compounds_removed":
        0,

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
print("DOCKING COMPATIBILITY AUDIT")
print("======================================")

print(
    f"Input compounds:             "
    f"{n_total}"
)

print(
    f"RDKit parse failures:        "
    f"{n_parse_failed}"
)

print(
    f"Multiple fragments:          "
    f"{n_multifragment}"
)

print(
    f"B/Si present:                "
    f"{n_special_elements}"
)

print(
    f"Other unusual elements:      "
    f"{n_other_elements}"
)

print(
    f"Rotatable bonds >32:         "
    f"{n_torsion_gt32}"
)

print(
    f"Manual-review compounds:     "
    f"{n_manual_review}"
)


print("\nDescriptor summary:")

print(
    f"MW median:                   "
    f"{descriptor_summary['molecular_weight']['median']:.1f}"
)

print(
    f"MW range:                    "
    f"{descriptor_summary['molecular_weight']['min']:.1f}"
    f" - "
    f"{descriptor_summary['molecular_weight']['max']:.1f}"
)

print(
    f"Rotatable bonds median:      "
    f"{descriptor_summary['rotatable_bonds']['median']:.1f}"
)

print(
    f"Rotatable bonds maximum:     "
    f"{descriptor_summary['rotatable_bonds']['max']}"
)

print(
    f"Formal charge range:         "
    f"{descriptor_summary['formal_charge']['min']}"
    f" to "
    f"{descriptor_summary['formal_charge']['max']}"
)


print("\nNo compounds were removed.")

print("\nFile generated:")
print(f"  {OUTPUT_CSV}")
print(f"  {OUTPUT_METADATA}")
