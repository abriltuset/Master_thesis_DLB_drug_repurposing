#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


# Main NAM-pocket residues using the 6FFI construct numbering
KNOWN_POCKET = {
    ("ILE", "625", "A"),
    ("PRO", "655", "A"),
    ("TYR", "659", "A"),
    ("TRP", "945", "A"),   # W785 en UniProt
    ("PHE", "948", "A"),   # F788 en UniProt
    ("VAL", "966", "A"),   # V806 en UniProt
    ("SER", "969", "A"),   # S809 en UniProt
    ("ALA", "970", "A"),   # A810 en UniProt
}


def distance(a, b):
    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


def read_receptor(pdb_file):
    """Read protein heavy atoms from receptor PDB."""
    atoms = []

    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            element = line[76:78].strip()

            if not element:
                element = atom_name[0]

            if element.upper() == "H":
                continue

            atoms.append(
                {
                    "atom": atom_name,
                    "resname": resname,
                    "resnum": resnum,
                    "chain": chain,
                    "xyz": (x, y, z),
                }
            )

    return atoms


def read_ligand(pdbqt_file):
    """Read ligand heavy atoms and Vina score from PDBQT."""
    atoms = []
    vina_score = None

    with open(pdbqt_file) as f:
        for line in f:

            if line.startswith("REMARK VINA RESULT"):
                try:
                    vina_score = float(line.split()[3])
                except (IndexError, ValueError):
                    pass

            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom_name = line[12:16].strip()

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            # AutoDock atom type is last field
            parts = line.split()
            autodock_type = parts[-1]

            # Exclude explicit hydrogens
            if autodock_type.upper().startswith("H"):
                continue

            atoms.append(
                {
                    "atom": atom_name,
                    "type": autodock_type,
                    "xyz": (x, y, z),
                }
            )

    return atoms, vina_score


def analyse_pose(receptor_atoms, ligand_file, contact_cutoff, clash_cutoff):

    ligand_atoms, vina_score = read_ligand(ligand_file)

    residue_min_dist = defaultdict(lambda: float("inf"))
    short_contacts = []

    global_min_dist = float("inf")

    for p in receptor_atoms:
        residue_key = (
            p["resname"],
            p["resnum"],
            p["chain"],
        )

        for l in ligand_atoms:

            d = distance(p["xyz"], l["xyz"])

            if d < global_min_dist:
                global_min_dist = d

            if d < residue_min_dist[residue_key]:
                residue_min_dist[residue_key] = d

            if d < clash_cutoff:
                short_contacts.append(
                    {
                        "distance": d,
                        "resname": p["resname"],
                        "resnum": p["resnum"],
                        "chain": p["chain"],
                        "protein_atom": p["atom"],
                        "ligand_atom": l["atom"],
                    }
                )

    contacting_residues = {
        residue: d
        for residue, d in residue_min_dist.items()
        if d <= contact_cutoff
    }

    known_contacts = {
        residue: d
        for residue, d in contacting_residues.items()
        if residue in KNOWN_POCKET
    }

    return {
        "ligand": Path(ligand_file).stem,
        "vina_score": vina_score,
        "n_ligand_heavy_atoms": len(ligand_atoms),
        "global_min_dist": global_min_dist,
        "contacts": contacting_residues,
        "known_contacts": known_contacts,
        "short_contacts": short_contacts,
    }


def residue_label(residue):
    resname, resnum, chain = residue
    return f"{resname}{resnum}.{chain}"


def main():

    parser = argparse.ArgumentParser(
        description="Analyse protein-ligand contacts for Vina poses."
    )

    parser.add_argument(
        "receptor",
        help="Receptor PDB file"
    )

    parser.add_argument(
        "ligands",
        nargs="+",
        help="One or more Vina output PDBQT files (single pose)"
    )

    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=4.0,
        help="Contact distance cutoff in Å (default: 4.0)"
    )

    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=2.0,
        help="Very short heavy-atom contact cutoff in Å (default: 2.0)"
    )

    parser.add_argument(
        "--csv",
        default="pose_review_summary.csv",
        help="Output summary CSV"
    )

    args = parser.parse_args()

    receptor_atoms = read_receptor(args.receptor)

    print(f"\nProtein heavy atoms loaded: {len(receptor_atoms)}")
    print(f"Contact cutoff: {args.contact_cutoff:.1f} Å")
    print(f"Short-contact cutoff: {args.clash_cutoff:.1f} Å")

    results = []

    for ligand in args.ligands:

        result = analyse_pose(
            receptor_atoms,
            ligand,
            args.contact_cutoff,
            args.clash_cutoff,
        )

        results.append(result)

        print("\n" + "=" * 70)
        print(result["ligand"])
        print("=" * 70)

        print(f"Vina score: {result['vina_score']}")
        print(f"Ligand heavy atoms: {result['n_ligand_heavy_atoms']}")
        print(f"Minimum protein-ligand distance: "
              f"{result['global_min_dist']:.2f} Å")

        print(
            f"\nResidues within {args.contact_cutoff:.1f} Å "
            f"({len(result['contacts'])}):"
        )

        for residue, d in sorted(
            result["contacts"].items(),
            key=lambda x: x[1]
        ):
            marker = " *NAM pocket*" if residue in KNOWN_POCKET else ""
            print(
                f"  {residue_label(residue):12s}"
                f" {d:5.2f} Å{marker}"
            )

        print(
            f"\nKnown NAM-pocket residues contacted: "
            f"{len(result['known_contacts'])}/{len(KNOWN_POCKET)}"
        )

        if result["known_contacts"]:
            print(
                "  "
                + ", ".join(
                    residue_label(r)
                    for r in sorted(
                        result["known_contacts"],
                        key=lambda x: int(x[1])
                    )
                )
            )

        print(
            f"\nVery short heavy-atom contacts "
            f"(<{args.clash_cutoff:.1f} Å):"
        )

        if result["short_contacts"]:
            for c in sorted(
                result["short_contacts"],
                key=lambda x: x["distance"]
            ):
                print(
                    f"  {c['distance']:.2f} Å   "
                    f"{c['resname']}{c['resnum']}.{c['chain']} "
                    f"{c['protein_atom']} -- ligand {c['ligand_atom']}"
                )
        else:
            print("  None")

    # Summary CSV
    with open(args.csv, "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow([
            "ligand",
            "vina_score_kcal_mol",
            "ligand_heavy_atoms",
            "min_protein_ligand_distance_A",
            "n_residues_within_contact_cutoff",
            "n_known_NAM_residues_contacted",
            "known_NAM_residues_contacted",
            "n_short_contacts_within_clash_cutoff",
        ])

        for r in results:

            writer.writerow([
                r["ligand"],
                r["vina_score"],
                r["n_ligand_heavy_atoms"],
                round(r["global_min_dist"], 3),
                len(r["contacts"]),
                len(r["known_contacts"]),
                "; ".join(
                    residue_label(res)
                    for res in sorted(
                        r["known_contacts"],
                        key=lambda x: int(x[1])
                    )
                ),
                len(r["short_contacts"]),
            ])

    print("\n" + "=" * 70)
    print(f"Summary written to: {args.csv}")
    print("=" * 70)


if __name__ == "__main__":
    main()
