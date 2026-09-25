# Summarize mode 1 (best-scoring Vina pose) across the four independent runs.

from rdkit import Chem
import math
import csv
import statistics
import re

reference_file = "6ffi_J_D8B.sdf"

runs = [
    ("2026", "MMPEP_redocking_dry.sdf",
             "MMPEP_redocking_dry.pdbqt"),
    ("2027", "MMPEP_redocking_seed2027.sdf",
             "MMPEP_redocking_seed2027.pdbqt"),
    ("2028", "MMPEP_redocking_seed2028.sdf",
             "MMPEP_redocking_seed2028.pdbqt"),
    ("2029", "MMPEP_redocking_seed2029.sdf",
             "MMPEP_redocking_seed2029.pdbqt"),
]

output_csv = "MMPEP_redocking_replicates_summary.csv"


def fixed_frame_rmsd(ref, pose):
    """Symmetry-corrected heavy-atom RMSD without superposition."""

    ref = Chem.RemoveHs(ref)
    pose = Chem.RemoveHs(pose)

    if ref.GetNumAtoms() != pose.GetNumAtoms():
        raise ValueError(
            f"Different heavy-atom counts: "
            f"{ref.GetNumAtoms()} vs {pose.GetNumAtoms()}"
        )

    matches = pose.GetSubstructMatches(
        ref,
        uniquify=False,
        useChirality=False,
        maxMatches=1000
    )

    if not matches:
        raise ValueError("No atom mapping found.")

    ref_conf = ref.GetConformer()
    pose_conf = pose.GetConformer()

    best_rmsd = float("inf")

    for match in matches:
        d2_sum = 0.0

        for ref_idx, pose_idx in enumerate(match):
            p1 = ref_conf.GetAtomPosition(ref_idx)
            p2 = pose_conf.GetAtomPosition(pose_idx)

            d2_sum += (
                (p1.x - p2.x)**2 +
                (p1.y - p2.y)**2 +
                (p1.z - p2.z)**2
            )

        rmsd = math.sqrt(d2_sum / len(match))
        best_rmsd = min(best_rmsd, rmsd)

    return best_rmsd


def first_affinity(pdbqt_file):
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])

    raise RuntimeError(f"No Vina score found in {pdbqt_file}")


# Reference crystallographic ligand
reference = Chem.SDMolSupplier(
    reference_file,
    removeHs=False,
    sanitize=True
)[0]

if reference is None:
    raise RuntimeError("Could not read reference ligand.")


results = []

for seed, sdf_file, pdbqt_file in runs:

    supplier = Chem.SDMolSupplier(
        sdf_file,
        removeHs=False,
        sanitize=True
    )

    # Only mode 1: best-scoring Vina pose
    pose1 = supplier[0]

    if pose1 is None:
        raise RuntimeError(f"Could not read mode 1 from {sdf_file}")

    rmsd = fixed_frame_rmsd(reference, pose1)
    affinity = first_affinity(pdbqt_file)

    results.append((seed, affinity, rmsd))

    print(
        f"Seed {seed} | "
        f"Affinity = {affinity:.3f} kcal/mol | "
        f"RMSD = {rmsd:.3f} Å"
    )


# Summary statistics
rmsds = [x[2] for x in results]
affinities = [x[1] for x in results]

print("\nSummary")
print("-------")
print(
    f"RMSD mean ± SD = "
    f"{statistics.mean(rmsds):.3f} ± "
    f"{statistics.stdev(rmsds):.3f} Å"
)
print(
    f"Affinity mean ± SD = "
    f"{statistics.mean(affinities):.3f} ± "
    f"{statistics.stdev(affinities):.3f} kcal/mol"
)

successful = sum(r <= 2.0 for r in rmsds)

print(
    f"Successful pose recoveries (RMSD <= 2 Å): "
    f"{successful}/{len(rmsds)}"
)


# CSV output
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "seed",
        "affinity_kcal_mol",
        "mode1_fixed_frame_symmetry_corrected_RMSD_A"
    ])
    writer.writerows(results)

print(f"\nWritten: {output_csv}")
