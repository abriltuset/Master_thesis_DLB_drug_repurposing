from rdkit import Chem
import math
import csv

reference_file = "6ffi_J_D8B.sdf"
poses_file = "MMPEP_redocking_dry.sdf"
pdbqt_file = "MMPEP_redocking_dry.pdbqt"
output_csv = "MMPEP_redocking_dry_rmsd.csv"


def fixed_frame_rmsd(ref, pose):
    """
    Symmetry-corrected heavy-atom RMSD WITHOUT translating,
    rotating or superimposing the molecules.

    Both molecules must already be in the same receptor
    coordinate frame.
    """

    ref = Chem.RemoveHs(ref)
    pose = Chem.RemoveHs(pose)

    if ref.GetNumAtoms() != pose.GetNumAtoms():
        raise ValueError(
            f"Different heavy-atom counts: "
            f"{ref.GetNumAtoms()} vs {pose.GetNumAtoms()}"
        )

    # Search all graph-equivalent atom mappings to account for molecular symmetry.
    # RMSD is calculated in the fixed receptor coordinate frame; no structural
    # alignment or superposition is performed.
    matches = pose.GetSubstructMatches(
        ref,
        uniquify=False,
        useChirality=False,
        maxMatches=1000
    )

    if not matches:
        raise ValueError("No atom mapping found between reference and pose.")

    ref_conf = ref.GetConformer()
    pose_conf = pose.GetConformer()

    best_rmsd = float("inf")

    for match in matches:
        squared_distances = []

        for ref_idx, pose_idx in enumerate(match):
            p_ref = ref_conf.GetAtomPosition(ref_idx)
            p_pose = pose_conf.GetAtomPosition(pose_idx)

            d2 = (
                (p_ref.x - p_pose.x) ** 2
                + (p_ref.y - p_pose.y) ** 2
                + (p_ref.z - p_pose.z) ** 2
            )

            squared_distances.append(d2)

        rmsd = math.sqrt(
            sum(squared_distances) / len(squared_distances)
        )

        best_rmsd = min(best_rmsd, rmsd)

    return best_rmsd


# ------------------------------------------------------------
# Read crystallographic reference
# ------------------------------------------------------------

ref_supplier = Chem.SDMolSupplier(
    reference_file,
    removeHs=False,
    sanitize=True
)

reference = ref_supplier[0]

if reference is None:
    raise RuntimeError("Could not read crystallographic M-MPEP.")


# ------------------------------------------------------------
# Read Vina affinities from PDBQT
# ------------------------------------------------------------

affinities = []

with open(pdbqt_file) as f:
    for line in f:
        if line.startswith("REMARK VINA RESULT:"):
            affinities.append(float(line.split()[3]))


# ------------------------------------------------------------
# Read docked poses
# ------------------------------------------------------------

pose_supplier = Chem.SDMolSupplier(
    poses_file,
    removeHs=False,
    sanitize=True
)

results = []

for mode, pose in enumerate(pose_supplier, start=1):

    if pose is None:
        print(f"Mode {mode}: could not be read")
        continue

    rmsd = fixed_frame_rmsd(reference, pose)

    affinity = (
        affinities[mode - 1]
        if mode <= len(affinities)
        else None
    )

    results.append((mode, affinity, rmsd))

    print(
        f"Mode {mode:2d} | "
        f"Affinity = {affinity:7.3f} kcal/mol | "
        f"RMSD = {rmsd:6.3f} Å"
    )


# ------------------------------------------------------------
# Save table
# ------------------------------------------------------------

with open(output_csv, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "mode",
        "affinity_kcal_mol",
        "fixed_frame_symmetry_corrected_RMSD_A"
    ])

    writer.writerows(results)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

if results:

    best_rmsd = min(results, key=lambda x: x[2])

    print("\nBest crystallographic-pose recovery:")
    print(
        f"Mode {best_rmsd[0]} | "
        f"Affinity = {best_rmsd[1]:.3f} kcal/mol | "
        f"RMSD = {best_rmsd[2]:.3f} Å"
    )

    print(f"\nWritten: {output_csv}")
