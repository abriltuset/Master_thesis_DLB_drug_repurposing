from rdkit import Chem

input_file = "6ffi_J_D8B.sdf"
output_file = "MMPEP_H.sdf"

# Read crystallographic ligand
supplier = Chem.SDMolSupplier(
    input_file,
    removeHs=False,
    sanitize=True
)

mol = supplier[0]

if mol is None:
    raise RuntimeError("RDKit could not read the M-MPEP SDF.")

print("Heavy atoms before H addition:", mol.GetNumHeavyAtoms())
print("Total atoms before H addition:", mol.GetNumAtoms())
print("Formal charge:", Chem.GetFormalCharge(mol))

# Store experimental heavy-atom coordinates
conf_before = mol.GetConformer()
coords_before = {
    atom.GetIdx(): conf_before.GetAtomPosition(atom.GetIdx())
    for atom in mol.GetAtoms()
    if atom.GetAtomicNum() > 1
}

# Add hydrogens while keeping existing 3D coordinates
mol_H = Chem.AddHs(mol, addCoords=True)

# Verify that heavy atoms have not moved
conf_after = mol_H.GetConformer()

max_displacement = 0.0

for idx, p_before in coords_before.items():
    p_after = conf_after.GetAtomPosition(idx)

    displacement = (
        (p_before.x - p_after.x)**2 +
        (p_before.y - p_after.y)**2 +
        (p_before.z - p_after.z)**2
    )**0.5

    max_displacement = max(max_displacement, displacement)

print("Total atoms after H addition:", mol_H.GetNumAtoms())
print(
    "Maximum heavy-atom displacement:",
    f"{max_displacement:.6f} Å"
)

# Write ligand with explicit hydrogens
writer = Chem.SDWriter(output_file)
writer.write(mol_H)
writer.close()

print("Written:", output_file)


import subprocess

# Convert the prepared crystallographic ligand to PDBQT using Meeko
subprocess.run(
    [
        "mk_prepare_ligand.py",
        "-i", output_file,
        "-o", "MMPEP.pdbqt"
    ],
    check=True
)

print("Written: MMPEP.pdbqt")