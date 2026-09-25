from rdkit import Chem

supplier = Chem.SDMolSupplier(
    "MMPEP_redocking_dry.sdf",
    removeHs=False,
    sanitize=True
)

pose = supplier[0]

if pose is None:
    raise RuntimeError("Could not read docking mode 1.")

Chem.MolToPDBFile(
    pose,
    "MMPEP_redocking_dry_mode1.pdb"
)

print("Written: MMPEP_redocking_dry_mode1.pdb")
