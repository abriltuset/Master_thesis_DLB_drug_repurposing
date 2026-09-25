from rdkit import Chem

# Read crystallographic SDF used to define the box
mol = Chem.SDMolSupplier(
    "6ffi_J_D8B.sdf",
    removeHs=False,
    sanitize=True
)[0]

conf = mol.GetConformer()

xs, ys, zs = [], [], []

for atom in mol.GetAtoms():
    p = conf.GetAtomPosition(atom.GetIdx())
    xs.append(p.x)
    ys.append(p.y)
    zs.append(p.z)

center = (
    (min(xs) + max(xs)) / 2,
    (min(ys) + max(ys)) / 2,
    (min(zs) + max(zs)) / 2,
)

dimensions = (
    max(xs) - min(xs),
    max(ys) - min(ys),
    max(zs) - min(zs),
)

print("Ligand bounding-box center:")
print(f"x = {center[0]:.3f}")
print(f"y = {center[1]:.3f}")
print(f"z = {center[2]:.3f}")

print("\nLigand dimensions:")
print(f"x = {dimensions[0]:.3f}")
print(f"y = {dimensions[1]:.3f}")
print(f"z = {dimensions[2]:.3f}")
