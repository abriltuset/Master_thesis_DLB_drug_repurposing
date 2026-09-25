from pathlib import Path
import statistics

pdb_file = Path("input/AF-P51693-F1-model_v6.pdb")

# Residues defining the two concordant fpocket cavities evaluated in APLP1:
# fpocket 42, corresponding to the P2Rank 2 / fpocket 42 pair
# fpocket 11, corresponding to the P2Rank 6 / fpocket 11 pair

pockets = {
    "pocket42": [49,106,108,109,110,111,136,139,141,145,146,147,148,150,170,171,174,191,192,193],
    "pocket11": [58,59,60,84,85,86,109,135,136,138]
}

plddt = {}

# In AlphaFold PDB files, per-residue pLDDT values are stored
# in the B-factor field. One value is extracted per CA atom.

with open(pdb_file) as f:
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            chain = line[21].strip()
            resid = int(line[22:26])
            score = float(line[60:66])

            if chain == "A":
                plddt[resid] = score


for pocket, residues in pockets.items():

    values = [plddt[r] for r in residues if r in plddt]

    print(f"\n=== {pocket.upper()} ===")

    for r in residues:
        print(f"{r}\t{plddt.get(r, 'NA')}")

    print(f"\nMean pLDDT: {statistics.mean(values):.2f}")
    print(f"Median pLDDT: {statistics.median(values):.2f}")
    print(f"Minimum pLDDT: {min(values):.2f}")
    print(f"Maximum pLDDT: {max(values):.2f}")
    print(f"Residues >=70: {sum(v >= 70 for v in values)}/{len(values)}")
    print(f"Residues <70: {sum(v < 70 for v in values)}/{len(values)}")
