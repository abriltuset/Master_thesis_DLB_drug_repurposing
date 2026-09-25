# "Consensus" regions correspond to residues shared between the
# fpocket and P2Rank predictions and are used as the structural
# core of each concordant pocket.

# AlphaFold PAE is directional. For each residue pair, the two
# directional PAE values are averaged to obtain a symmetric summary.

import json
import statistics
from pathlib import Path

pae_file = Path("input/AF-P51693-F1-predicted_aligned_error_v6.json")

regions = {
    "pocket42_full": [
        49,106,108,109,110,111,136,139,141,145,
        146,147,148,150,170,171,174,191,192,193
    ],

    # 15 residues shared exactly with P2Rank pocket 2
    "pocket42_consensus": [
        106,108,109,111,141,145,146,147,148,
        170,171,174,191,192,193
    ],

    "pocket11_full": [
        58,59,60,84,85,86,109,135,136,138
    ],

    # Residues shared with P2Rank pocket 6
    "pocket11_consensus": [
        58,59,85,135,136,138
    ]
}

with open(pae_file) as f:
    data = json.load(f)

if isinstance(data, list):
    pae = data[0]["predicted_aligned_error"]
else:
    pae = data["predicted_aligned_error"]

for name, residues in regions.items():

    values = []
    problematic = []

    for i, r1 in enumerate(residues):
        for r2 in residues[i+1:]:

            p12 = pae[r1 - 1][r2 - 1]
            p21 = pae[r2 - 1][r1 - 1]

            sym = (p12 + p21) / 2

            values.append(sym)

            if sym > 10:
                problematic.append((r1, r2, sym))

    print(f"\n=== {name.upper()} ===")
    print(f"Residues: {len(residues)}")
    print(f"Pairs: {len(values)}")
    print(f"Mean symmetric PAE: {statistics.mean(values):.2f} Å")
    print(f"Median symmetric PAE: {statistics.median(values):.2f} Å")
    print(f"Maximum symmetric PAE: {max(values):.2f} Å")
    print(f"Pairs >10 Å: {len(problematic)}/{len(values)}")

    if problematic:
        print("\nPairs with symmetric PAE >10 Å:")
        for r1, r2, value in problematic:
            print(f"{r1}-{r2}: {value:.2f} Å")
