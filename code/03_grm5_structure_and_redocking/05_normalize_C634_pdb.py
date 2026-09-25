# Normalize PDB records after restoring the native CYS634 residue in Chimera.
# Atomic coordinates are preserved; only residue annotations and record types
# associated with the former YCM modification are corrected.

input_file = "6FFI_complex_clean_C634.pdb"
output_file = "6FFI_complex_clean_C634_fixed.pdb"

with open(input_file) as fin, open(output_file, "w") as fout:
    for line in fin:

        # Correct sequence annotation: native residue is CYS634
        if line.startswith("SEQRES") and "YCM" in line:
            line = line.replace("YCM", "CYS")

        # Remove obsolete annotations for the modified residue YCM634
        if line.startswith("MODRES") and "YCM A  634" in line:
            continue

        if line.startswith("HET ") and "YCM  A 634" in line:
            continue

        if line.startswith("HETNAM") and "YCM" in line:
            continue

        if line.startswith("HETSYN") and "YCM" in line:
            continue

        if line.startswith("FORMUL") and "YCM" in line:
            continue

        # Peptide connectivity will be represented by standard protein records
        if line.startswith("LINK") and "YCM A 634" in line:
            continue

        # C634 is now a standard amino-acid residue
        if (
            line.startswith("HETATM")
            and line[17:20] == "CYS"
            and line[21] == "A"
            and line[22:26].strip() == "634"
        ):
            line = "ATOM  " + line[6:]

        fout.write(line)

print(f"Written: {output_file}")