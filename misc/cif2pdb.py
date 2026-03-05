import sys
import gemmi

if len(sys.argv) != 3:
    print("Usage: python cif_to_pdb.py input.cif output.pdb")
    sys.exit(1)

input_cif = sys.argv[1]
output_pdb = sys.argv[2]

# Read CIF structure
structure = gemmi.read_structure(input_cif)

# Write PDB
structure.write_pdb(output_pdb)

print(f"Converted {input_cif} -> {output_pdb}")
