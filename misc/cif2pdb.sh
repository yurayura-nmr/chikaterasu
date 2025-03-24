#!/bin/bash

# sudo apt-get install openbabel

# Check if the correct number of arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.cif> <output.pdb>"
    exit 1
fi

# Input and output file paths
input_cif="$1"
output_pdb="$2"

# Check if the input CIF file exists
if [ ! -f "$input_cif" ]; then
    echo "Error: CIF file '$input_cif' does not exist."
    exit 1
fi

# Convert CIF to PDB using OpenBabel
obabel "$input_cif" -O "$output_pdb"

# Check if the conversion was successful
if [ $? -eq 0 ]; then
    echo "Successfully converted '$input_cif' to '$output_pdb'."
else
    echo "Error during conversion."
    exit 1
fi

