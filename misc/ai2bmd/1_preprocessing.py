"""
Run this within pymol itself, not python.

conda install conda-forge::python
"""

cmd.load("1UBQ.pdb","molecule")
cmd.h_add("molecule")

cmd.wizard("mutagenesis")
cmd.get_wizard().set_n_cap("acet")
# Click N-terminal amino acid to select it
cmd.get_wizard().apply()

cmd.get_wizard().set_c_cap("nmet")
# Click C-terminal amino acid to select it
cmd.get_wizard().apply()

cmd.set_wizard()

