# For use with ai2bmd
# https://github.com/microsoft/AI2BMD

from ase.io import read, write

# Load entire trajectory
traj = read('example.traj', index=':')

# Write to PDB format (suitable for GROMACS input)
write('example.pdb', traj)


# Better for pymol - convert to dcd format
# traj2dcd.py is found in ai2bmbd github /src/utils
python3 traj2dcd.py --input ../Logs-chig/chig-traj.traj --output test.dcd --pdb ../chig_preprocessed/chig-preeq.pdb --stride 1
