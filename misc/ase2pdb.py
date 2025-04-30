# For use with ai2bmd
# https://github.com/microsoft/AI2BMD

from ase.io import read, write

# Load entire trajectory
traj = read('example.traj', index=':')

# Write to PDB format (suitable for GROMACS input)
write('example.pdb', traj)
