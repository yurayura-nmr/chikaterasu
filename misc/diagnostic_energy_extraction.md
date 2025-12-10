# Chikaterasu — Optional Diagnostic Energy Extraction
This sidecar file collects optional `gmx energy` commands that can be activated
when additional thermodynamic or box-dimension diagnostics are needed.
These commands remain disabled by default to keep `chikaterasu.sh`
clean and focused on the standard MD workflow.

Typical use cases:
- verifying pressure stabilization during NPT
- checking anisotropic pressure tensor components
- inspecting box dimension fluctuations
- confirming temperature or density convergence
- debugging unexpected jumps in the simulation cell

To enable any diagnostic, simply uncomment the corresponding command block
in this file and copy it into your workflow (or uncomment the lines inside
`chikaterasu.sh` if you prefer inline diagnostics).


## Disabled-by-default diagnostics (copy/paste if needed)

```bash
# Temperature
#printf "Temperature" | gmx energy -f npt.edr -o ./temperature.xvg

# Density
#printf "Density" | gmx energy -f npt.edr -o ./density.xvg

# Pressure tensor components
#printf "Pres-XX" | gmx energy -f npt.edr -o ./pressure_XX.xvg
#printf "Pres-YY" | gmx energy -f npt.edr -o ./pressure_YY.xvg
#printf "Pres-ZZ" | gmx energy -f npt.edr -o ./pressure_ZZ.xvg
#printf "Pres-XY" | gmx energy -f npt.edr -o ./pressure_XY.xvg
#printf "Pres-XZ" | gmx energy -f npt.edr -o ./pressure_XZ.xvg
#printf "Pres-YX" | gmx energy -f npt.edr -o ./pressure_YX.xvg
#printf "Pres-YZ" | gmx energy -f npt.edr -o ./pressure_YZ.xvg
#printf "Pres-ZX" | gmx energy -f npt.edr -o ./pressure_ZX.xvg
#printf "Pres-ZY" | gmx energy -f npt.edr -o ./pressure_ZY.xvg

# Box vectors
#printf "Box-XX" | gmx energy -f npt.edr -o ./box_XX.xvg
#printf "Box-YY" | gmx energy -f npt.edr -o ./box_YY.xvg
#printf "Box-ZZ" | gmx energy -f npt.edr -o ./box_ZZ.xvg

# Volume
#printf "Volume" | gmx energy -f npt.edr -o ./volume.xvg
