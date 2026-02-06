import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import numpy as np

# --- Load system ---
u = mda.Universe(
    "../md_target.tpr",
    "../md_target_centered_no_PBC.xtc"
)

# >>> PUT YOUR RESID NUMBERS HERE <<<
ile_resids = [3, 13, 23, 30, 36, 44, 61]   # replace zeros later
n_res = len(ile_resids)

# --- Time axis (shared for all plots) ---
# Will be rebuilt per residue, but this defines dt cleanly
dt_ns = u.trajectory.dt / 1000.0

# --- Figure ---
fig, axes = plt.subplots(
    n_res, 1,
    figsize=(10, 2.2 * n_res),
    sharex=True
)

# If only one subplot, make it iterable
if n_res == 1:
    axes = [axes]

# --- Loop over residues ---
for ax, resid in zip(axes, ile_resids):

    roi = u.select_atoms(f"protein and resname ILE and resid {resid}")

    if len(roi.residues) != 1:
        ax.set_title(f"ILE {resid} (selection error)")
        continue

    janin = dihedrals.Janin(roi)
    janin.run()

    angles = janin.results.angles
    chi1 = angles[:, 0, 0]
    chi2 = angles[:, 0, 1]

    # Time axis (must match analysis length)
    time_ns = np.linspace(
        0.0,
        (len(chi1) - 1) * dt_ns,
        len(chi1)
    )

    ax.plot(time_ns, chi1, label="χ1", lw=1.1)
    ax.plot(time_ns, chi2, label="χ2", lw=1.1)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 360)
    ax.set_ylabel("deg")
    ax.set_title(f"ILE {resid}")

    ax.legend(loc="upper right", fontsize=9)

# --- Common labels ---
axes[-1].set_xlabel("Time (ns)")

fig.tight_layout()
fig.savefig("ILE_chi_vs_time_all.png", dpi=300, bbox_inches="tight")
plt.show()

