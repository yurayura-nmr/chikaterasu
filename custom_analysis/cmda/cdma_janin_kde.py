import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import numpy as np

# --- Load system ---
u = mda.Universe(
    "results/md_1/md_target.tpr",
    "results/md_1/md_target_centered_no_PBC.xtc"
)

# --- Residue numbers ---
ile_resids = [3, 13, 23, 30, 36, 44, 61]  
n_res = len(ile_resids)

# --- Grid layout 3x3 ---
n_rows, n_cols = 3, 3
fig, axes = plt.subplots(
    n_rows, n_cols,
    figsize=(12, 12),
    sharex=True,
    sharey=True
)
axes = axes.flatten()  # flatten to 1D array for easy indexing

# --- Loop over residues ---
for i, resid in enumerate(ile_resids):

    ax = axes[i]

    roi = u.select_atoms(f"protein and resname ILE and resid {resid}")

    if len(roi.residues) != 1:
        ax.set_title(f"ILE {resid} (selection error)")
        ax.axis("off")
        continue

    janin = dihedrals.Janin(roi)
    janin.run()
    angles = janin.results.angles

    chi1 = angles[:, 0, 0]
    chi2 = angles[:, 0, 1]

    # --- 2D KDE "mountain" style ---
    sns.kdeplot(
        x=chi1,
        y=chi2,
        fill=True,
        ax=ax,
        #cmap="viridis",
        thresh=0.05,
        levels=25,         # number of contour lines
        bw_method=0.3,     # adjust smoothness
        alpha=0.9
    )

    # --- Formatting ---
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 360)
    ax.set_aspect('equal')  # square plots
    ax.set_xlabel("χ1 (deg)")
    ax.set_ylabel("χ2 (deg)")
    ax.set_title(f"ILE {resid}", fontsize=12)

# --- Turn off unused subplots ---
for j in range(n_res, n_rows * n_cols):
    axes[j].axis("off")

fig.tight_layout()
fig.savefig("ILE_chi1_chi2_kde_3x3.png", dpi=300, bbox_inches="tight")
plt.show()

