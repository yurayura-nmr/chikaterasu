import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals

u = mda.Universe(
    "../md_target.tpr",
    "../md_target_centered_no_PBC.xtc"
)

# Residues of interest
roi3 = u.select_atoms("protein and resname ILE and resid 3")   # ILE
roi1 = u.select_atoms("protein and resname ILE and resid 13")  # ILE
roi2 = u.select_atoms("protein and resname ILE and resid 23")  # ILE
roi4 = u.select_atoms("protein and resname ILE and resid 30")  # ILE
roi5 = u.select_atoms("protein and resname ILE and resid 36")  # ILE
roi6 = u.select_atoms("protein and resname ILE and resid 44")  # ILE
roi7 = u.select_atoms("protein and resname ILE and resid 61")  # ILE

rois = [
    ("ILE 13", roi1, "blue", "o"),
    ("ILE 23", roi2, "orange", "o"),
    ("ILE 3",  roi3, "green", "o"),
    ("ILE 30", roi4, "cyan", "o"),
    ("ILE 36", roi5, "black", "o"),
    ("ILE 44", roi6, "magenta", "o"),
    ("ILE 61", roi7, "grey", "o"),
]

fig, ax = plt.subplots(figsize=plt.figaspect(1))

"""
for label, roi, color, marker in rois:
    janin = dihedrals.Janin(roi, verbose=True).run()
    janin.plot(
        ax=ax,
        color=color,
        marker=marker,
        label=label,
        #ref=True
        ref=False
        #show_allowed=True
    )
"""

# --- Draw reference FIRST (background) ---
janin_ref = dihedrals.Janin(roi1).run()
janin_ref.plot(
    ax=ax,
    ref=True,
    color="0.7",      # light gray
    alpha=0.6
)

# --- Now draw data ON TOP ---
for label, roi, color, marker in rois:
    janin = dihedrals.Janin(roi).run()
    janin.plot(
        ax=ax,
        ref=False,     # important: don't redraw ref
        color=color,
        marker=marker,
        label=label,
        zorder=3       # force points to front
    )


ax.legend()
ax.set_title("Janin plot (Ile residues)")
plt.show()
plt.savefig("janin.pdf")
