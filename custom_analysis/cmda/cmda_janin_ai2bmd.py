import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals

# DCD reader for ai2bmd

u = mda.Universe(
    "ub-preeq.pdb",
    "test.dcd"
)

# Universe already loaded as u
dt = u.trajectory.dt           # time between frames (ps)
n_frames = len(u.trajectory)  # number of frames
total_time = dt * n_frames     # total simulation time in ps

print(f"Timestep per frame: {dt} ps")
print(f"Number of frames: {n_frames}")
print(f"Total trajectory time: {total_time} ps")


# Residues of interest
roi1 = u.select_atoms("protein and resname ILE and resid 14")  # ILE 13
roi2 = u.select_atoms("protein and resname ILE and resid 24")  # ILE 23
roi3 = u.select_atoms("protein and resname ILE and resid 4")   # ILE 3
roi4 = u.select_atoms("protein and resname ILE and resid 31")  # ILE 30
roi5 = u.select_atoms("protein and resname ILE and resid 37")  # ILE 36
roi6 = u.select_atoms("protein and resname ILE and resid 45")  # ILE 44
roi7 = u.select_atoms("protein and resname ILE and resid 62")  # ILE 61

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
