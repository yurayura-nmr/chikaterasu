# MDAnalysis + NGLView Quick Visualization Tutorial

This short guide shows how to quickly launch a Jupyter Lab session and visualize an MD trajectory using **MDAnalysis** and **nglview**.  
Useful for rapid trajectory inspection and generating interactive molecular views.

---

## 1. Launching Jupyter Lab

From your trajectory directory:

```bash
# Install required dependencies
pip install mdanalysis nglview numpy

cd results/md_1/traj
jupyter lab
````

Then open the link printed in the terminal (usually `http://localhost:8888/lab`).

---

## 2. Minimal Notebook to Visualize Your Trajectory

Create a new notebook in Jupyter Lab and run the following:

```python
import MDAnalysis as mda, nglview as nv, numpy as np

u = mda.Universe("../md_target.tpr", "../md_target_centered_no_PBC.xtc")

# Load trajectory into memory (optional but speeds up playback)
u.transfer_to_memory(step=5)

# Create viewer
view = nv.show_mdanalysis(u)

# Add unit cell display
view.add_unitcell()

# Apply rotation (example from MDAnalysis tutorial)
view.control.rotate(mda.lib.transformations.quaternion_from_euler(-np.pi/2, np.pi/3, np.pi/6, 'rzyz').tolist())

# Zoom out slightly
view.control.zoom(-0.3)

# Display widget (must be last line of cell)
view
```

This will open an interactive NGL viewer inside the notebook, showing your trajectory with PBC box and rotation applied.

---

## 3. Minimal Notebook to Visualize Your Trajectory with Water Displayed

```python
import MDAnalysis as mda, nglview as nv, numpy as np

u = mda.Universe("../md.tpr", "../md_full.xtc")

# Load trajectory into memory (optional but speeds up playback)
u.transfer_to_memory(step=10)

# Create viewer
view = nv.show_mdanalysis(u)

# Add unit cell display
view.add_unitcell()

# Apply rotation (example from MDAnalysis tutorial)
view.control.rotate(mda.lib.transformations.quaternion_from_euler(-np.pi/2, np.pi/3, np.pi/6, 'rzyz').tolist())

# Show water atoms
view.add_representation("spacefill", selection="water", radius=0.2, opacity=0.3, color="skyblue")

# Zoom out slightly
view.control.zoom(-0.3)

# Display widget (must be last line of cell)
view
```

This will open an interactive NGL viewer inside the notebook, showing your trajectory with PBC box and rotation applied.
