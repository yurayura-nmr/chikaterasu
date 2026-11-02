import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# 1. Map DSSP codes to numerical values, full secondary structure names, and colors
dssp_definitions = {
    '~': (0, 'Coil / Loop (~)', '#ffffff'),     # White
    ' ': (0, 'Coil / Loop (~)', '#ffffff'),
    'H': (1, r'$\alpha$-Helix (H)', '#ff0000'),   # Red
    'E': (2, r'$\beta$-Strand (E)', '#ffff00'),   # Yellow
    'G': (3, r'$3_{10}$-Helix (G)', '#ff00ff'),  # Magenta
    'I': (4, r'$\pi$-Helix (I)', '#800080'),     # Purple
    'T': (5, 'Turn (T)', '#0000ff'),           # Blue
    'S': (6, 'Bend (S)', '#00ffff'),           # Cyan
    'B': (7, r'$\beta$-Bridge (B)', '#00ff00'),  # Green
}

# Extract unique numerical codes, labels, and color palette
char_to_num = {k: v[0] for k, v in dssp_definitions.items()}
unique_entries = sorted(
    list({code: (label, color) for code, label, color in dssp_definitions.values()}.items())
)

num_codes = [e[0] for e in unique_entries]
labels = [e[1][0] for e in unique_entries]
palette = [e[1][1] for e in unique_entries]

# 2. Parse GROMACS dssp.dat file
with open('dssp.dat', 'r') as f:
    lines = [line.strip() for line in f if not line.startswith(('#', '@'))]

# Build 2D grid (rows = frames, cols = residues)
matrix = [[char_to_num.get(c, 0) for c in line] for line in lines]

# Transpose matrix so Y-axis = Residue index, X-axis = Trajectory Progress
data = np.array(matrix).T
num_residues, num_frames = data.shape

# 3. Render Plot
fig, ax = plt.subplots(figsize=(12, 6))

cmap = mcolors.ListedColormap(palette)
bounds = np.arange(-0.5, len(num_codes), 1)
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Set extent: X-axis spans 0% to 100%, Y-axis spans 1 to total residue count
img = ax.imshow(
    data,
    aspect='auto',
    origin='lower',
    cmap=cmap,
    norm=norm,
    extent=[0, 100, 1, num_residues]
)

# Axis Configuration
ax.set_xlabel('Trajectory Completion (%)', fontsize=12, fontweight='bold')
ax.set_ylabel('Residue Index', fontsize=12, fontweight='bold')
ax.set_title('Secondary Structure Timeline', fontsize=14, fontweight='bold')
ax.set_xlim(0, 100)
ax.set_ylim(1, num_residues)

# Custom colorbar mapped to DSSP names
cbar = fig.colorbar(img, ax=ax, ticks=num_codes, pad=0.02)
cbar.ax.set_yticklabels(labels, fontsize=10)
cbar.ax.tick_params(size=0)

# Export directly to vector PDF
plt.savefig('dssp_percentage_plot.pdf', dpi=300, bbox_inches='tight')
plt.show()