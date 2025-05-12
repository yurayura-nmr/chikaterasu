import numpy as np
import matplotlib.pyplot as plt

def read_xvg(filename):
    """Reads a GROMACS .xvg file and returns the data (excluding comments and metadata)."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    time = float(parts[0])
                    value = float(parts[1])
                    data.append((time, value))
                except ValueError:
                    continue
    return np.array(data)

# Replace with your actual file paths
file1 = 'results/md_1/distance/distance_1.xvg'
file2 = '/media/kentayamada/data/Q40E_1000ns/results/md_1/distance/distance_1.xvg'

data1 = read_xvg(file1)
data2 = read_xvg(file2)

# Extract only the distance values (second column)
distances1 = data1[:, 1]
distances2 = data2[:, 1]

print(distances1)
print(distances2)

def overlay():
    # Plot overlaid histograms
    plt.figure(figsize=(10, 6))
    plt.hist(distances1, bins=50, alpha=0.8, label='Distance 1', color='black', density=True)
    plt.hist(distances2, bins=50, alpha=0.8, label='Distance 2', color='red', density=True)

    plt.xlabel('Distance (nm)')
    plt.ylabel('Probability Density')
    plt.title('Overlaid Histograms of Distances from GROMACS')
    #plt.legend()
    #plt.grid(True)
    #plt.tight_layout()
    plt.show()

def sidebyside():
    # Create side-by-side plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Histogram for distance 1
    axes[0].hist(distances1, bins=50, color='blue', alpha=0.7, density=True)
    axes[0].set_title('Distance 1')
    axes[0].set_xlabel('Distance (nm)')
    axes[0].set_ylabel('Probability Density')
    axes[0].grid(True)

    # Histogram for distance 2
    axes[1].hist(distances2, bins=50, color='orange', alpha=0.7, density=True)
    axes[1].set_title('Distance 2')
    axes[1].set_xlabel('Distance (nm)')
    axes[1].grid(True)

    plt.suptitle('Side-by-Side Histograms of GROMACS Distances')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

sidebyside()

