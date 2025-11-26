#!/bin/bash

# PC Setup Script for chikaterasu MD Analysis Repository [untested]
# This script automates the initial setup of a new PC for molecular dynamics analysis

set -e  # Exit on any error

echo "Starting chikaterasu PC setup..."

# ---------------------------------------------------------------------------
# System Configuration Checks
# ---------------------------------------------------------------------------
echo "1. Checking system configuration..."
nvidia-smi  # Confirm GPU is correctly recognized (driver working?)

# Note: Manual steps required:
# - In PopOS GUI turn off display dimming
# - Also turn off auto suspend

# ---------------------------------------------------------------------------
# System Package Updates
# ---------------------------------------------------------------------------
echo "2. Updating system packages..."
sudo apt-get update
sudo apt-get upgrade -y

# ---------------------------------------------------------------------------
# Essential Package Installation
# ---------------------------------------------------------------------------
echo "3. Installing essential packages..."
sudo apt install -y plocate
sudo apt-get install -y ssh
sudo apt-get install -y libssh-dev  # required for cmake
sudo apt-get install -y htop
sudo apt-get install -y libfftw3-dev # can also build from scratch 
sudo apt-get install -y cmake-curses-gui
sudo apt install -y libflame-dev # may be needed on some hardware 
sudo apt-get install -y libopenmpi-dev
sudo apt-get install -y vim
sudo apt-get install -y cmake
sudo apt-get install -y wget
sudo apt-get install -y grace # displaying gromacs output xvg

# ---------------------------------------------------------------------------
# CUDA Toolkit Installation
# ---------------------------------------------------------------------------
echo "4. Installing CUDA toolkit..."
# Remove existing installation if needed
# sudo apt-get remove nvidia-cuda-toolkit

cd /tmp
wget https://developer.download.nvidia.com/compute/cuda/12.8.0/local_installers/cuda_12.8.0_570.86.10_linux.run
chmod +x cuda_12.8.0_570.86.10_linux.run
sudo ./cuda_12.8.0_570.86.10_linux.run --toolkit --silent --override

# Add CUDA to PATH (temporary for this session)
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

# ---------------------------------------------------------------------------
# GROMACS Installation
# ---------------------------------------------------------------------------
echo "5. Building and installing GROMACS..."

# Note: GROMACS 2025.3 requires CMake 4.0 or newer
# This was confirmed to work with the Chikaterasu repository
# If CMake is not already installed, you can easily install it with:
#   wget https://github.com/Kitware/CMake/releases/download/v4.2.0/cmake-4.2.0.tar.gz
#   tar xvf cmake-4.2.0.tar.gz
#   cd cmake-4.2.0
#   ./configure
#   make -j$(nproc)
#   sudo make install

mkdir -p ~/chikaterasu/setup
cd ~/chikaterasu/setup

# Download and build GROMACS
# Note: Use the highest/newest GROMACS version compatible with your hardware.
# However, newer GROMACS versions often require a recent CMake version that may
# not be available via system package managers (apt-get). If the build fails,
# you may need to manually build and install a newer CMake version first.
wget https://ftp.gromacs.org/gromacs/gromacs-2025.3.tar.gz
tar xvf gromacs-2025.3.tar.gz
cd gromacs-2025.3

mkdir build
cd build

cmake .. \
    -DREGRESSIONTEST_DOWNLOAD=ON \
    -DGMX_GPU=CUDA \
    -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
    -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs

# Note: If the downloaded regression test version differs significantly from your
# GROMACS version, some tests may fail. In that case, manually download the 
# regression tests that match your GROMACS version, extract the tarball, and use:
#
# cmake .. \
#     -DREGRESSIONTEST_DOWNLOAD=OFF \
#     -DREGRESSIONTEST_PATH=/path/to/regressiontests-xxxx \
#     -DGMX_GPU=CUDA \
#     -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
#     -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs

make -j$(nproc)
make check -j$(nproc)
sudo make install

# ---------------------------------------------------------------------------
# Environment Configuration
# ---------------------------------------------------------------------------
echo "6. Configuring environment..."

# Add to bashrc for permanent setup
echo "Adding GROMACS and CUDA to environment..."

cat >> ~/.bashrc << 'EOF'

# chikaterasu MD Analysis Environment
export PATH=/usr/local/gromacs/bin:$PATH
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
export CUDA_HOME=/usr/local/cuda

EOF

# Source bashrc to apply changes immediately
source ~/.bashrc

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
echo "7. Verifying installation..."

echo "CUDA version:"
nvcc --version || echo "CUDA not in PATH, may need to restart terminal"

echo "GROMACS version:"
gmx --version || echo "GROMACS not in PATH, may need to restart terminal"

echo ""
echo "Setup complete! Notes:"
echo "- Manual steps: Disable display dimming and auto-suspend in PopOS GUI"
echo "- Test ubiquitin with chikaterasu 5 to verify GPU functionality"
echo "- You may need to restart your terminal for all paths to take effect"
