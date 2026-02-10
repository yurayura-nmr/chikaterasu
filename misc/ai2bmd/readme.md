# 🐳 Installing Docker + NVIDIA Container Toolkit (Pop!_OS / Ubuntu)

```bash
sudo apt update
sudo apt install docker.io -y
```

### Start and enable Docker:

```bash
sudo systemctl enable docker
sudo systemctl start docker
```

### (Optional) Add your user to the Docker group:

```bash
sudo usermod -aG docker $USER
newgrp docker
```

---

Make sure your system has the NVIDIA driver installed:

```bash
nvidia-smi
```

If not installed, do:

```bash
sudo apt install nvidia-driver-535  # Replace with the latest version if needed
sudo reboot
```

---

## Install NVIDIA Container Toolkit

### Add the NVIDIA repository:

```bash
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/libnvidia-container/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
```

### Install the toolkit:

```bash
sudo apt update
sudo apt install -y nvidia-container-toolkit
```

### Configure Docker to use the NVIDIA runtime:

```bash
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

---

```bash
docker run --rm --gpus all nvidia/cuda:12.2.0-base-ubuntu20.04 nvidia-smi
```

To test the installation, run the example from: https://github.com/microsoft/AI2BMD
