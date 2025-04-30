# 🐳 Installing Docker + NVIDIA Container Toolkit (Pop!_OS / Ubuntu)

This guide provides steps to install Docker along with the NVIDIA Container Toolkit to enable GPU-accelerated containers (e.g., for molecular dynamics simulations).

---

## 📦 1. Install Docker Engine

### Update packages and install Docker:

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

## 🖥️ 2. Verify NVIDIA Driver

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

## ⚙️ 3. Install NVIDIA Container Toolkit

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

## ✅ 4. Test GPU Access in Docker

```bash
docker run --rm --gpus all nvidia/cuda:12.2.0-base-ubuntu20.04 nvidia-smi
```

You should see your GPU listed just like on the host system.

