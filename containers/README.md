
### Introduction

This README provides essential information for installing and
downloading Singularity containers used in the paper titled “Deep Whole
Genome Analysis of 494 Hepatocellular Carcinomas”.

The utilization of Singularity containers ensures the reproducibility of
research findings by encapsulating the entire software stack and
dependencies, facilitating sharing and deployment of the research
environment.

### Prerequisites

To install and utilize Singularity, a Linux system is required.

The disk space needed for Singularity is about 140Mb once it is
installed.

To leverage the full functionality of Singularity, the following kernel
support is required:

- **OverlayFS mounts** (minimum kernel \>=3.18): To ensure comprehensive
  flexibility in bind mounts to containers and to enable support for
  persistent overlays in writable containers.
- **Unprivileged user namespaces** (minimum kernel \>=3.8, \>=3.18
  recommended): To execute containers without the need for root or
  setuid privileges.

### Installation

Various methods are available for installing Singularity, with one
convenient approach being the installation from Linux distribution
repositories. Below, an example is provided that demonstrates the
installation process of Singularity specifically on CentOS/RHEL.

- Install the epel-release package and then install Singularity

``` bash
sudo yum update -y && \
    sudo yum install -y epel-release && \
    sudo yum update -y && \
    sudo yum install -y singularity-ce
```

- Check the Singularity version

``` bash
singularity --version
```

For further information on alternative methods of installing
Singularity, it is recommended to refer to the comprehensive guide
available at
<https://docs.sylabs.io/guides/3.0/user-guide/installation.html>.

### Downloading Singularity containers

The containers have been published on Zenodo at
<https://doi.org/10.5281/zenodo.7260221>.

To run the analysis using the containers, users are required to download
the `.sif` files and place them inside the `containers` folder located in
the top directory of a clone of this repository. This folder structure
is utilized by the analysis code to invoke the containers.

For a comprehensive list of software installed within the containers,
please refer to
[CLCA_containers_info.xlsx](https://github.com/ChongJenniferZhang/CLCA_WGS/raw/main/containers/CLCA_containers_info.xlsx).
