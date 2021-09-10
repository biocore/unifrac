# Compiling a GPU-enabled version of UniFrac

Note: The GPU-enabled version is currenlty only supported on Linux systems.
One can however run on Windows systems, too, using [CUDA-enabled WSL2](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).


## System toools and libraries

While we will install most of the tools we need, there these instrcutions assume you have installed the following packages
```
wget make git glibc-devel
```

## Anaconda 

UniFrac has several dependencies, which we assume come via [Anaconda](https://www.anaconda.com/products/individual).

The instructions below has been tested with version [2021.05](https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh).

In case you have never used Anaconda below, here are the installation instruction:

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
chmod a+x Anaconda3-2021.05-Linux-x86_64.sh
./Anaconda3-2021.05-Linux-x86_64.sh
#log out and back in
```

## Create a dedicated environment

While it is possible to build a GPU-enabled UniFrac in any Anaconda environment, we assume a dedicated one in this document.
We call it **unifrac-gpu**.

Note: If you decide to change the used environment, you will have to make the appropriate changes to the scripts below. 

To create our **unifrac-gpu** with all the needed dependencies, run:

```
# create and activate unifrac-gpu Anaconda environment
conda create --name unifrac-gpu -c conda-forge -c bioconda unifrac
conda activate unifrac-gpu
conda install -c conda-forge -c bioconda gxx_linux-64=7.5.0 
conda install -c conda-forge -c bioconda hdf5-static mkl-include
```

## Installing the NVIDIA HPC SDK

Currently, the only supported GPU-enabled compiler is the freely available [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk).

Note that internally the NVIDIA HPC SDK relies on GCC, which makes it possible for the resulting objects to link with the libraries provided through Anaconda. 

We provide a helper script to install it, so you should fetch the code from github first, and then run that script.
You are encouraged to look into the script and improved as needed (e.g. updating the NVIDIA SDK version).
```
git clone https://github.com/biocore/unifrac.git
# the following will install and configure the NVIDIA SDK
./unifrac/scripts/install_hpc_sdk.sh
```

For convenience, we also create a setup script that will be used to properly setup the environment anytime needed.


## Compiling the GPU-enabled UniFrac with the NVIDIA HPC SDK

In order to compile UniFrac, you will need both the Anaconda dependencies, the NVIDIA HPC SDK and the UniFrac source code.

```
source setup_nv_h5.sh

# save original version of binaries
mkdir $CONDA_PREFIX/bin/org
mkdir $CONDA_PREFIX/lib/org

mv $CONDA_PREFIX/bin/ssu $CONDA_PREFIX/bin/org/
mv $CONDA_PREFIX/bin/faithpd $CONDA_PREFIX/bin/org/
mv $CONDA_PREFIX/lib/libssu*.so $CONDA_PREFIX/lib/org/

(cd $CONDA_PREFIX/lib/python*/site-packages && mv unifrac unifrac.org)

(cd unifrac/ && export USE_CYTHON=True && python setup.py build && python setup.py install)
```

And you are all done.
The UniFrac binary and libraries in the Anaconda environment are now the GPU-enabled ones.

*Note:* The produced binary has both GPU and CPU unifrac code paths. 
The GPU path will be used if a GPU is detected, and use the CPU code path else.
You can forecefully disable the GPU with:
```
export UNIFRAC_USE_GPU=N
```

If you want to know which path is being used, enabe info messages with:
```
export UNIFRAC_GPU_INFO=Y
```

## Compiling an older version of Unifrac for GPUs

If you do not want the cutting edge UniFrac from git, you will have to use version-specific instructions:
* [0.20.1 GPU compile instructions](https://github.com/sfiligoi/unifrac/blob/v0.20.1-docs/docs/compile_gpu.README.txt)
* [0.20.2 GPU compile instructions](https://github.com/sfiligoi/unifrac/blob/v0.20.2-docs/docs/compile_gpu.README.md)

