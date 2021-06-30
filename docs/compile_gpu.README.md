# Compiling a GPU-enabled version of UniFrac

Note: The GPU-enabled version is currenlty only supported on Linux systems.
One can however run on Windows systems, too, using [CUDA-enabled WSL2](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).


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
conda create --name unifrac-gpu -c conda-forge -c bioconda python=3.6 unifrac
conda activate unifrac-gpu
conda install -c conda-forge -c bioconda gxx_linux-64=9.3
conda install -c conda-forge -c bioconda hdf5-static mkl-include
```

## Installing the NVIDIA HPC SDK

Currently, the only supported GPU-enabled compiler is the freely available [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk).

Note that internally the NVIDIA HPC SDK relies on GCC, which makes it possible for the resulting objects to link with the libraries provided through Anaconda. 

Our Anaconda environment provides GCC 9.3, but the executable names are mangled. In order to make it usable by the NVIDIA HPC SDK, we have to create a few symbolic links:

```
# Create GCC symbolic links
mkdir conda_nv_bins
(cd conda_nv_bins && for f in \
  ar as c++ cc cpp g++ gcc ld nm ranlib strip; \
  do \
    ln -s $CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-${f} ${f}; \
  done )

mkdir setup_scripts
echo "PATH=${PWD}/conda_nv_bins:\$PATH" \
  >> setup_scripts/setup_conda_nv_bins.sh
```

We are now ready to install the NVIDIA HPC SDK proper. Make sure you do this on a GPU-enabled node, as the installer checks for an existing CUDA driver.

Fell free to download the latest version froom the [NVIDIA official site](https://developer.nvidia.com/hpc-sdk). 

The following instructions will download and unpack version 21.5:

```
wget https://developer.download.nvidia.com/hpc-sdk/21.5/nvhpc_2021_215_Linux_x86_64_cuda_multi.tar.gz
tar xpzf nvhpc_2021_215_Linux_x86_64_cuda_multi.tar.gz
rm -f nvhpc_2021_215_Linux_x86_64_cuda_multi.tar.gz
```

Once you have the install directory unpacked, you need to patch the installer to use the right GCC version; then you are ready to run the actual installer:

```
source setup_scripts/setup_conda_nv_bins.sh

# must patch the  install scripts to find the right gcc
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/install 
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/*/*/compilers/bin/makelocalrc
sed -i -e "s#print_line 'set LOCALRC=YES;'#print_line 'set LOCALRC=YES;';print_line 'set DEFSTDOBJDIR=$CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/lib64;'#g" \
  nvhpc_*/install_components/*/*/compilers/bin/makelocalrc
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/install_cuda

export NVHPC_INSTALL_DIR=$PWD/hpc_sdk
export NVHPC_SILENT=true

(cd nvhpc_*; ./install)

# If you prefer an interactive install, do not export the above NVHPC env and
# Select "Single system install"
# Expand $PWD/hpc_sdk as install dir

cat > setup_nv.sh  << EOF
source $PWD/setup_scripts/setup_conda_nv_bins.sh

export PATH=\$PATH:`ls -d $PWD/hpc_sdk/*/202*/compilers/bin`
export ACC_CXX=pgc++

EOF
```

For convenience, we also create a setup script that will be used to properly setup the environment anytime needed.


## Compiling the GPU-enabled UniFrac with the NVIDIA HPC SDK

In order to compile UniFrac, you will need both the Anaconda dependencies, the NVIDIA HPC SDK and the UniFrac source code.
The first two you setup above, and can enable with the helper script; the later can be imported (once) with git.

```
source setup_nv.sh

# save original version of binaries
mkdir $CONDA_PREFIX/bin/org
mkdir $CONDA_PREFIX/lib/org

mv $CONDA_PREFIX/bin/ssu $CONDA_PREFIX/bin/org/
mv $CONDA_PREFIX/bin/faithpd $CONDA_PREFIX/bin/org/
mv $CONDA_PREFIX/lib/libssu*.so $CONDA_PREFIX/lib/org/

mv $CONDA_PREFIX/lib/python3.6/site-packages/unifrac $CONDA_PREFIX/lib/python3.6/site-packages/unifrac.org

git clone https://github.com/biocore/unifrac.git
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

