# Compiling a GPU-enabled version of UniFrac

Note: The GPU-enabled version is currenlty only supported on Linux systems.
One can however run on Windoows systems, too, using [CUDA-enabled WSL2](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).


## Anaconda 

UniFrac has several dependencies, which we assume come via [Anaconda](https://www.anaconda.com/products/individual).

The instructions below has been tested with version [2020.07](https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh).

In case you have never used Anaconda below, here are the installation instruction:

```
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
chmod a+x Anaconda3-2020.07-Linux-x86_64.sh
./Anaconda3-2020.07-Linux-x86_64.sh
#log out and back in
```

## Create a dedicated environment

While it is possible to build a GPU-enabled UniFrac in any Anaconda environment, we asume a dedicated one in this document.
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

## Installing the NVIDIA HPC DK

Currently, the only supported GPU-enabled compiler is the freely available [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk).

Note that internally the NVIDIA HPC SDK relies on GCC, which makes it possible for the resulting objects to link with the libraries provided through Anaconda. 

Our Anaconda environment provides GCC 7.5, but the executable names are mangled. In order to make it usable by the NVIDIA HPC SDK, we have to create a few symbolic links:

```
# Create GCC symbolic links
mkdir conda_nv_bins
(cd conda_nv_bins && for f in \
  ar as c++ cc cpp g++ gcc ld nm ranlib strip; \
  do \
    ln -s $CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-${f} ${f}; \
  done )

mkdir setup_scripts
echo "conda activate unifrac-gpu " \
  > setup_scripts/setup_conda_nv_bins.sh
echo "PATH=${PWD}/conda_nv_bins:\$PATH" \
  >> setup_scripts/setup_conda_nv_bins.sh
```

We are now ready to install the NVIDIA HPC SDK proper. Make sure you do this on a GPU-enabled node, as the installer checks for an existing CUDA driver.

Fell free to download the latest version froom the [NVIDIA official site](https://developer.nvidia.com/hpc-sdk). 

The following instructions will download and unpack version 20.9:

```
wget https://developer.download.nvidia.com/hpc-sdk/20.9/nvhpc_2020_209_Linux_x86_64_cuda_11.0.tar.gz
tar xpzf nvhpc_2020_209_Linux_x86_64_cuda_11.0.tar.gz
rm -f nvhpc_2020_209_Linux_x86_64_cuda_11.0.tar.gz
```

Once you have the install directory unpacked, you need to patch the installer to use the right GCC version; then you are ready to run the actual installer:

```
source setup_scripts/setup_conda_nv_bins.sh

# must patch the  install scripts to find the right gcc
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/install 
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/*/*/compilers/bin/makelocalrc
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/*/*/compilers/bin/addlocalrc
sed -i -e "s#PATH=/#PATH=$PWD/conda_nv_bins:/#g" \
  nvhpc_*/install_components/install_cuda

(cd nvhpc_*; ./install)
# Select "Single system install"
# Expand $PWD/hpc_sdk as install dir

echo "PATH=\$PATH:`ls -d $PWD/hpc_sdk/*/202*/compilers/bin`" \
  > setup_scripts/setup_nv_hpc_bins.sh

# h5c++ patch
mkdir conda_h5
cp $CONDA_PREFIX/bin/h5c++ conda_h5/

# This works on linux with gcc ..
sed -i \
  "s#x86_64-conda.*-linux-gnu-c++#pgc++ -I`ls -d $PWD/hpc_sdk/*/202*/compilers/include`#g" \
  conda_h5/h5c++ 
sed -i \
  's#H5BLD_CXXFLAGS=".*"#H5BLD_CXXFLAGS=" -fvisibility-inlines-hidden -std=c++17 -fPIC -O2 -I${includedir}"#g'  \
  conda_h5/h5c++
sed -i \
  's#H5BLD_CPPFLAGS=".*"#H5BLD_CPPFLAGS=" -I${includedir} -DNDEBUG -D_FORTIFY_SOURCE=2 -O2"#g' \
  conda_h5/h5c++
sed -i \
  's#H5BLD_LDFLAGS=".*"#H5BLD_LDFLAGS=" -L${prefix}/x86_64-conda-linux-gnu/sysroot/usr/lib64/ -L${libdir} -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,-rpath,\\\\\\$ORIGIN/../x86_64-conda-linux-gnu/sysroot/usr/lib64/ -Wl,-rpath,\\\\\\$ORIGIN/../lib -Wl,-rpath,${prefix}/x86_64-conda-linux-gnu/sysroot/usr/lib64/ -Wl,-rpath,${libdir}"#g' \
 conda_h5/h5c++

cat > setup_nv_h5.sh  << EOF
source $PWD/setup_scripts/setup_conda_nv_bins.sh
source $PWD/setup_scripts/setup_nv_hpc_bins.sh

PATH=${PWD}/conda_h5:\$PATH

# pgc++ does not define it, but gcc libraries expect it
export CPPFLAGS=-D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=0

EOF
```

For convenience, we also create a setup script that will be used to properly setup the environment anytime needed.


## Compiling the GPU-enabled UniFrac with the NVIDIA HPC SDK

In order to compile UniFrac, you will need both the Anaconda dependencies, the NVIDIA HPC SDK and the UniFrac source code.
The first two you setup above, and can enable with the helper script; the later can be imported (once) with git.

```
source setup_nv_h5.sh

git clone https://github.com/biocore/unifrac.git

# save CPU version of binaries
mv $CONDA_PREFIX/bin/ssu $CONDA_PREFIX/bin/ssu.cpu
mv $CONDA_PREFIX/bin/faithpd $CONDA_PREFIX/bin/faithpd.cpu
mv $CONDA_PREFIX/lib/libssu.so $CONDA_PREFIX/bin/libssu.so.cpu

(cd unifrac/sucpp/ && make && make main && make api)
```

And you are all done.
The UniFrac binary and libraries in the Anaconda environment are now the GPU-enabled ones.

Note: If you do not want the cutting edge UniFrac from git, you can get a tarball of the released versions. The first fully GPU-enabled version was [0.20.1](https://codeload.github.com/biocore/unifrac/tar.gz/0.20.1).
 
