# Compiling an optimized version of UniFrac

The binaries provided through conda should work on all major platforms.
However, you may be able to get slightly faster results by locally compiling your own binaries.

Note: UniFrac supports only Linux and MacOS builds.
One can however run on Windows systems, too, using [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install-win10).


## Anaconda 

UniFrac has several dependencies, which we assume come via [Anaconda](https://www.anaconda.com/products/individual).

The instructions below has been tested with version 2021.05.

In case you have never used Anaconda below, here are the installation instruction:

### On Linux:
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
chmod a+x Anaconda3-2021.05-Linux-x86_64.sh
./Anaconda3-2021.05-Linux-x86_64.sh
#log out and back in
```

### On MacOS:
```
curl -o Anaconda3-2021.05-MacOSX-x86_64.sh  https://repo.anaconda.com/archive/Anaconda3-2021.05-MacOSX-x86_64.sh
chmod a+x Anaconda3-2021.05-MacOSX-x86_64.sh
./Anaconda3-2021.05-MacOSX-x86_64.sh
#log out and back in
```


## Create a dedicated environment

While it is possible to build your own UniFrac binaries in any Anaconda environment, we assume a dedicated one in this document.
We call it **unifrac-cpu**.

Note: If you decide to change the used environment, you will have to make the appropriate changes to the commands below. 

To create our **unifrac-cpu** with all the needed dependencies, run:

```
# create and activate unifrac-gpu Anaconda environment
conda create --name unifrac-cpu -c conda-forge -c bioconda unifrac
conda activate unifrac-cpu
```

On linux, you will need the conda provided gcc compiler:
```
# Linux only
conda install -c conda-forge -c bioconda gxx_linux-64=9.3
```

On MacOs, you will need the conda provided clang compiler:
```
# MacOS only
conda install -c conda-forge -c bioconda clangxx_osx-64=10.0.0
# Also add library which is not imported automatically
conda install -c conda-forge -c bioconda liblapacke
```

Finally, on both Linux and MacOS, add two additinal compile-only dependencies:
```
conda install -c conda-forge -c bioconda hdf5-static mkl-include
```

*Note:* On MacOS, the conda-provided clang is not compatible with XCode 12. 
        If you have it installed, you must either remove it or [downgrade to XCode 11](https://developer.apple.com/download/more/?=command%20line%20tools). 

## Compile UniFrac

Assuming you are in the environment above, fetch the UniFrac source code and perform the build.
The binaries will be automatically put in the conda path.

```
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
The UniFrac binary and libraries in the Anaconda environment are now the locally compiled ones.

Note: If you do not want the cutting edge UniFrac from git, you can get a tarball of the released versions. These instructions have been tested with version [0.20.2](https://codeload.github.com/biocore/unifrac/tar.gz/0.20.2).

