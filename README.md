# q2-state-unifrac

[![Build Status](https://travis-ci.org/wasade/q2-state-unifrac.svg?branch=master)](https://travis-ci.org/wasade/q2-state-unifrac)

An implementation of the strided state UniFrac algorithm (manuscript in prep).

# Install

To install, first the binary needs to be compiled. This assumes that the HDF5 
toolchain and libraries are available. More information about how to setup the
stack can be found [here](https://support.hdfgroup.org/HDF5/Tutor/compile.html). 

Assuming `h5c++` is in your path, the following should work:

    pip install -e . 

**Note**: if you are using `conda` we recommend installing HDF5 using the
`conda-forge` channel, for example:

    conda install -c conda-forge hdf5

On OSX, it may be necessary to update the development target if c++11x is not being found.

    export MACOSX_DEPLOYMENT_TARGET=10.12
