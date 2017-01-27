# q2-state-unifrac

An implementation of the strided state UniFrac algorithm (manuscript in prep).

# Install

To install, first the binary needs to be compiled. This assumes that the HDF5 
toolchain and libraries are available. More information about how to setup the
stack can be found [here](https://support.hdfgroup.org/HDF5/Tutor/compile.html). 

Assuming `h5c++` is in your path, the following should work:

    cd sucpp
    make main  # this will copy the su binary to the parent directory
    pip install -e . 
