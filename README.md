# UniFrac

[![Build Status](https://travis-ci.org/biocore/unifrac.svg?branch=master)](https://travis-ci.org/biocore/unifrac)

The *de facto* repository for UniFrac, based on an implementation of the Strided State UniFrac algorithm (manuscript in prep) which is faster, and uses less memory than [Fast UniFrac](http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html). Strided State UniFrac supports [Unweighted UniFrac](http://aem.asm.org/content/71/12/8228.abstract), [Weighted UniFrac](http://aem.asm.org/content/73/5/1576), [Generalized UniFrac](https://academic.oup.com/bioinformatics/article/28/16/2106/324465/Associating-microbiome-composition-with), [Variance Adjusted UniFrac](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-118) and [meta UniFrac](http://www.pnas.org/content/105/39/15076.short).

This repository produces a C API exposed via a shared library which can be linked against by any programming language. 

# Install

To install, first the binary needs to be compiled. This assumes that the HDF5 
toolchain and libraries are available. More information about how to setup the
stack can be found [here](https://support.hdfgroup.org/HDF5/Tutor/compile.html). 

Assuming `h5c++` is in your path, the following should work:

    pip install -e . 

**Note**: if you are using `conda` we recommend installing HDF5 using the
`conda-forge` channel, for example:

    conda install -c conda-forge hdf5
