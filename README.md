# UniFrac

[![Build Status](https://travis-ci.org/biocore/unifrac.svg?branch=master)](https://travis-ci.org/biocore/unifrac)

The *de facto* repository for UniFrac, based on an implementation of the Strided State UniFrac algorithm (manuscript in prep) which is faster, and uses less memory than [Fast UniFrac](http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html). Strided State UniFrac supports [Unweighted UniFrac](http://aem.asm.org/content/71/12/8228.abstract), [Weighted UniFrac](http://aem.asm.org/content/73/5/1576), [Generalized UniFrac](https://academic.oup.com/bioinformatics/article/28/16/2106/324465/Associating-microbiome-composition-with), [Variance Adjusted UniFrac](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-118) and [meta UniFrac](http://www.pnas.org/content/105/39/15076.short).

This repository produces a C API exposed via a shared library which can be linked against by any programming language. 

# Install

At this time, there are two primary ways to install the library. The first is through QIIME2, and the second is via `pip`. It is also possible to clone the repository and install using either the `sucpp/Makefile` or `setup.py` but those are not covered at this time. 

## Install (QIIME2)

The easiest way to use this library is through [QIIME2](https://docs.qiime2.org/2018.2/install/). The implementation of this algorithm is installed by default and is available under `qiime diversity beta-phylogenetic-alt`.

## Install (native)

To install, first the binary needs to be compiled. This assumes that the HDF5 
toolchain and libraries are available. More information about how to setup the
stack can be found [here](https://support.hdfgroup.org/HDF5/Tutor/compile.html). 

Assuming `h5c++` is in your path, the following should work:

    pip install -e . 

**Note**: if you are using `conda` we recommend installing HDF5 using the
`conda-forge` channel, for example:

    conda install -c conda-forge hdf5
        
# Examples of use

Below are a few light examples of different ways to use this library.

## QIIME2 

The to use this library through QIIME2, you need to provide a `FeatureTable[Frequency]` and a `Phylogeny[Rooted]` artifacts. An example of use is:

    qiime diversity beta-phylogenetic-alt --i-table table-evenly-samples.qza \
                                          --i-phylogeny a-tree.qza \
                                          --o-distance-matrix resulting-distance-matrix.qza \
                                          --p-metric unweighted_unifrac
                                          
## Python

The library can be accessed directly from within Python. If operating in this mode, the API methods are expecting a filepath to a BIOM-Format V2.1.0 table, and a filepath to a Newick formatted phylogeny.

    $ python
    Python 3.5.4 | packaged by conda-forge | (default, Aug 10 2017, 01:41:15)
    [GCC 4.2.1 Compatible Apple LLVM 6.1.0 (clang-602.0.53)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import unifrac
    >>> dir(unifrac)
    ['__all__', '__builtins__', '__cached__', '__doc__', '__file__', '__loader__', '__name__', '__package__', '__path__', '__spec__', '__version__', '_api', '_meta', '_methods', 'generalized', 'meta', 'pkg_resources', 'ssu', 'unweighted', 'weighted_normalized', 'weighted_unnormalized']
    >>> print(unifrac.unweighted.__doc__)
    Compute Unweighted UniFrac

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        The number of threads to use. Default of 1.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.

    Returns
    -------
    skbio.DistanceMatrix
        The resulting distance matrix.

    Raises
    ------
    IOError
        If the tree file is not found
        If the table is not found
    ValueError
        If the table does not appear to be BIOM-Format v2.1.
        If the phylogeny does not appear to be in Newick format.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).

## Command line

The methods can also be used directly through the command line after install:

    $ which ssu
    /Users/mcdonadt/miniconda3/envs/qiime2-2018.2/bin/ssu
    (qiime2-2018.2) 09:02:29 (mcdonadt@codinator):~$ ssu --help
    usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-n threads] [-a alpha] [--vaw]

        -i		The input BIOM table.
        -t		The input phylogeny in newick.
        -m		The method, [unweighted | weighted_normalized | weighted_unnormalized | generalized].
        -o		The output distance matrix.
        -n		[OPTIONAL] The number of threads, default is 1.
        -a		[OPTIONAL] Generalized UniFrac alpha, default is 1.
        -f		[OPTIONAL] Bypass tips, reduces compute by about 50%.
        --vaw	[OPTIONAL] Variance adjusted, default is to not adjust for variance.

    Citations:
        For UniFrac, please see:
            Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005
            Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06
            Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97
            Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133
        For Generalized UniFrac, please see:
            Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342
        For Variance Adjusted UniFrac, please see:
            Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118
            
## Shared library access

In addition to the above methods to access UniFrac, it is also possible to link against the shared library. The C API is described in `sucpp/api.hpp`, and examples of linking against this API can be found in `examples/`. 
