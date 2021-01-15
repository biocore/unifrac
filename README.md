# UniFrac
##### Canonically pronounced *yew-nih-frak*

[![Build Status](https://travis-ci.com/biocore/unifrac.svg?branch=master)](https://travis-ci.com/biocore/unifrac)

The *de facto* repository for high-performance phylogenetic diversity calculations. The methods in this repository are based on an implementation of the [Strided State UniFrac](https://www.nature.com/articles/s41592-018-0187-8) algorithm which is faster, and uses less memory than [Fast UniFrac](http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html). Strided State UniFrac supports [Unweighted UniFrac](http://aem.asm.org/content/71/12/8228.abstract), [Weighted UniFrac](http://aem.asm.org/content/73/5/1576), [Generalized UniFrac](https://academic.oup.com/bioinformatics/article/28/16/2106/324465/Associating-microbiome-composition-with), [Variance Adjusted UniFrac](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-118) and [meta UniFrac](http://www.pnas.org/content/105/39/15076.short), in both double and single precision (fp32).
This repository also includes Stacked Faith (manuscript in preparation), a method for calculating Faith's PD that is faster and uses less memory than the Fast UniFrac-based [reference implementation](http://scikit-bio.org/).

This repository produces a C API exposed via a shared library which can be linked against by any programming language. 

# Citation

A detailed description of the Strided State UniFrac algorithm can be found in [McDonald et al. 2018 Nature Methods](https://www.nature.com/articles/s41592-018-0187-8). Please note that this package implements multiple UniFrac variants, which may have their own citation. Details can be found in the help output from the command line interface in the citations section, and is included immediately below:

    ssu
    For UniFrac, please see:
        McDonald et al. Nature Methods 2018; DOI: 10.1038/s41592-018-0187-8
        Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005
        Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06
        Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97
        Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133
    For Generalized UniFrac, please see: 
        Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342
    For Variance Adjusted UniFrac, please see: 
        Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118

    faithpd
    For Faith's PD, please see:
        Faith Biological Conservation 1992; DOI: 10.1016/0006-3207(92)91201-3

# Install

At this time, there are three primary ways to install the library. The first is through QIIME2, the second is through `bioconda`, and the third is via `pip`. It is also possible to clone the repository and install the C++ API with `sucpp/Makefile` or python bindings with `setup.py`. 

Compilation has been performed on both LLVM 9.0.0 (OS X >= 10.12) or GCC 4.9.2 (Centos >= 6) and HDF5 >= 1.8.17. Python installation requires Python >= 3.5, NumPy >= 1.12.1, scikit-bio >= 0.5.1, and Cython >= 0.28.3. 

Installation time should be a few minutes at most.

## Install (QIIME2)

The easiest way to use this library is through [QIIME2](https://docs.qiime2.org/2019.7/install/). This library is installed by default with the QIIME 2 Core Distribution. Currently, this module is used for phylogenetic diversity calculations in `qiime diversity beta-phylogenetic` for UniFrac and `qiime diversity alpha-phylogenetic-alt` for Faith's PD.

## Install (bioconda)

This library can also be installed via a combination of `conda-forge` and `bioconda`:

```
conda install -c conda-forge -c bioconda unifrac
```

Note: Only the CPU version of the binaries is currently available in conda. 
The GPU version must either be [locally compiled using freely-available NVIDIA HPC SDK](docs/compile_gpu.README.md) 
or obtained [from a github branch](https://github.com/sfiligoi/unifrac/blob/v0.20.1-docs/docs/install_gpu.README.txt).

Note: If you desire a fully optimized the binaries for your CPU, you can [compile them locally](docs/compile_cpu.README.md).

## Install (pip)

```
pip install unifrac
```

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

To use Strided State UniFrac through QIIME2, you need to provide a `FeatureTable[Frequency]` and a `Phylogeny[Rooted]` artifacts. An example of use is:

    qiime diversity beta-phylogenetic --i-table table-evenly-sampled.qza \
                                      --i-phylogeny a-tree.qza \
                                      --o-distance-matrix resulting-distance-matrix.qza \
                                      --p-metric unweighted_unifrac

To use Stacked Faith through QIIME2, given similar artifacts, you can use:

    qiime diversity alpha-phylogenetic-alt --i-table table-evenly-sampled.qza \
                                           --i-phylogeny a-tree.qza \
                                           --o-alpha-diversity resulting-diversity-series.qza \
                                           --p-metric faith_Pd
                                          
## Python

The library can be accessed directly from within Python. If operating in this mode, the API methods are expecting a filepath to a BIOM-Format V2.1.0 table, and a filepath to a Newick formatted phylogeny.

    $ python
    Python 3.7.8 | packaged by conda-forge | (default, Nov 27 2020, 19:24:58) 
    [GCC 9.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import unifrac
    >>> dir(unifrac)
    ['__all__', '__builtins__', '__cached__', '__doc__', '__file__', '__loader__', '__name__', 
     '__package__', '__path__', '__spec__', '__version__', '_api', '_meta', '_methods', 
     'faith_pd', 
     'generalized', 'generalized_fp32', 'generalized_fp32_to_file', 'generalized_to_file', 
     'h5pcoa', 'h5unifrac', 'meta', 'pkg_resources', 'ssu', 'ssu_to_file', 
     'unweighted', 'unweighted_fp32', 'unweighted_fp32_to_file', 'unweighted_to_file', 
     'weighted_normalized', 'weighted_normalized_fp32', 'weighted_normalized_fp32_to_file', 'weighted_normalized_to_file', 
     'weighted_unnormalized', 'weighted_unnormalized_fp32', 'weighted_unnormalized_fp32_to_file', 'weighted_unnormalized_to_file']
    >>> print(unifrac.unweighted_fp32.__doc__)
    Compute Unweighted UniFrac using fp32 math

        Parameters
        ----------
        table : str
            A filepath to a BIOM-Format 2.1 file.
        phylogeny : str
            A filepath to a Newick formatted tree.
        threads : int, optional
            The number of threads to split it into. Default of 1.
        variance_adjusted : bool, optional
            Adjust for varianace or not. Default is False.
        bypass_tips : bool, optional
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

    >>> print(unifrac.unweighted_fp32_to_file.__doc__)
    Compute Unweighted UniFrac using fp32 math and write to file

        Parameters
        ----------
        table : str
            A filepath to a BIOM-Format 2.1 file.
        phylogeny : str
            A filepath to a Newick formatted tree.
        out_filename : str
            A filepath to the output file.
        pcoa_dims : int, optional
            Number of dimensions to use for PCoA compute.
            if set to 0, no PCoA is computed.
            Defaults of 10.
        threads : int, optional
            The number of threads to split it into. Default of 1.
        variance_adjusted : bool, optional
            Adjust for varianace or not. Default is False.
        bypass_tips : bool, optional
            Bypass the tips of the tree in the computation. This reduces compute
            by about 50%, but is an approximation.
        format : str, optional
            Output format to use. Defaults to "hdf5".
        buf_dirname : str, optional
            If set, the directory where the disk buffer is hosted,
            can be used to reduce the amount of memory needed.

        Returns
        -------
        str
            A filepath to the output file.

        Raises
        ------
        IOError
            If the tree file is not found
            If the table is not found
            If the output file cannot be created
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

    >>> print(unifrac.h5unifrac.__doc__)
    Read UniFrac from a hdf5 file

        Parameters
        ----------
        h5file : str
            A filepath to a hdf5 file.

        Returns
        -------
        skbio.DistanceMatrix
            The distance matrix.

        Raises
        ------
        OSError
            If the hdf5 file is not found
        KeyError
            If the hdf5 does not have the necessary fields

        References
        ----------
        .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
           comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
           (2005).
        .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
           powerful beta diversity measure for comparing communities based on
           phylogeny. BMC Bioinformatics 12:118 (2011).

	>>> print(unifrac.faith_pd.__doc__)
	Execute a call to the Stacked Faith API in the UniFrac package

		Parameters
		----------
		biom_filename : str
			A filepath to a BIOM 2.1 formatted table (HDF5)
		tree_filename : str
			A filepath to a Newick formatted tree

		Returns
		-------
		pd.Series
			Series of Faith's PD for each sample in `biom_filename`

		Raises
		------
		IOError
			If the tree file is not found
			If the table is not found
			If the table is empty
	

## Command line

The methods can also be used directly through the command line after install:

    $ which ssu
    /Users/<username>/miniconda3/envs/qiime2-20xx.x/bin/ssu
    $ ssu --help
    usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-n threads] [-a alpha] [-f]  [--vaw]
        [--mode [MODE]] [--start starting-stripe] [--stop stopping-stripe] [--partial-pattern <glob>]
        [--n-partials number_of_partitions] [--report-bare] [--format|-r out-mode]

        -i		The input BIOM table.
        -t		The input phylogeny in newick.
        -m		The method, [unweighted | weighted_normalized | weighted_unnormalized | generalized | 
                                 unweighted_fp32 | weighted_normalized_fp32 | weighted_unnormalized_fp32 | generalized_fp32].
        -o		The output distance matrix.
        -n		[OPTIONAL] The number of threads, default is 1.
        -a		[OPTIONAL] Generalized UniFrac alpha, default is 1.
        -f		[OPTIONAL] Bypass tips, reduces compute by about 50%.
        --vaw	[OPTIONAL] Variance adjusted, default is to not adjust for variance.
        --mode	[OPTIONAL] Mode of operation:
                                one-off : [DEFAULT] compute UniFrac.
                                partial : Compute UniFrac over a subset of stripes.
                                partial-report : Start and stop suggestions for partial compute.
                                merge-partial : Merge partial UniFrac results.
        --start	[OPTIONAL] If mode==partial, the starting stripe.
        --stop	[OPTIONAL] If mode==partial, the stopping stripe.
        --partial-pattern	[OPTIONAL] If mode==merge-partial, a glob pattern for partial outputs to merge.
        --n-partials	[OPTIONAL] If mode==partial-report, the number of partitions to compute.
        --report-bare	[OPTIONAL] If mode==partial-report, produce barebones output.
        --format|-r	[OPTIONAL]  Output format:
                                 ascii : [DEFAULT] Original ASCII format.
                                 hfd5 : HFD5 format.  May be fp32 or fp64, depending on method.
                                 hdf5_fp32 : HFD5 format, using fp32 precision.
                                 hdf5_fp64 : HFD5 format, using fp64 precision.
        --pcoa	[OPTIONAL] Number of PCoA dimensions to compute (default: 10, do not compute if 0)
        --diskbuf	[OPTIONAL] Use a disk buffer to reduce memory footprint. Provide path to a fast partition (ideally NVMe).

    Citations: 
        For UniFrac, please see:
            McDonald et al. Nature Methods 2018; DOI: 10.1038/s41592-018-0187-8
            Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005
            Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06
            Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97
            Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133
        For Generalized UniFrac, please see: 
            Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342
        For Variance Adjusted UniFrac, please see: 
            Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118

    $ which faithpd
    /Users/<username>/miniconda3/envs/qiime2-20xx.x/bin/faithpd
    $ faithpd --help
	usage: faithpd -i <biom> -t <newick> -o <out.txt>

		-i          The input BIOM table.
		-t          The input phylogeny in newick.
		-o          The output series.

	Citations: 
		For Faith's PD, please see:
			Faith Biological Conservation 1992; DOI: 10.1016/0006-3207(92)91201-3

            
## Shared library access

In addition to the above methods to access UniFrac, it is also possible to link against the shared library. The C API is described in `sucpp/api.hpp`, and examples of linking against this API can be found in `examples/`. 

## Minor test dataset

A small test `.biom` and `.tre` can be found in `sucpp/`. An example with expected output is below, and should execute in 10s of milliseconds:

    $ ssu -i sucpp/test.biom -t sucpp/test.tre -m unweighted -o test.out
    $ cat test.out
    	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
    Sample1	0	0.2	0.5714285714285714	0.6	0.5	0.2
    Sample2	0.2	0	0.4285714285714285	0.6666666666666666	0.6	0.3333333333333333
    Sample3	0.5714285714285714	0.4285714285714285	0	0.7142857142857143	0.8571428571428571	0.4285714285714285
    Sample4	0.6	0.6666666666666666	0.7142857142857143	0	0.3333333333333333	0.4
    Sample5	0.5	0.6	0.8571428571428571	0.3333333333333333	0	0.6
    Sample6	0.2	0.3333333333333333	0.4285714285714285	0.4	0.6	0
