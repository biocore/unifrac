# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, UniFrac development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from warnings import warn
from functools import reduce
from operator import or_
from typing import Union
import collections.abc
import tempfile
import sys

import numpy as np
import pandas as pd
import skbio
import h5py
from bp import BP, write_newick
from skbio import TreeNode
from skbio.stats.distance._base import _build_results as _build_stat
from biom import Table
from biom.util import biom_open

import unifrac as qsu
from unifrac._meta import CONSOLIDATIONS


def is_biom_v210(f, ids=None):
    if not h5py.is_hdf5(f):
        return False
    with h5py.File(f, 'r') as fp:
        if 'format-version' not in fp.attrs:
            return False

        version = fp.attrs.get('format-version', None)

        if version is None:
            return False

        if tuple(version) != (2, 1):
            return False

        if ids is not None:
            for idel in fp['sample/ids']:
                if isinstance(idel, bytes):
                    ids.append(idel.decode('ascii'))
                else:
                    ids.append(idel)

    return True


def has_samples_biom_v210(f):
    # assumes this is already checked to be biom v210
    with h5py.File(f, 'r') as fp:
        if len(fp['sample/ids']) == 0:
            return False

    return True


def is_newick(f):
    sniffer = skbio.io.format.newick.newick.sniffer_function
    return sniffer(f)[0]


def has_atleast_two_samples(f):
    with h5py.File(f, 'r') as fp:
        assert 'sample/ids' in fp
        return len(fp['sample/ids']) >= 2


def has_negative(f):
    with h5py.File(f, 'r') as fp:
        assert 'sample/matrix/data' in fp
        return (fp['sample/matrix/data'][:] < 0).any()


def _validate(table, phylogeny, ids=None):
    if not is_biom_v210(table, ids):
        raise ValueError("Table does not appear to be a BIOM-Format v2.1")
    if not has_samples_biom_v210(table):
        raise ValueError("Table does not contain any samples")
    if not has_atleast_two_samples(table):
        raise ValueError("Table must have at least two samples")
    if has_negative(table):
        raise ValueError("Table has negatives")
    if not is_newick(phylogeny):
        raise ValueError("The phylogeny does not appear to be newick")


def _call_ssu(table, phylogeny, *args):
    if isinstance(table, Table) and isinstance(phylogeny, (TreeNode, BP)):
        if table.is_empty():
            raise ValueError("Table does not contain any samples")
        return qsu.ssu_inmem(table, phylogeny, *args)
    elif isinstance(table, str) and isinstance(phylogeny, str):
        ids = []
        _validate(table, phylogeny, ids)
        return qsu.ssu_fast(table, phylogeny, ids, *args)
    else:
        table_type = type(table)
        tree_type = type(phylogeny)
        raise ValueError(f"table ('{table_type}') and tree ('{tree_type}') "
                         f"are incompatible with the library call")


def _call_ssu_dense_pair(ids, sample1, sample2,
                         phylogeny, *args):
    if isinstance(phylogeny, (TreeNode, BP)):
        tree = phylogeny
    elif isinstance(phylogeny, str):
        tree = TreeNode.read(phylogeny)
    else:
        tree_type = type(phylogeny)
        raise ValueError(f"tree ('{tree_type}') "
                         f"is incompatible with the library call")
    return qsu.ssu_dense_pair(ids, sample1, sample2, tree, *args)


class _using_temp:
    def __init__(self, table, phylogeny):
        self._inmem_table = table
        self._inmem_phy = phylogeny
        self._tmpd = None
        self.table = None
        self.phylogeny = None

    def __enter__(self):
        if self._tmpd is not None:
            raise IOError("already have tmp reference")

        kwargs = {}
        if sys.version_info.minor >= 10:
            kwargs.update({"ignore_cleanup_errors": True})

        self._tmpd = tempfile.TemporaryDirectory(**kwargs)
        self.table = f"{self._tmpd.name}/table.biom"
        self.phylogeny = f"{self._tmpd.name}/phylogeny.nwk"

        with biom_open(self.table, 'w') as fp:
            self._inmem_table.to_hdf5(fp, 'asd')

        with open(self.phylogeny, 'w') as fp:
            if isinstance(self._inmem_phy, TreeNode):
                self._inmem_phy.write(fp)
            elif isinstance(self._inmem_phy, BP):
                write_newick(self._inmem_phy, fp, False)
            else:
                raise IOError("unknown phylogey type")

        return self

    def __exit__(self, *args, **kwargs):
        self._cleanup()

    def _cleanup(self):
        if self._tmpd is not None:
            self._tmpd.cleanup()
            self._tmpd = None
            self.table = None
            self.phylogeny = None

    def __del__(self):
        self._cleanup()


def _call_faith_pd(table, phylogeny):
    if isinstance(table, Table) and isinstance(phylogeny, (TreeNode, BP)):
        if table.is_empty():
            raise ValueError("Table does not contain any samples")
        with _using_temp(table, phylogeny) as tmp:
            return qsu._faith_pd(tmp.table, tmp.phylogeny)

    elif isinstance(table, str) and isinstance(phylogeny, str):
        ids = []
        _validate(table, phylogeny, ids)
        return qsu._faith_pd(table, phylogeny)

    else:
        table_type = type(table)
        tree_type = type(phylogeny)
        raise ValueError(f"table ('{table_type}') and tree ('{tree_type}') "
                         f"are incompatible with the library call")


def _call_ssu_to_file(table, phylogeny, *args):
    if isinstance(table, Table) and isinstance(phylogeny, (TreeNode, BP)):
        raise NotImplementedError("Direct to file support from in memory "
                                  "objects has not been implemented yet")
    elif isinstance(table, str) and isinstance(phylogeny, str):
        _validate(table, phylogeny)
        return qsu.ssu_to_file_v2(table, phylogeny, *args)
    else:
        table_type = type(table)
        tree_type = type(phylogeny)
        raise ValueError(f"table ('{table_type}') and tree ('{tree_type}') "
                         f"are incompatible with the library call")


#
# Functions that compute Unifrac and return a memory object
#
def unweighted(table: Union[str, Table],
               phylogeny: Union[str, TreeNode, BP],
               threads: int = 1,
               variance_adjusted: bool = False,
               bypass_tips: bool = False,
               n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute unweighted UniFrac

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'unweighted', variance_adjusted, 1.0,
                     bypass_tips, n_substeps)


def unweighted_fp64(table: Union[str, Table],
                    phylogeny: Union[str, TreeNode, BP],
                    threads: int = 1,
                    variance_adjusted: bool = False,
                    bypass_tips: bool = False,
                    n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute unweighted UniFrac using fp64 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'unweighted_fp64', variance_adjusted,
                     1.0, bypass_tips, n_substeps)


def unweighted_fp32(table: Union[str, Table],
                    phylogeny: Union[str, TreeNode, BP],
                    threads: int = 1,
                    variance_adjusted: bool = False,
                    bypass_tips: bool = False,
                    n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute unweighted UniFrac using fp32 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'unweighted_fp32', variance_adjusted,
                     1.0, bypass_tips, n_substeps)


def unweighted_unnormalized(table: Union[str, Table],
                            phylogeny: Union[str, TreeNode, BP],
                            threads: int = 1,
                            variance_adjusted: bool = False,
                            bypass_tips: bool = False,
                            n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute unweighted unnormalized UniFrac

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'unweighted_unnormalized',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def unweighted_unnormalized_fp64(table: Union[str, Table],
                                 phylogeny: Union[str, TreeNode, BP],
                                 threads: int = 1,
                                 variance_adjusted: bool = False,
                                 bypass_tips: bool = False,
                                 n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute unweighted unnormalized UniFrac using fp64 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'unweighted_unnormalized_fp64',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def unweighted_unnormalized_fp32(table: Union[str, Table],
                                 phylogeny: Union[str, TreeNode, BP],
                                 threads: int = 1,
                                 variance_adjusted: bool = False,
                                 bypass_tips: bool = False,
                                 n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute unweighted unnormalized UniFrac using fp32 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'unweighted_unnormalized_fp32',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def weighted_normalized(table: Union[str, Table],
                        phylogeny: Union[str, TreeNode, BP],
                        threads: int = 1,
                        variance_adjusted: bool = False,
                        bypass_tips: bool = False,
                        n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute weighted normalized UniFrac

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'weighted_normalized',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def weighted_normalized_fp64(table: Union[str, Table],
                             phylogeny: Union[str, TreeNode, BP],
                             threads: int = 1,
                             variance_adjusted: bool = False,
                             bypass_tips: bool = False,
                             n_substeps: int = 1
                             ) -> skbio.DistanceMatrix:
    """Compute weighted normalized UniFrac using fp64 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'weighted_normalized_fp64',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def weighted_normalized_fp32(table: Union[str, Table],
                             phylogeny: Union[str, TreeNode, BP],
                             threads: int = 1,
                             variance_adjusted: bool = False,
                             bypass_tips: bool = False,
                             n_substeps: int = 1
                             ) -> skbio.DistanceMatrix:
    """Compute weighted normalized UniFrac using fp32 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'weighted_normalized_fp32',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def weighted_unnormalized(table: Union[str, Table],
                          phylogeny: Union[str, TreeNode, BP],
                          threads: int = 1,
                          variance_adjusted: bool = False,
                          bypass_tips: bool = False,
                          n_substeps: int = 1) -> skbio.DistanceMatrix:
    # noqa
    """Compute weighted unnormalized UniFrac

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op..
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'weighted_unnormalized',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def weighted_unnormalized_fp64(table: Union[str, Table],
                               phylogeny: Union[str, TreeNode, BP],
                               threads: int = 1,
                               variance_adjusted: bool = False,
                               bypass_tips: bool = False,
                               n_substeps: int = 1
                               ) -> skbio.DistanceMatrix:
    # noqa
    """Compute weighted unnormalized UniFrac using fp64 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        TDeprecated, no-op..
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'weighted_unnormalized_fp64',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def weighted_unnormalized_fp32(table: Union[str, Table],
                               phylogeny: Union[str, TreeNode, BP],
                               threads: int = 1,
                               variance_adjusted: bool = False,
                               bypass_tips: bool = False,
                               n_substeps: int = 1
                               ) -> skbio.DistanceMatrix:
    # noqa
    """Compute weighted unnormalized UniFrac using fp32 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        TDeprecated, no-op..
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu(table, phylogeny, 'weighted_unnormalized_fp32',
                     variance_adjusted, 1.0, bypass_tips, n_substeps)


def generalized(table: Union[str, Table],
                phylogeny: Union[str, TreeNode, BP],
                threads: int = 1,
                alpha: float = 1.0,
                variance_adjusted: bool = False,
                bypass_tips: bool = False,
                n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute Generalized UniFrac

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return weighted_normalized(table, phylogeny, threads,
                                   variance_adjusted, bypass_tips, n_substeps)
    else:
        return _call_ssu(table, phylogeny, 'generalized',
                         variance_adjusted, alpha, bypass_tips, n_substeps)


def generalized_fp64(table: Union[str, Table],
                     phylogeny: Union[str, TreeNode, BP],
                     threads: int = 1,
                     alpha: float = 1.0,
                     variance_adjusted: bool = False,
                     bypass_tips: bool = False,
                     n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute Generalized UniFrac using fp64 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return weighted_normalized_fp64(table, phylogeny, threads,
                                        variance_adjusted, bypass_tips,
                                        n_substeps)
    else:
        return _call_ssu(table, phylogeny, 'generalized_fp64',
                         variance_adjusted, alpha, bypass_tips, n_substeps)


def generalized_fp32(table: Union[str, Table],
                     phylogeny: Union[str, TreeNode, BP],
                     threads: int = 1,
                     alpha: float = 1.0,
                     variance_adjusted: bool = False,
                     bypass_tips: bool = False,
                     n_substeps: int = 1) -> skbio.DistanceMatrix:
    """Compute Generalized UniFrac using fp32 math

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.
    threads : int, optional
        Deprecated, no-op.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return weighted_normalized_fp32(table, phylogeny, threads,
                                        variance_adjusted, bypass_tips,
                                        n_substeps)
    else:
        return _call_ssu(table, phylogeny, 'generalized_fp32',
                         variance_adjusted, alpha, bypass_tips, n_substeps)


METHODS = {'unweighted': unweighted,
           'unweighted_unnormalized': unweighted_unnormalized,
           'weighted_normalized': weighted_normalized,
           'weighted_unnormalized': weighted_unnormalized,
           'generalized': generalized,
           'unweighted_fp64': unweighted_fp64,
           'unweighted_unnormalized_fp64': unweighted_unnormalized_fp64,
           'weighted_normalized_fp64': weighted_normalized_fp64,
           'weighted_unnormalized_fp64': weighted_unnormalized_fp64,
           'generalized_fp64': generalized_fp64,
           'unweighted_fp32': unweighted_fp32,
           'unweighted_unnormalized_fp32': unweighted_unnormalized_fp32,
           'weighted_normalized_fp32': weighted_normalized_fp32,
           'weighted_unnormalized_fp32': weighted_unnormalized_fp32,
           'generalized_fp32': generalized_fp32}


def meta(tables: tuple, phylogenies: tuple, weights: tuple = None,
         consolidation: str = None, method: str = None,
         threads: int = 1, variance_adjusted: bool = False,
         alpha: float = None, bypass_tips: bool = False,
         n_substeps: int = 1) -> \
         skbio.DistanceMatrix:
    """Compute meta UniFrac

    Parameters
    ----------
    tables : tuple of str
        Filepaths to BIOM-Format 2.1 files. This tuple is expected to be in
        index order with phylogenies.
    phylogenies : tuple of str
        Filepaths to Newick formatted trees. This tuple is expected to be in
        index order with tables.
    weights : tuple of float, optional
        The weight applied to each tree/table pair. This tuple is expected to
        be in index order with tables and phylogenies. Default is to weight
        each tree/table pair evenly.
    consolidation : str, optional
        The matrix consolidation method. The available choices are:
        'skipping_missing_matrices', 'missing_zero', 'missing_one',
        'skipping_missing_values'. The default is 'skipping_missing_values'.
    method : str
        The UniFrac method to use. The available choices are:
        'unweighted', 'unweighted_fp64', 'unweighted_fp32',
        'unweighted_unnormalized', 'unweighted_unnormalized_fp64',
        'unweighted_unnormalized_fp32',
        'weighted_unnormalized', 'weighted_unnormalized_fp64',
        'weighted_unnormalized_fp32',
        'weighted_normalized', 'weighted_normalized_fp64',
        'weighted_normalized_fp32',
        'generalized', 'generalized_fp64' and 'generalized_fp32'.
    threads : int, optional
        TDeprecated, no-op.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    UniFrac can be adapted to account for multiple genes, as originally
    done in [1]_.

    Generalized UniFrac was originally described in [2]_. Variance Adjusted
    UniFrac was originally described in [3]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    References
    ----------
    .. [1] Lozupone C. A., Hamady M., Cantarel B. L., Coutinho P. M.,
       Henrissat B., Gordon J. I. & Knight R. The convergence of carbohydrate
       active gene repertoires in human gut microbes. PNAS 105(39):15076-81
       (2008).
    .. [2] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [3] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if not len(tables):
        raise ValueError("No tables specified.")

    if not len(phylogenies):
        raise ValueError("No trees specified.")

    if len(tables) != len(phylogenies):
        raise ValueError("Number of trees and tables must be the same.")

    if weights is None:
        weights = tuple(1 for _ in phylogenies)
    else:
        if len(weights) != len(phylogenies):
            raise ValueError("Number of weights does not match number of "
                             "trees and tables.")

    if method is None:
        raise ValueError("No method specified.")
    method_ = METHODS.get(method.replace('-', '_'))
    if method_ is None:
        raise ValueError("Method (%s) unrecognized. Available methods are: %s"
                         % (method, ', '.join(METHODS.keys())))

    if consolidation is None:
        consolidation = 'skipping_missing_values'
    consolidation_ = CONSOLIDATIONS.get(consolidation.replace('-', '_'))
    if consolidation_ is None:
        raise ValueError("Consolidation (%s) unrecognized. Available "
                         "consolidations are: %s"
                         % (consolidation, ', '.join(CONSOLIDATIONS.keys())))

    if alpha is not None and method_ is not generalized:
        raise ValueError("The alpha parameter can only be set when the method "
                         "is set as 'generalized', the selected method is "
                         "'%s'." % method)

    kwargs = {'n_substeps': n_substeps,
              'bypass_tips': bypass_tips,
              'variance_adjusted': variance_adjusted}
    if alpha is not None:
        kwargs['alpha'] = alpha

    weights = np.array(weights, float)/sum(weights)
    dms = [method_(table, tree, **kwargs) for table, tree in zip(tables,
                                                                 phylogenies)]
    all_ids = sorted(reduce(or_, [set(dm.ids) for dm in dms]))
    dm = consolidation_(dms, [dm.ids for dm in dms], weights, all_ids)

    return skbio.DistanceMatrix(dm, ids=all_ids)


#
# Functions that compute Unifrac on a single dense pair
#
def unweighted_dense_pair(ids,
                          sample1,
                          sample2,
                          phylogeny: Union[str, TreeNode, BP],
                          variance_adjusted: bool = False,
                          bypass_tips: bool = False) -> float:
    """Compute unweighted UniFrac

    Parameters
    ----------
    ids : tuple or list of str
        Obeservation IDss
    sample1 : tuple or list of float
        First sample counts
    sample2 : tuple or list of float
        Second sample counts
    phylogeny : str
        A filepath to a Newick formatted tree.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.

    Returns
    -------
    float
        The resulting distance value.

    Raises
    ------
    IOError
        If the tree file is not found
    ValueError
        If the phylogeny does not appear to be in Newick format.

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    # Use the FP64 version, since the
    # other overheads hide any difference in performance
    return _call_ssu_dense_pair(ids, sample1, sample2,
                                phylogeny, 'unweighted_fp64',
                                variance_adjusted, 1.0,
                                bypass_tips)


def unweighted_unnormalized_dense_pair(ids,
                                       sample1,
                                       sample2,
                                       phylogeny: Union[str, TreeNode, BP],
                                       variance_adjusted: bool = False,
                                       bypass_tips: bool = False) -> float:
    """Compute unweighted unnormalized UniFrac

    Parameters
    ----------
    ids : tuple or list of str
        Obeservation IDss
    sample1 : tuple or list of float
        First sample counts
    sample2 : tuple or list of float
        Second sample counts
    phylogeny : str
        A filepath to a Newick formatted tree.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.

    Returns
    -------
    float
        The resulting distance value.

    Raises
    ------
    IOError
        If the tree file is not found
    ValueError
        If the phylogeny does not appear to be in Newick format.

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    # Use the FP64 version, since the
    # other overheads hide any difference in performance
    return _call_ssu_dense_pair(ids, sample1, sample2,
                                phylogeny, 'unweighted_unnormalized_fp64',
                                variance_adjusted, 1.0,
                                bypass_tips)


def weighted_normalized_dense_pair(ids,
                                   sample1,
                                   sample2,
                                   phylogeny: Union[str, TreeNode, BP],
                                   variance_adjusted: bool = False,
                                   bypass_tips: bool = False) -> float:
    """Compute weighted normalized UniFrac

    Parameters
    ----------
    ids : tuple or list of str
        Obeservation IDss
    sample1 : tuple or list of float
        First sample counts
    sample2 : tuple or list of float
        Second sample counts
    phylogeny : str
        A filepath to a Newick formatted tree.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.

    Returns
    -------
    float
        The resulting distance value.

    Raises
    ------
    IOError
        If the tree file is not found
    ValueError
        If the phylogeny does not appear to be in Newick format.

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    # Use the FP64 version, since the
    # other overheads hide any difference in performance
    return _call_ssu_dense_pair(ids, sample1, sample2,
                                phylogeny, 'weighted_normalized_fp64',
                                variance_adjusted, 1.0,
                                bypass_tips)


def weighted_unnormalized_dense_pair(ids,
                                     sample1,
                                     sample2,
                                     phylogeny: Union[str, TreeNode, BP],
                                     variance_adjusted: bool = False,
                                     bypass_tips: bool = False) -> float:
    """Compute weighted unnormalized UniFrac

    Parameters
    ----------
    ids : tuple or list of str
        Obeservation IDss
    sample1 : tuple or list of float
        First sample counts
    sample2 : tuple or list of float
        Second sample counts
    phylogeny : str
        A filepath to a Newick formatted tree.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.

    Returns
    -------
    float
        The resulting distance value.

    Raises
    ------
    IOError
        If the tree file is not found
    ValueError
        If the phylogeny does not appear to be in Newick format.

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    # Use the FP64 version, since the
    # other overheads hide any difference in performance
    return _call_ssu_dense_pair(ids, sample1, sample2,
                                phylogeny, 'weighted_unnormalized_fp64',
                                variance_adjusted, 1.0,
                                bypass_tips)


def generalized_dense_pair(ids,
                           sample1,
                           sample2,
                           phylogeny: Union[str, TreeNode, BP],
                           alpha: float = 1.0,
                           variance_adjusted: bool = False,
                           bypass_tips: bool = False) -> float:
    """Compute unweighted UniFrac

    Parameters
    ----------
    ids : tuple or list of str
        Obeservation IDss
    sample1 : tuple or list of float
        First sample counts
    sample2 : tuple or list of float
        Second sample counts
    phylogeny : str
        A filepath to a Newick formatted tree.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.

    Returns
    -------
    float
        The resulting distance value.

    Raises
    ------
    IOError
        If the tree file is not found
    ValueError
        If the phylogeny does not appear to be in Newick format.

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    # Use the FP64 version, since the
    # other overheads hide any difference in performance
    return _call_ssu_dense_pair(ids, sample1, sample2,
                                phylogeny, 'generalized_fp64',
                                variance_adjusted, alpha,
                                bypass_tips)


#
# Functions that compute Unifrac and write into a file
#

def unweighted_to_file(table: str,
                       phylogeny: str,
                       out_filename: str,
                       pcoa_dims: int = 10,
                       threads: int = 1,
                       variance_adjusted: bool = False,
                       bypass_tips: bool = False,
                       format: str = "",
                       buf_dirname: str = "",
                       n_substeps: int = 1,
                       n_subsamples: int = 1,
                       subsample_depth: int = 0,
                       subsample_with_replacement: bool = True,
                       permanova_perms: int = 0,
                       grouping_filename: str = "",
                       grouping_columns: str = "") -> str:
    """Compute unweighted UniFrac and write to file

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
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    format : str, optional
        Output format to use.
        Defaults to "hdf5" if n_subsamples<=1 else "hdf5_nodist"
    buf_dirname : str, optional
        If set, the directory where the disk buffer is hosted,
        can be used to reduce the amount of memory needed.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'unweighted',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def unweighted_fp64_to_file(table: str,
                            phylogeny: str,
                            out_filename: str,
                            pcoa_dims: int = 10,
                            threads: int = 1,
                            variance_adjusted: bool = False,
                            bypass_tips: bool = False,
                            format: str = "hdf5",
                            buf_dirname: str = "",
                            n_substeps: int = 1,
                            n_subsamples: int = 1,
                            subsample_depth: int = 0,
                            subsample_with_replacement: bool = True,
                            permanova_perms: int = 0,
                            grouping_filename: str = "",
                            grouping_columns: str = "") -> str:
    """Compute unweighted UniFrac using fp64 math and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'unweighted_fp64',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def unweighted_fp32_to_file(table: str,
                            phylogeny: str,
                            out_filename: str,
                            pcoa_dims: int = 10,
                            threads: int = 1,
                            variance_adjusted: bool = False,
                            bypass_tips: bool = False,
                            format: str = "hdf5",
                            buf_dirname: str = "",
                            n_substeps: int = 1,
                            n_subsamples: int = 1,
                            subsample_depth: int = 0,
                            subsample_with_replacement: bool = True,
                            permanova_perms: int = 0,
                            grouping_filename: str = "",
                            grouping_columns: str = "") -> str:
    """Compute unweighted UniFrac using fp32 math and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'unweighted_fp32',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def unweighted_unnormalized_to_file(table: str,
                                    phylogeny: str,
                                    out_filename: str,
                                    pcoa_dims: int = 10,
                                    threads: int = 1,
                                    variance_adjusted: bool = False,
                                    bypass_tips: bool = False,
                                    format: str = "",
                                    buf_dirname: str = "",
                                    n_substeps: int = 1,
                                    n_subsamples: int = 1,
                                    subsample_depth: int = 0,
                                    subsample_with_replacement: bool = True,
                                    permanova_perms: int = 0,
                                    grouping_filename: str = "",
                                    grouping_columns: str = "") -> str:
    """Compute unweighted unnormalized UniFrac and write to file

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
        Deprecated, no-op.
    variance_adjusted : bool, optional
        Adjust for varianace or not. Default is False.
    bypass_tips : bool, optional
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    format : str, optional
        Output format to use.
        Defaults to "hdf5" if n_subsamples<=1 else "hdf5_nodist"
    buf_dirname : str, optional
        If set, the directory where the disk buffer is hosted,
        can be used to reduce the amount of memory needed.
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'unweighted_unnormalized',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def unweighted_unnormalized_fp64_to_file(table: str,
                                         phylogeny: str,
                                         out_filename: str,
                                         pcoa_dims: int = 10,
                                         threads: int = 1,
                                         variance_adjusted: bool = False,
                                         bypass_tips: bool = False,
                                         format: str = "hdf5",
                                         buf_dirname: str = "",
                                         n_substeps: int = 1,
                                         n_subsamples: int = 1,
                                         subsample_depth: int = 0,
                                         subsample_with_replacement:
                                         bool = True,
                                         permanova_perms: int = 0,
                                         grouping_filename: str = "",
                                         grouping_columns: str = "") -> str:
    """Compute unweighted unnormalized UniFrac using fp64 math
    and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'unweighted_unnormalized_fp64',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def unweighted_unnormalized_fp32_to_file(table: str,
                                         phylogeny: str,
                                         out_filename: str,
                                         pcoa_dims: int = 10,
                                         threads: int = 1,
                                         variance_adjusted: bool = False,
                                         bypass_tips: bool = False,
                                         format: str = "hdf5",
                                         buf_dirname: str = "",
                                         n_substeps: int = 1,
                                         n_subsamples: int = 1,
                                         subsample_depth: int = 0,
                                         subsample_with_replacement:
                                         bool = True,
                                         permanova_perms: int = 0,
                                         grouping_filename: str = "",
                                         grouping_columns: str = "") -> str:
    """Compute unweighted unnormalized UniFrac using fp32 math
    and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, and while its application to
    Unweighted UniFrac was not described, factoring in the variance adjustment
    is still feasible and so it is exposed. Current implementation is
    described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'unweighted_unnormalized_fp32',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def weighted_normalized_to_file(table: str,
                                phylogeny: str,
                                out_filename: str,
                                pcoa_dims: int = 10,
                                threads: int = 1,
                                variance_adjusted: bool = False,
                                bypass_tips: bool = False,
                                format: str = "hdf5",
                                buf_dirname: str = "",
                                n_substeps: int = 1,
                                n_subsamples: int = 1,
                                subsample_depth: int = 0,
                                subsample_with_replacement: bool = True,
                                permanova_perms: int = 0,
                                grouping_filename: str = "",
                                grouping_columns: str = "") -> str:
    """Compute weighted normalized UniFrac and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'weighted_normalized',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def weighted_normalized_fp64_to_file(table: str,
                                     phylogeny: str,
                                     out_filename: str,
                                     pcoa_dims: int = 10,
                                     threads: int = 1,
                                     variance_adjusted: bool = False,
                                     bypass_tips: bool = False,
                                     format: str = "hdf5",
                                     buf_dirname: str = "",
                                     n_substeps: int = 1,
                                     n_subsamples: int = 1,
                                     subsample_depth: int = 0,
                                     subsample_with_replacement: bool = True,
                                     permanova_perms: int = 0,
                                     grouping_filename: str = "",
                                     grouping_columns: str = "") -> str:
    """Compute weighted normalized UniFrac using fp64 math and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'weighted_normalized_fp64',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def weighted_normalized_fp32_to_file(table: str,
                                     phylogeny: str,
                                     out_filename: str,
                                     pcoa_dims: int = 10,
                                     threads: int = 1,
                                     variance_adjusted: bool = False,
                                     bypass_tips: bool = False,
                                     format: str = "hdf5",
                                     buf_dirname: str = "",
                                     n_substeps: int = 1,
                                     n_subsamples: int = 1,
                                     subsample_depth: int = 0,
                                     subsample_with_replacement: bool = True,
                                     permanova_perms: int = 0,
                                     grouping_filename: str = "",
                                     grouping_columns: str = "") -> str:
    """Compute weighted normalized UniFrac using fp32 math and write to file

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
        Deprecated, no-op.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'weighted_normalized_fp32',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def weighted_unnormalized_to_file(table: str,
                                  phylogeny: str,
                                  out_filename: str,
                                  pcoa_dims: int = 10,
                                  threads: int = 1,
                                  variance_adjusted: bool = False,
                                  bypass_tips: bool = False,
                                  format: str = "hdf5",
                                  buf_dirname: str = "",
                                  n_substeps: int = 1,
                                  n_subsamples: int = 1,
                                  subsample_depth: int = 0,
                                  subsample_with_replacement: bool = True,
                                  permanova_perms: int = 0,
                                  grouping_filename: str = "",
                                  grouping_columns: str = "") -> str:
    """Compute weighted unnormalized UniFrac and write it to file

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
        TDeprecated, no-op..
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'weighted_unnormalized',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def weighted_unnormalized_fp64_to_file(table: str,
                                       phylogeny: str,
                                       out_filename: str,
                                       pcoa_dims: int = 10,
                                       threads: int = 1,
                                       variance_adjusted: bool = False,
                                       bypass_tips: bool = False,
                                       format: str = "hdf5",
                                       buf_dirname: str = "",
                                       n_substeps: int = 1,
                                       n_subsamples: int = 1,
                                       subsample_depth: int = 0,
                                       subsample_with_replacement: bool = True,
                                       permanova_perms: int = 0,
                                       grouping_filename: str = "",
                                       grouping_columns: str = "") -> str:
    """Compute weighted unnormalized UniFrac using fp64 math and write to file

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
        TDeprecated, no-op..
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'weighted_unnormalized_fp64',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def weighted_unnormalized_fp32_to_file(table: str,
                                       phylogeny: str,
                                       out_filename: str,
                                       pcoa_dims: int = 10,
                                       threads: int = 1,
                                       variance_adjusted: bool = False,
                                       bypass_tips: bool = False,
                                       format: str = "hdf5",
                                       buf_dirname: str = "",
                                       n_substeps: int = 1,
                                       n_subsamples: int = 1,
                                       subsample_depth: int = 0,
                                       subsample_with_replacement: bool = True,
                                       permanova_perms: int = 0,
                                       grouping_filename: str = "",
                                       grouping_columns: str = "") -> str:
    """Compute weighted unnormalized UniFrac using fp32 math and write to file

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
        TDeprecated, no-op..
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Weighted UniFrac was originally described in [1]_. Variance Adjusted
    Weighted UniFrac was originally described in [2]_. Current implementation
    is described in [3]_.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """
    return _call_ssu_to_file(table, phylogeny, out_filename,
                             'weighted_unnormalized_fp32',
                             variance_adjusted, 1.0, bypass_tips,
                             n_substeps, format,
                             n_subsamples,
                             subsample_depth, subsample_with_replacement,
                             pcoa_dims,
                             permanova_perms,
                             grouping_filename, grouping_columns,
                             buf_dirname)


def generalized_to_file(table: str,
                        phylogeny: str,
                        out_filename: str,
                        pcoa_dims: int = 10,
                        threads: int = 1,
                        alpha: float = 1.0,
                        variance_adjusted: bool = False,
                        bypass_tips: bool = False,
                        format: str = "hdf5",
                        buf_dirname: str = "",
                        n_substeps: int = 1,
                        n_subsamples: int = 1,
                        subsample_depth: int = 0,
                        subsample_with_replacement: bool = True,
                        permanova_perms: int = 0,
                        grouping_filename: str = "",
                        grouping_columns: str = "") -> str:
    """Compute Generalized UniFrac and write to file

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
        TDeprecated, no-op.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return _call_ssu_to_file(table, phylogeny, out_filename,
                                 'weighted_normalized',
                                 variance_adjusted, 1.0, bypass_tips,
                                 n_substeps, format,
                                 n_subsamples,
                                 subsample_depth, subsample_with_replacement,
                                 pcoa_dims,
                                 permanova_perms,
                                 grouping_filename, grouping_columns,
                                 buf_dirname)
    else:
        return _call_ssu_to_file(table, phylogeny, out_filename,
                                 'generalized',
                                 variance_adjusted, alpha, bypass_tips,
                                 n_substeps, format,
                                 n_subsamples,
                                 subsample_depth, subsample_with_replacement,
                                 pcoa_dims,
                                 permanova_perms,
                                 grouping_filename, grouping_columns,
                                 buf_dirname)


def generalized_fp64_to_file(table: str,
                             phylogeny: str,
                             out_filename: str,
                             pcoa_dims: int = 10,
                             threads: int = 1,
                             alpha: float = 1.0,
                             variance_adjusted: bool = False,
                             bypass_tips: bool = False,
                             format: str = "hdf5",
                             buf_dirname: str = "",
                             n_substeps: int = 1,
                             n_subsamples: int = 1,
                             subsample_depth: int = 0,
                             subsample_with_replacement: bool = True,
                             permanova_perms: int = 0,
                             grouping_filename: str = "",
                             grouping_columns: str = "") -> str:
    """Compute Generalized UniFrac using fp64 math and write to file

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
        TDeprecated, no-op.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return _call_ssu_to_file(table, phylogeny, out_filename,
                                 'weighted_normalized_fp64',
                                 variance_adjusted, 1.0, bypass_tips,
                                 n_substeps, format,
                                 n_subsamples,
                                 subsample_depth, subsample_with_replacement,
                                 pcoa_dims,
                                 permanova_perms,
                                 grouping_filename, grouping_columns,
                                 buf_dirname)
    else:
        return _call_ssu_to_file(table, phylogeny, out_filename,
                                 'generalized_fp64',
                                 variance_adjusted, alpha, bypass_tips,
                                 n_substeps, format,
                                 n_subsamples,
                                 subsample_depth, subsample_with_replacement,
                                 pcoa_dims,
                                 permanova_perms,
                                 grouping_filename, grouping_columns,
                                 buf_dirname)


def generalized_fp32_to_file(table: str,
                             phylogeny: str,
                             out_filename: str,
                             pcoa_dims: int = 10,
                             threads: int = 1,
                             alpha: float = 1.0,
                             variance_adjusted: bool = False,
                             bypass_tips: bool = False,
                             format: str = "hdf5",
                             buf_dirname: str = "",
                             n_substeps: int = 1,
                             n_subsamples: int = 1,
                             subsample_depth: int = 0,
                             subsample_with_replacement: bool = True,
                             permanova_perms: int = 0,
                             grouping_filename: str = "",
                             grouping_columns: str = "") -> str:
    """Compute Generalized UniFrac using fp32 math and write to file

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
        TDeprecated, no-op.
    alpha : float, optional
        The level of contribution of high abundance branches. Higher alpha
        increases the contribution of from high abundance branches while lower
        alpha reduces the contribution. Alpha was originally defined over the
        range [0, 1]. Default is 1.0.
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
    n_substeps : int, optional
        Internally split the problem in substeps for reduced memory footprint.
    n_subsamples : int
        If >1, perform multiple subsamples.
    subsample_depth : int
        Depth of subsampling, if >0
    subsample_with_replacement : bool
        Use subsampling with replacement? (only True supported in 1.3)
    permanova_perms : int
        If not 0, compute PERMANOVA using that many permutations
    grouping_filename : str
        The TSV filename containing grouping information
    grouping_columns : str
        The columns to use for grouping

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

    Environment variables
    ---------------------
    OMP_NUM_THREADS
        Number of CPU cores to use. If not defined, use all detected cores.
    UNIFRAC_USE_GPU
        Enable or disable GPU offload. If not defined, autodetect.
    ACC_DEVICE_NUM
        The GPU to use. If not defined, the first GPU will be used.

    Notes
    -----
    Generalized UniFrac was originally described in [1]_. Variance Adjusted
    UniFrac was originally described in [2]_, but was not described in as
    applied to Generalized UniFrac. It is feasible to do, so it is exposed
    here.

    An alpha of 1.0 is Weighted normalized UniFrac. An alpha of 0.0 is
    approximately Unweighted UniFrac, and is if the proportions are
    dichotomized.

    References
    ----------
    .. [1] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann C., Lewis, J.,
       Wu, G. D., Collman R. G., Bushman, F. D. & Hongzhe L. Associating
       microbiome composition with environmental covariates using generalized
       UniFrac distances. Bioinformatics 28(16), 2106–2113 (2012).
    .. [2] Chang, Q., Luan, Y. & Sun, F. Variance adjusted weighted UniFrac: a
       powerful beta diversity measure for comparing communities based on
       phylogeny. BMC Bioinformatics 12:118 (2011).
    """
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return _call_ssu_to_file(table, phylogeny, out_filename,
                                 'weighted_normalized_fp32',
                                 variance_adjusted, 1.0, bypass_tips,
                                 n_substeps, format,
                                 n_subsamples,
                                 subsample_depth, subsample_with_replacement,
                                 pcoa_dims,
                                 permanova_perms,
                                 grouping_filename, grouping_columns,
                                 buf_dirname)
    else:
        return _call_ssu_to_file(table, phylogeny, out_filename,
                                 'generalized_fp32',
                                 variance_adjusted, alpha, bypass_tips,
                                 n_substeps, format,
                                 n_subsamples,
                                 subsample_depth, subsample_with_replacement,
                                 pcoa_dims,
                                 permanova_perms,
                                 grouping_filename, grouping_columns,
                                 buf_dirname)

#
# Functions that read Unifrac from hdf5 files
#


def h5unifrac(h5file: str) -> skbio.DistanceMatrix:
    """Read UniFrac distance matrix from a hdf5 file

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
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """

    with h5py.File(h5file, "r") as f_u:
        if 'matrix:0' in f_u.keys():
            # multi format
            dm = skbio.DistanceMatrix(
               f_u['matrix:0'][:, :],
               [c.decode('ascii') for c in f_u['order'][:]])
        else:
            # single format
            dm = skbio.DistanceMatrix(
               f_u['matrix'][:, :],
               [c.decode('ascii') for c in f_u['order'][:]])

    return dm


class H5UnifracTuple(collections.abc.Sequence):
    """Read all UniFrac distance matrices from a hdf5 file"""

    def __init__(self, h5file: str):
        self.f_u = h5py.File(h5file, "r")
        self.order = [c.decode('ascii') for c in self.f_u['order'][:]]
        # cache some often used values
        self.nels = None
        self.cached_idx = None
        self.cached_el = None

    def __getitem__(self, i: int) -> skbio.DistanceMatrix:
        if i == self.cached_idx:
            return self.cached_el
        i_str = 'matrix:%i' % i
        if i == 0:
            if 'matrix' in self.f_u.keys():
                # single format
                i_str = 'matrix'
        el = skbio.DistanceMatrix(self.f_u[i_str][:, :],
                                  self.order)
        # if it did not throw, cache
        self.cached_idx = i
        self.cached_el = el
        return self.cached_el

    def __len__(self) -> int:
        if self.nels is None:
            i = 0
            if 'matrix' in self.f_u.keys():
                # single format
                i = 1
            else:
                # multi format
                while 'matrix:%i' % i in self.f_u.keys():
                    i = i + 1
            self.nels = i
        return self.nels

    def close(self):
        """Explicitly close the underlying file descriptor"""
        self.f_u.close()
        # invalidate all other cache values
        self.order = None
        self.nels = 0
        self.cached_idx = None
        self.cached_el = None


def h5unifrac_all(h5file: str) -> H5UnifracTuple:
    """Read all UniFrac distance matrices from a hdf5 file

    Parameters
    ----------
    h5file : str
        A filepath to a hdf5 file.

    Returns
    -------
    H5UnifracTuple
        A collection of distance matrices.

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
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """

    return H5UnifracTuple(h5file)


def _build_pcoa(f_u, long_method_name, order_index,
                eigval_key, samples_key, prop_key):
    axis_labels = ["PC%d" % i for i in
                   range(1, len(f_u[eigval_key][:]) + 1)]

    pc = skbio.OrdinationResults(
              short_method_name="PCoA",
              long_method_name=long_method_name,
              eigvals=pd.Series(f_u[eigval_key][:], index=axis_labels),
              samples=pd.DataFrame(f_u[samples_key][:, :],
                                   index=order_index,
                                   columns=axis_labels),
              proportion_explained=pd.Series(
                                     f_u[prop_key][:],
                                     index=axis_labels))

    return pc


def h5pcoa(h5file: str) -> skbio.OrdinationResults:
    """Read PCoA from a hdf5 file

    Parameters
    ----------
    h5file : str
        A filepath to a hdf5 file.

    Returns
    -------
    skbio.OrdinationResults
        The PCoA of the distance matrix

    Raises
    ------
    OSError
        If the hdf5 file is not found
    KeyError
        If the hdf5 does not have the necessary fields
    """

    with h5py.File(h5file, "r") as f_u:
        pcoa_method = f_u['pcoa_method'][0].decode('ascii')
        if 'FSVD' == pcoa_method:
            long_method_name = "Approximate Principal Coordinate Analysis" + \
                               " using FSVD"
        else:
            long_method_name = "Possibly Approximate Principal " + \
                               "Coordinate Analysis " + \
                               "using " + pcoa_method
        order_index = [c.decode('ascii')
                       for c in f_u['order'][:]]

        if 'pcoa_eigvals:0' in f_u.keys():
            # multi interface
            pc = _build_pcoa(f_u, long_method_name, order_index,
                             'pcoa_eigvals:0', 'pcoa_samples:0',
                             'pcoa_proportion_explained:0')
        else:
            # single interface
            pc = _build_pcoa(f_u, long_method_name, order_index,
                             'pcoa_eigvals', 'pcoa_samples',
                             'pcoa_proportion_explained')

    return pc


def h5pcoa_all(h5file: str) -> tuple:
    """Read all PCoAs from a hdf5 file

    Parameters
    ----------
    h5file : str
        A filepath to a hdf5 file.

    Returns
    -------
    tuple(skbio.OrdinationResults)
        The PCoAs of the distance matrix

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
    .. [3] Sfiligoi, I. et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
    """

    with h5py.File(h5file, "r") as f_u:
        pcoa_method = f_u['pcoa_method'][0].decode('ascii')
        if 'FSVD' == pcoa_method:
            long_method_name = "Approximate Principal Coordinate Analysis" + \
                               " using FSVD"
        else:
            long_method_name = "Possibly Approximate Principal " + \
                               "Coordinate Analysis " + \
                               "using " + pcoa_method
        order_index = [c.decode('ascii')
                       for c in f_u['order'][:]]

        if 'pcoa_eigvals' in f_u.keys():
            # single matrix single PCoA version
            pcs = [_build_pcoa(f_u, long_method_name, order_index,
                               'pcoa_eigvals', 'pcoa_samples',
                               'pcoa_proportion_explained')]
        else:
            # multi-matrix version
            pcs = []
            i = 0
            while 'pcoa_eigvals:%i' % i in f_u.keys():
                pcs.append(_build_pcoa(f_u, long_method_name, order_index,
                                       'pcoa_eigvals:%i' % i,
                                       'pcoa_samples:%i' % i,
                                       'pcoa_proportion_explained:%i' % i))
                i = i + 1

    return pcs


def h5permanova(h5file: str) -> pd.Series:
    """Read first PERMANOVA statistical test from a hdf5 file

    As describe in scikit-bio skbio.stats.distance.permanova.py,
    Permutational Multivariate Analysis of Variance (PERMANOVA) is a
    non-parametric method that tests whether two or more groups of objects
    are significantly different based on a categorical factor.

    Parameters
    ----------
    h5file : str
        A filepath to a hdf5 file.

    Returns
    -------
    pandas.Series
        Results of the statistical test, including ``test statistic`` and
        ``p-value``.

    Raises
    ------
    OSError
        If the hdf5 file is not found
    KeyError
        If the hdf5 does not have the necessary fields

    References
    ----------
    .. [1] Anderson, Marti J. "A new method for non-parametric multivariate
       analysis of variance." Austral Ecology 26.1 (2001): 32-46.
    """

    found = False
    with h5py.File(h5file, "r") as f_u:
        methods = f_u['stat_methods'][:]
        test_names = f_u['stat_test_names'][:]
        values = f_u['stat_values'][:]
        pvalues = f_u['stat_pvalues'][:]
        n_permutations = f_u['stat_n_permutations'][:]
        num_groups = f_u['stat_n_groups'][:]

        sample_size = len(f_u['order'][:])

        n_stats = len(methods)

        for i in range(n_stats):
            if (methods[i] == b'PERMANOVA') and (test_names[i] == b'pseudo-F'):
                found = True
                pmn = _build_stat('PERMANOVA', 'pseudo-F',
                                  sample_size, num_groups[i],
                                  values[i], pvalues[i], n_permutations[i])
                break

    if (not found):
        raise KeyError("PERMANOVA not found")

    return pmn


def h5permanova_dict(h5file: str) -> dict:
    """Read PERMANOVA statistical tests from a hdf5 file

    As describe in scikit-bio skbio.stats.distance.permanova.py,
    Permutational Multivariate Analysis of Variance (PERMANOVA) is a
    non-parametric method that tests whether two or more groups of objects
    are significantly different based on a categorical factor.

    Parameters
    ----------
    h5file : str
        A filepath to a hdf5 file.

    Returns
    -------
    dict[str]=pandas.Series
        Results of the statistical test, including ``test statistic`` and
        ``p-value``.

    Raises
    ------
    OSError
        If the hdf5 file is not found
    KeyError
        If the hdf5 does not have the necessary fields

    References
    ----------
    .. [1] Anderson, Marti J. "A new method for non-parametric multivariate
       analysis of variance." Austral Ecology 26.1 (2001): 32-46.
    """

    pmns = {}
    with h5py.File(h5file, "r") as f_u:
        methods = f_u['stat_methods'][:]
        test_names = f_u['stat_test_names'][:]
        grouping_names = f_u['stat_grouping_names'][:]
        values = f_u['stat_values'][:]
        pvalues = f_u['stat_pvalues'][:]
        n_permutations = f_u['stat_n_permutations'][:]
        num_groups = f_u['stat_n_groups'][:]

        sample_size = len(f_u['order'][:])

        n_stats = len(methods)

        for i in range(n_stats):
            if (methods[i] == b'PERMANOVA') and (test_names[i] == b'pseudo-F'):
                kname = grouping_names[i].decode('ascii')
                pmns[kname] = _build_stat('PERMANOVA', 'pseudo-F',
                                          sample_size, num_groups[i],
                                          values[i], pvalues[i],
                                          n_permutations[i])

    return pmns


def faith_pd(table: Union[str, Table],
             phylogeny: Union[str, TreeNode, BP]) -> pd.Series:
    """Compute Faith's Phylogenetic Diversity

    Parameters
    ----------
    table : str
        A filepath to a BIOM-Format 2.1 file.
    phylogeny : str
        A filepath to a Newick formatted tree.

    Returns
    -------
    pd.Series
        The resulting diversity calculations

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
    Faith's PD was originally described in [1]_. The implementation used here
    was originally described in [2]_.

    References
    ----------
    .. [1] Faith DP (1992) Conservation evaluation and phylogenetic diversity.
       Biol Conserv 61:1–10
    .. [2] Armstrong G, et al (2021) Efficient computation of Faith's
       phylogenetic diversity with applications in characterizing microbiomes.
       Genome Research 31(11):2131–2137
    """
    return _call_faith_pd(table, phylogeny)
