# cython: boundscheck=False
import skbio
import numpy as np
cimport numpy as np
import bp
import pandas as pd
from cython.parallel import prange
from libc.stdlib cimport malloc, free
from libc.string cimport strcpy


def check_status(compute_status status):
    if status != okay:
        if status == tree_missing:
            raise IOError("Tree file not found.")
        elif status == table_missing:
            raise IOError("Table file not found.")
        elif status == table_empty:
            raise ValueError("Table file is empty.")
        elif status == table_and_tree_do_not_overlap:
            raise ValueError("The table does not appear to be completely "
                             "represented by the phylogeny.")
        elif status == unknown_method:
            raise ValueError("Unknown method.")
        else:
            raise Exception("Unknown Error: {}".format(status))

# 
# Functions that compute Unifrac and return a memory object
#
def ssu_inmem(object table, object tree,
              str unifrac_method, bool variance_adjust, double alpha,
              bool bypass_tips, unsigned int threads):
    """Execute a call to Strided State UniFrac via the direct API

    Parameters
    ----------
    table : biom.Table
        An instance of a biom.Table object
    tree : bp.BP or skbio.TreeNode
        A phylogeny corresponding to the table
    unifrac_method : str
        The requested UniFrac method, one of {unweighted,
        weighted_normalized, weighted_unnormalized, generalized,
        unweighted_fp32, weighted_normalized_fp32, 
        weighted_unnormalized_fp32, generalized_fp32}
    variance_adjust : bool
        Whether to perform Variance Adjusted UniFrac
    alpha : float
        The value of alpha for Generalized UniFrac; only applies to
        Generalized UniFraca
    bypass_tips : bool
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    threads : int
        The number of threads to use.

    Returns
    -------
    skbio.DistanceMatrix
        The resulting distance matrix

    Raises
    ------
    ValueError
        If the table is empty
        If the table is not completely represented by the phylogeny
        If an unknown method is requested.
    Exception
        If an unkown error is experienced
    """
    cdef:
        mat *result;
        compute_status status;
        bytes met_py_bytes
        char* met_c_string
        support_biom* inmem_biom
        support_bptree* inmem_tree

    inmem_biom = construct_support_biom(table)
    inmem_tree = construct_support_bptree(tree)

    met_py_bytes = unifrac_method.encode()
    met_c_string = met_py_bytes

    status = one_off_inmem(inmem_biom,
                           inmem_tree,
                           met_c_string,
                           variance_adjust,
                           alpha,
                           bypass_tips,
                           threads,
                           &result)
    check_status(status)
   
    destroy_support_biom(inmem_biom)
    destroy_support_bptree(inmem_tree)

    return result_to_skbio_distance_matrix(result)


def ssu(str biom_filename, str tree_filename,
        str unifrac_method, bool variance_adjust, double alpha,
        bool bypass_tips, unsigned int threads):
    """Execute a call to Strided State UniFrac via the direct API

    Parameters
    ----------
    biom_filename : str
        A filepath to a BIOM 2.1 formatted table (HDF5)
    tree_filename : str
        A filepath to a Newick formatted tree
    unifrac_method : str
        The requested UniFrac method, one of {unweighted,
        weighted_normalized, weighted_unnormalized, generalized,
        unweighted_fp32, weighted_normalized_fp32, 
        weighted_unnormalized_fp32, generalized_fp32}
    variance_adjust : bool
        Whether to perform Variance Adjusted UniFrac
    alpha : float
        The value of alpha for Generalized UniFrac; only applies to
        Generalized UniFraca
    bypass_tips : bool
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    threads : int
        The number of threads to use.

    Returns
    -------
    skbio.DistanceMatrix
        The resulting distance matrix

    Raises
    ------
    IOError
        If the tree file is not found
        If the table is not found
    ValueError
        If the table is empty
        If the table is not completely represented by the phylogeny
        If an unknown method is requested.
    Exception
        If an unkown error is experienced
    """
    cdef:
        mat *result;
        compute_status status;
        bytes biom_py_bytes
        bytes tree_py_bytes
        bytes met_py_bytes
        char* biom_c_string
        char* tree_c_string
        char* met_c_string

    biom_py_bytes = biom_filename.encode()
    tree_py_bytes = tree_filename.encode()
    met_py_bytes = unifrac_method.encode()
    biom_c_string = biom_py_bytes
    tree_c_string = tree_py_bytes
    met_c_string = met_py_bytes

    status = one_off(biom_c_string,
                     tree_c_string,
                     met_c_string,
                     variance_adjust,
                     alpha,
                     bypass_tips,
                     threads,
                     &result)
    check_status(status)
    
    return result_to_skbio_distance_matrix(result)


cdef object result_to_skbio_distance_matrix(mat *result):
    cdef:
        np.ndarray[np.double_t, ndim=1] numpy_arr
        list ids
        int i

    ids = []
    numpy_arr = np.zeros(result.cf_size, dtype=np.double)
    numpy_arr[:] = <np.double_t[:result.cf_size]> result.condensed_form

    for i in range(result.n_samples):
        ids.append(result.sample_ids[i].decode('utf-8'))

    destroy_mat(&result)

    return skbio.DistanceMatrix(numpy_arr, ids)


def faith_pd(str biom_filename, str tree_filename):
    """Execute a call to the Stacked Faith API in the UniFrac package

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
    ValueError
        If the table is empty
        If the table is not completely represented by the phylogeny
    Exception
        If an unkown error is experienced
    """
    cdef:
        results_vec *result;
        compute_status status;
        np.ndarray[np.double_t, ndim=1] numpy_arr
        bytes biom_py_bytes
        bytes tree_py_bytes
        char* biom_c_string
        char* tree_c_string
        list ids


    biom_py_bytes = biom_filename.encode()
    tree_py_bytes = tree_filename.encode()
    biom_c_string = biom_py_bytes
    tree_c_string = tree_py_bytes

    status = faith_pd_one_off(biom_c_string, tree_c_string, &result)
    check_status(status)

    numpy_arr = np.zeros(result.n_samples, dtype=np.double)
    numpy_arr[:] = <np.double_t[:result.n_samples]> result.values

    ids = []
    for i in range(result.n_samples):
        ids.append(result.sample_ids[i].decode('utf-8'))

    faith_pd_series = pd.Series(numpy_arr, index=ids)
    faith_pd_series.rename("faith_pd", inplace=True)

    destroy_results_vec(&result)

    return faith_pd_series

# 
# Functions that compute Unifrac and write into a file
#

def ssu_to_file(str biom_filename, str tree_filename, str out_filename,
                str unifrac_method, bool variance_adjust, double alpha,
                bool bypass_tips, unsigned int threads, str format,
                unsigned int pcoa_dims, str buf_dirname):
    """Execute a call to Strided State UniFrac to file via the direct API

    Parameters
    ----------
    biom_filename : str
        A filepath to a BIOM 2.1 formatted table (HDF5)
    tree_filename : str
        A filepath to a Newick formatted tree
    out_filename : str
        A filepath to the output file.
    unifrac_method : str
        The requested UniFrac method, one of {unweighted,
        weighted_normalized, weighted_unnormalized, generalized,
        unweighted_fp32, weighted_normalized_fp32, 
        weighted_unnormalized_fp32, generalized_fp32}
    variance_adjust : bool
        Whether to perform Variance Adjusted UniFrac
    alpha : float
        The value of alpha for Generalized UniFrac; only applies to
        Generalized UniFraca
    bypass_tips : bool
        Bypass the tips of the tree in the computation. This reduces compute
        by about 50%, but is an approximation.
    threads : int
        The number of threads to use.
    format : str
        Onput format to use; one of {hdf5, hdf5_fp32, hdf5_fp64}
    pcoa_dims : int
        Number of dimension for PCoA, if 0, no PCoA is computed.
    out_filename : str
        If using a disk buffer for saving memory is desired, a dirpath.

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
        If the table is empty
        If the table is not completely represented by the phylogeny
        If an unknown method is requested.
    Exception
        If an unkown error is experienced
    """
    cdef:
        compute_status status;
        int i
        bytes biom_py_bytes
        bytes tree_py_bytes
        bytes out_py_bytes
        bytes met_py_bytes
        bytes format_py_bytes
        bytes dirbuf_py_bytes
        char* biom_c_string
        char* tree_c_string
        char* out_c_string
        char* met_c_string
        char* format_c_string
        char* dirbuf_c_string
        list ids

    biom_py_bytes = biom_filename.encode()
    tree_py_bytes = tree_filename.encode()
    out_py_bytes = out_filename.encode()
    met_py_bytes = unifrac_method.encode()
    format_py_bytes = format.encode()
    dirbuf_py_bytes = buf_dirname.encode()
    biom_c_string = biom_py_bytes
    tree_c_string = tree_py_bytes
    out_c_string = out_py_bytes
    met_c_string = met_py_bytes
    format_c_string = format_py_bytes
    dirbuf_c_string = dirbuf_py_bytes

    status = unifrac_to_file(biom_c_string, tree_c_string, out_c_string,
                             met_c_string,
                             variance_adjust, alpha, bypass_tips,
                             threads, format_c_string, 
                             pcoa_dims, dirbuf_c_string)
    check_status(status)

    return out_filename


cdef support_biom* construct_support_biom(object table):
    cdef: 
        char** obs_ids
        char** sample_ids
        int32_t* indices
        int32_t* indptr
        double* data
        support_biom* sp_biom

        object matrix_data = table.matrix_data.tocsr()
        int i, nsamples, nobs, nnz, nindptr, n

        np.ndarray[object, ndim=1] table_obs_ids
        np.ndarray[object, ndim=1] table_samp_ids
        np.ndarray[np.int32_t, ndim=1] table_indices 
        np.ndarray[np.int32_t, ndim=1] table_indptr 
        np.ndarray[np.float64_t, ndim=1] table_data

    table_obs_ids = table.ids(axis='observation')
    table_samp_ids = table.ids(axis='sample')
    table_indices = matrix_data.indices
    table_indptr = matrix_data.indptr
    table_data = matrix_data.data

    nsamples = table_samp_ids.size
    nobs = table_obs_ids.size
    nnz = table_data.size
    nindptr = table_indptr.size
   
    obs_ids = <char**>malloc(nobs * sizeof(char*))
    if not obs_ids:
        return NULL

    sample_ids = <char**>malloc(nsamples * sizeof(char*))
    if not sample_ids:
        return NULL

    indices = <int32_t*>malloc(nnz * sizeof(int32_t))
    if not indices:
        return NULL

    indptr = <int32_t*>malloc(nindptr * sizeof(int32_t))
    if not indptr:
        return NULL

    data = <double*>malloc(nnz * sizeof(double))
    if not data:
        return NULL

    # cannot use prange for strings as the GIL is required for indexing
    # arrays of object dtype
    for i in range(nsamples):
        n = len(table_samp_ids[i]) + 1  # for \0
        sample_ids[i] = <char*>malloc(n * sizeof(char))
        if not sample_ids[i]:
            return NULL
        strcpy(sample_ids[i], table_samp_ids[i].encode('ascii'))

    for i in range(nobs):
        n = len(table_obs_ids[i]) + 1  # for \0
        obs_ids[i] = <char*>malloc(n * sizeof(char))
        if not obs_ids[i]:
            return NULL
        strcpy(obs_ids[i], table_obs_ids[i].encode('ascii'))
    
    for i in prange(nindptr, nogil=True, schedule='static'):
        indptr[i] = table_indptr[i]

    for i in prange(nnz, nogil=True, schedule='static'):
        indices[i] = table_indices[i]
        data[i] = table_data[i]

    sp_biom = <support_biom*>malloc(sizeof(support_biom))
    if not sp_biom:
        return NULL

    sp_biom.obs_ids = obs_ids
    sp_biom.sample_ids = sample_ids
    sp_biom.indices = indices
    sp_biom.indptr = indptr
    sp_biom.data = data
    sp_biom.n_obs = nobs
    sp_biom.n_samples = nsamples
    sp_biom.nnz = nnz

    return sp_biom

cdef void destroy_support_biom(support_biom *sp_biom):
    cdef:
        int i

    for i in range(sp_biom.n_obs):
        free(sp_biom.obs_ids[i])
    free(sp_biom.obs_ids)

    for i in range(sp_biom.n_samples):
        free(sp_biom.sample_ids[i])
    free(sp_biom.sample_ids)

    free(sp_biom.indices)
    free(sp_biom.indptr)
    free(sp_biom.data)
    free(sp_biom)


cdef support_bptree* construct_support_bptree(object tree):
    cdef:
        bool* structure
        double* lengths
        char** names
        object tree_as_bp
        int length, i, n
        support_bptree* sp_bptree

    if isinstance(tree, skbio.TreeNode):
        tree_as_bp = bp.from_skbio_treenode(tree)
    else:
        tree_as_bp = tree

    length = tree_as_bp.B.size

    # malloc these things
    structure = <bool*>malloc(length * sizeof(bool))
    if not structure:
        return NULL

    lengths = <double*>malloc(length * sizeof(double))
    if not lengths:
        return NULL

    names = <char**>malloc(length * sizeof(char*))
    if not names:
        return NULL

    # cannot prange as we need to do attribute access
    # TODO: expose these structures directly from BP
    for i in range(length):
        structure[i] = tree_as_bp.B[i]
        lengths[i] = tree_as_bp.length(i)
        
        name = tree_as_bp.name(i)
        if name is None:
            names[i] = <char*>malloc(sizeof(char))
            if not names[i]:
                return NULL
            names[i][0] = b'\0'
        else:
            n = len(name) + 1  # for \0
            
            names[i] = <char*>malloc(n * sizeof(char))
            if not names[i]:
                return NULL

            strcpy(names[i], name.encode('ascii'))
  
    sp_bptree = <support_bptree*>malloc(sizeof(support_bptree))
    if not sp_bptree:
        return NULL

    sp_bptree.structure = structure
    sp_bptree.lengths = lengths
    sp_bptree.names = names
    sp_bptree.n_parens = length
    
    return sp_bptree


cdef void destroy_support_bptree(support_bptree* sp_bptree):
    cdef:
        int i

    for i in range(sp_bptree.n_parens):
        free(sp_bptree.names[i])

    free(sp_bptree.names)
    free(sp_bptree.lengths)
    free(sp_bptree.structure)
    free(sp_bptree)
