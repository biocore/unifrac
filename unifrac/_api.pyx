# cython: boundscheck=False
# distutils: language = c++
import skbio
import numpy as np
cimport numpy as np
import bp
import pandas as pd
from cython.parallel import prange
from libcpp.string cimport string
from libcpp.vector cimport vector

# 
# Functions that compute Unifrac and return a memory object
#

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
        np.ndarray[np.double_t, ndim=1] numpy_arr
        double *cf
        int i
        bytes biom_py_bytes
        bytes tree_py_bytes
        bytes met_py_bytes
        char* biom_c_string
        char* tree_c_string
        char* met_c_string
        list ids

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
        else:
            raise Exception("Unknown Error: {}".format(status))

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

    if status != okay:
        if status == tree_missing:
            raise IOError("Tree file not found.")
        elif status == table_missing:
            raise IOError("Table file not found.")
        elif status == output_error:
            raise IOError("Failed to generate outputfile.")
        elif status == table_empty:
            raise ValueError("Table file is empty.")
        elif status == table_and_tree_do_not_overlap:
            raise ValueError("The table does not appear to be completely "
                             "represented by the phylogeny.")
        elif status == unknown_method:
            raise ValueError("Unknown method.")
        else:
            raise Exception("Unknown Error: {}".format(status))

    return out_filename


def cppbptree_constructor(object tree):
    cdef:
        object tree_as_bp
        np.ndarray[np.uint8_t, ndim=1] structure
        np.ndarray[np.float64_t, ndim=1] lengths
        np.ndarray[object, ndim=1] names
        int i, size

    if isinstance(tree, skbio.TreeNode):
        tree_as_bp = bp.from_skbio_treenode(tree)
    else:
        tree_as_bp = tree

    return cppbptree(tree_as_bp)


cdef class cppbptree:
    cdef BPTree* obj

    def __cinit__(self, object tree):
        cdef:
            vector[bool] input_structure
            vector[double] input_lengths
            vector[string] input_names
            int i, length
            unicode name

        # TODO: enable direct access to names and lengths from BP
        length = tree.B.size

        input_structure.resize(length)
        input_lengths.resize(length)
        input_names.resize(length)
   
        for i in range(length):
            input_structure[i] = tree.B[i]
            input_lengths[i] = tree.length(i)
            
            name = tree.name(i)
            if name is not None:
                input_names[i] = name.encode('utf8')
    
        self.obj = new BPTree(input_structure, input_lengths, input_names)

    def __del__(self):
        del self.obj


cdef class cppbiom:
    cdef biom* obj

    def __cinit__(self, object table):
        cdef: 
            vector[string] obs_ids
            vector[string] samp_ids
            vector[uint32_t] indices
            vector[uint32_t] indptr
            vector[double] data

            np.ndarray[object, ndim=1] table_obs_ids
            np.ndarray[object, ndim=1] table_samp_ids
            np.ndarray[np.int32_t, ndim=1] table_indices 
            np.ndarray[np.int32_t, ndim=1] table_indptr 
            np.ndarray[np.float64_t, ndim=1] table_data
           
            int i, nsamples, nobs, nnz, nindptr
            object matrix_data = table.matrix_data.tocsr()

        ### I think these need to be expressed as memoryviews for prange ###
        table_obs_ids = table.ids(axis='observation')
        table_samp_ids = table.ids(axis='sample')
        table_indices = matrix_data.indices
        table_indptr = matrix_data.indptr
        table_data = matrix_data.data

        nsamples = table_samp_ids.size
        nobs = table_obs_ids.size
        nnz = table_data.size
        nindptr = table_indptr.size
        
        obs_ids.resize(nobs)
        samp_ids.resize(nsamples)
        indices.resize(nnz)
        indptr.resize(nindptr)
        data.resize(nnz)

        # cannot use prange for strings as the GIL is required for indexing
        # arrays of object dtype
        for i in range(nsamples):
            samp_ids[i] = table_samp_ids[i].encode('utf8')

        for i in range(nobs):
            obs_ids[i] = table_obs_ids[i].encode('utf8')
        
        for i in prange(nindptr, nogil=True, schedule='static'):
            indptr[i] = table_indptr[i]

        for i in prange(nnz, nogil=True, schedule='static'):
            indices[i] = table_indices[i]
            data[i] = table_data[i]

        self.obj = new biom(obs_ids, samp_ids, indices, indptr, data)

    def __del__(self):
        del self.obj
