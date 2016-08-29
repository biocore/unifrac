# cython: boundscheck=False, wraparound=False, cdivision=True, linetrace=False
from libc.stdlib cimport malloc, free
#from cython.parallel cimport parallel, prange
cimport numpy as np
import numpy as np

from bp._bp cimport BP
from q2_su._ba cimport *

np.import_array()


cdef unicode tounicode(char* s):
    # from http://docs.cython.org/en/latest/src/tutorial/strings.html
    return s.decode('UTF-8', 'strict')


cdef print_ba(BIT_ARRAY* bitarr):
    cdef char* str_
    cdef object result

    str_ = <char*>malloc(bit_array_length(bitarr) + 1)
    result = tounicode(bit_array_to_str(bitarr, str_))
    free(str_)
    print(result)


#cpdef np.ndarray _unweighted_unifrac(object table, BP tree):
def _unweighted_unifrac(object table, BP tree):
    cdef dict otu_index
    cdef np.ndarray otu_ids
    cdef int i
    cdef int j
    cdef int k
    cdef int node
    cdef int shift_idx
    cdef int shift
    
    cdef np.double_t[:, :] agg_total
    cdef np.double_t[:, :] agg_unique
    cdef np.double_t[:, :] result
    cdef np.double_t length

    cdef BIT_ARRAY** state_arrays
    cdef BIT_ARRAY* current_state_p
    cdef BIT_ARRAY* rolled
    cdef BIT_ARRAY* total
    cdef BIT_ARRAY* unique
    cdef BitArrayStack stack
    cdef BIT_ARRAY** occupancy_map  # node -> state_array, wonder if this can be done smarter
    cdef bit_index_t* set_bits;
    cdef bit_index_t n_set_bits;

    cdef Py_ssize_t state_array_idx = 0
    cdef Py_ssize_t state_arrays_in_use = 0
    cdef Py_ssize_t n_state_arrays = 1024
    
    cdef Py_ssize_t n_samples;
    cdef int[:] obs_nz
    cdef int indices_start
    cdef object indices = table['observation/matrix/indices']
    cdef object indptr = table['observation/matrix/indptr']
    cdef np.double_t observed_bit

    n_samples = len(table['sample/ids'])
    agg_total = np.zeros((n_samples, n_samples), dtype=np.double)
    agg_unique = np.zeros((n_samples, n_samples), dtype=np.double)
    result = np.zeros((n_samples, n_samples), dtype=np.double)

    # mapping between a tip ID and its corresponding sparse vector in table
    otu_index = {}
    otu_ids = table['observation/ids'][:]
    for i in range(otu_ids.size):
        otu_index[otu_ids[i]] = i

    # reduce the tree to just the OTUs of interest
    tree = tree.shear(set(otu_ids)).collapse()

    # preallocate a bunch of state vectors, and push them into our memory 
    # management stack
    init_stack(&stack, n_state_arrays)
    state_arrays = <BIT_ARRAY**>malloc(sizeof(BIT_ARRAY*) * n_state_arrays)
    if not state_arrays:
        raise MemoryError()

    for i in range(n_state_arrays):
        state_arrays[i] = bit_array_create(n_samples)
        push(&stack, state_arrays[i])

    rolled = bit_array_create(n_samples)
    total = bit_array_create(n_samples)
    unique = bit_array_create(n_samples)
    
    # 2n on nodes as indices from a BP span closing parentheses as well
    occupancy_map = <BIT_ARRAY**>malloc(len(tree) * 2 * sizeof(BIT_ARRAY*))
    if not occupancy_map:
        raise MemoryError()


 #   with nogil, parallel():
    # local set_bits
    set_bits = <bit_index_t*>malloc(n_samples * sizeof(bit_index_t))

    #for k in prange(1, len(tree), schedule='guided'):  # do not include root as there is no parent
    for k in range(1, len(tree)):  # do not include root as there is no parent
        node = tree.postorderselect(k)
        length = <np.double_t>tree.length(node)
        
        current_bit_array_p = pop(&stack)
        bit_array_clear_all(current_bit_array_p)

        # associate this nodes index with the index of the bit array
        occupancy_map[node] = current_bit_array_p

        if tree.isleaf(node):
            _set_leaf_state(current_bit_array_p, otu_index[tree.name(node)],
                            indices, indptr)
        else:
            _set_nonleaf_state(&stack, current_bit_array_p, tree, node, 
                               occupancy_map)
        
        bit_array_copy_all(rolled, current_bit_array_p)
        for shift in range(1, n_samples):
            bit_array_cycle_right(rolled, 1)
            bit_array_or(total, current_bit_array_p, rolled)
            bit_array_xor(unique, current_bit_array_p, rolled)
           
            
            n_set_bits = bit_array_get_bits(total, n_samples - shift, set_bits, n_samples)
            for i in range(n_set_bits):
                j = set_bits[i]
                shift_idx = (j + shift) % n_samples
                agg_total[j, shift_idx] += length
                agg_unique[j, shift_idx] += (length * bit_array_get_bit(unique, j))

            #for i in range(n_samples - shift):
                #if bit_array_get_bit(total, i):
                #    shift_idx = (i + shift) % n_samples
                #    agg_total[i, shift_idx] += length
                #    agg_unique[i, shift_idx] += (length * bit_array_get_bit(unique, i))
    free(set_bits)

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            result[i, j] = agg_unique[i, j] / agg_total[i, j]
            result[j, i] = result[i, j]

    # free all the memory
    for i in range(n_state_arrays):
        bit_array_free(state_arrays[i])
    free(state_arrays)
    free(rolled)
    free(total)
    free(unique)
    free(stack.the_stack)  # note, pointers in stack freed by state_arrays
    free(occupancy_map)

    return np.asarray(result)


cdef void _set_leaf_state(BIT_ARRAY* state, 
                          np.int32_t indptr_start,
                          object indices,
                          object indptr): 
    """Set bits in a vector where a sample is nonzero

    Parameters
    ----------
    state : BIT_ARRAY*
        The bit array to modify. This array will be zero'd out entirely.
    indptr_start : np.int32_t
        The offset in indptr corresponding to the observation of interest
    indices : object (h5py.Group)
        A Group corresponding to 'observation/matrix/indices' in the h5py file.
    indptr : object (h5py.Group)
        A Group corresponding to 'observation/matrix/indptr' in the h5py file.
    """
    cdef int i
    cdef int indices_start
    cdef int indices_stop
    cdef np.int32_t[:] indices_slice 

    # TODO: can this be done outside of the GIL? Answer is yes, but easily?
    # issue is that the h5py parameters need to be represented as objects
    # unfortunately. So we can at the minimum limited exposure...

    # obtain the nonzero data positions 
    indices_start = indptr[indptr_start]
    indices_stop  = indptr[indptr_start + 1]
    indices_slice = indices[indices_start:indices_stop]

    # set the nonzero indices in the state array
    for i in range(indices_stop - indices_start):
        bit_array_set_bit(state, indices_slice[i])


cdef void _set_nonleaf_state(BitArrayStack* stack, BIT_ARRAY* node_state, 
                             BP tree, int node, BIT_ARRAY** occupancy):
    """Reduce the state of the children of node into node

    Parameters
    ----------
    ### STALE
    state_arrays : BIT_ARRAY**
        All of the state arrays
    node_state : int
        The index position in state_arrays assigned to the current node
    tree : BP
        The tree object
    node : int
        The node to reduce
    occupancy : memoryview of np.int32
        A mapping of node index -> state array
    ###

    Returns
    -------
    int
        The number of children reduced
    """
    cdef int last
    cdef int current
    cdef int dropped = 0
    cdef BIT_ARRAY* child_state

    current = tree.fchild(node)
    last = tree.lchild(node)

    # 0 is what we'd receive if we "nsibling"ed off the tree
    # TODO: I don't think we need to check that current != 0
    while current <= last and current != 0:
        if occupancy[current] == NULL:
            # in a parallel context, it might be possible for a processor
            # to evaluate a parent before the children are evaluated in the 
            # event of unevent processing
            # TODO: block or release if parallel until occupancy[current] != -1
            pass

        child_state = occupancy[current]
        bit_array_or(node_state, child_state, node_state)
        
        occupancy[current] = NULL
        push(stack, child_state)
        
        current = tree.nsibling(current)
