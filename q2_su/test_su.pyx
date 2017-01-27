from unittest import TestCase
from libc.stdlib cimport malloc, free

cimport numpy as np
import numpy as np
import numpy.testing as npt
import scipy.sparse as ss
import h5py
from bp import parse_newick, BP

from q2_su._su cimport (_set_leaf_state, _set_nonleaf_state, BitArrayStack,
                        init_stack, push, pop) 
from q2_su._su import _unweighted_unifrac
from q2_su._ba cimport (bit_array_create, bit_array_free, bit_array_set_bit,
                        bit_array_get_bit, bit_array_cmp, BIT_ARRAY, 
                        bit_array_set_bits, bit_array_length, bit_array_to_str)


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

class StateUniFracTests(TestCase):
    def test_set_leaf_state(self):
        cdef BIT_ARRAY* obs = bit_array_create(32)
        cdef BIT_ARRAY* exp = bit_array_create(32)
        cdef np.int32_t[:] indices
        cdef np.int32_t[:] indptr
        
        indices = np.arange(0, 64, 2, dtype=np.int32)
        indptr = np.array([0, 4, 10, 11, 33], dtype=np.int32)

        bit_array_set_bit(exp, 8)
        bit_array_set_bit(exp, 10)
        bit_array_set_bit(exp, 12)
        bit_array_set_bit(exp, 14)
        bit_array_set_bit(exp, 16)
        bit_array_set_bit(exp, 18)
        
        _set_leaf_state(obs, 1, indices, indptr)
        self.assertEqual(bit_array_cmp(obs, exp), 0)
       
        bit_array_free(obs)
        bit_array_free(exp)

    def test_set_nonleaf_state(self):
        cdef BitArrayStack stack
        cdef BIT_ARRAY** occupancy
        cdef BIT_ARRAY* a_state
        cdef BIT_ARRAY* b_state
        cdef BIT_ARRAY* c_state

        cdef int i
        cdef int node
        tree = parse_newick("(((a,b)c,(d,e)f)g,h)root;")
        # TODO: add find_by_name method for BP
        for i in range(len(tree)):
            if tree.name(i) == 'c':
                node = i
                break

        cdef BIT_ARRAY* exp
        cdef int exp_n_reduced = 2
        cdef int obs_n_reduced = 0

        init_stack(&stack, 3)
        push(&stack, bit_array_create(4))
        push(&stack, bit_array_create(4))
        push(&stack, bit_array_create(4))
        exp = bit_array_create(4)

        bit_array_set_bits(stack.the_stack[0], 2, 0, 2)
        bit_array_set_bits(stack.the_stack[1], 2, 0, 3)
        bit_array_set_bits(exp, 3, 0, 2, 3)
        
        a_state = pop(&stack) # node a
        b_state = pop(&stack) # node b
        c_state = pop(&stack) # node c
       
        occupancy = <BIT_ARRAY**>malloc(len(tree) * 2 * sizeof(BIT_ARRAY*))
        occupancy[2] = c_state
        occupancy[3] = a_state
        occupancy[4] = a_state
        occupancy[5] = b_state
        occupancy[6] = b_state
        occupancy[7] = c_state
       
        _set_nonleaf_state(&stack, c_state, tree, node, occupancy)
        self.assertEqual(bit_array_cmp(c_state, exp), 0)
        
        bit_array_free(a_state)
        bit_array_free(b_state)
        bit_array_free(c_state)
        free(stack.the_stack)
        bit_array_free(exp)
        free(occupancy)

    def test_unweighted_unifrac(self):
        # expected result pulled from skbio:
        # https://github.com/biocore/scikit-bio/blob/0.5.0/skbio/diversity/tests/test_driver.py#L514
        tree = parse_newick('((O1:0.25, O2:0.50):0.25, O3:0.75)root;')
        table_data = ss.csc_matrix([[1, 5],
                                    [2, 3], 
                                    [0, 1]]).astype(float)
        sids = list('ABC')
        oids = ['O1', 'O2']
        
        # create an in-memory HDF5 file so we don't have to muck with tempfiles
        table = h5py.File('ignored', driver='core', backing_store=False)
        table['observation/ids'] = oids
        table['sample/ids'] = sids
        table['observation/matrix/data'] = table_data.data
        table['observation/matrix/indices'] = table_data.indices
        table['observation/matrix/indptr'] = table_data.indptr
        
        expected_data = np.asarray([[0.0, 0.0, 0.25],
                                    [0.0, 0.0, 0.25],
                                    [0.25, 0.25, 0.0]])

        obs = _unweighted_unifrac(table, tree)
        npt.assert_almost_equal(obs, expected_data)
