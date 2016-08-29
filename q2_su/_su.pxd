from libc.stdlib cimport malloc, free
cimport numpy as np
from bp._bp cimport BP

from q2_su._ba cimport *


cdef struct s_BitArrayStack:
    BIT_ARRAY** the_stack
    Py_ssize_t position
    Py_ssize_t stack_size
ctypedef s_BitArrayStack BitArrayStack


cdef inline void init_stack(BitArrayStack* stack, Py_ssize_t n):
    stack.the_stack = <BIT_ARRAY**>malloc(sizeof(BIT_ARRAY*) * n)
    if not stack.the_stack:
        raise MemoryError()

    stack.position = n
    stack.stack_size = n


cdef inline void push(BitArrayStack* stack, BIT_ARRAY* item):
    if stack.position == 0:
        # resize
        raise MemoryError()
    stack.position -= 1
    stack.the_stack[stack.position] = item

    
cdef inline BIT_ARRAY* pop(BitArrayStack* stack):
    cdef BIT_ARRAY* item

    assert stack.position <= stack.stack_size
    # TODO: barrier for parallel
    item = stack.the_stack[stack.position]
    stack.position += 1
    return item


cdef void _set_leaf_state(BIT_ARRAY* state, 
                          np.int32_t indptr_start,
                          object indices,
                          object indptr)

cdef void _set_nonleaf_state(BitArrayStack* stack, BIT_ARRAY* node_state, 
                            BP tree, int node, BIT_ARRAY** occupancy)
