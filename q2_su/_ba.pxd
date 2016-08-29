cdef extern from "<inttypes.h>":
    ctypedef unsigned int uint64_t

cdef extern from "../BitArray/bit_array.h":
    struct BIT_ARRAY:
        pass
    ctypedef uint64_t bit_index_t

    # allocations
    BIT_ARRAY* bit_array_create(bit_index_t nbits)
    void bit_array_free(BIT_ARRAY* bitarray)
    BIT_ARRAY* bit_array_clone(const BIT_ARRAY* bitarr)
    void bit_array_copy_all(BIT_ARRAY* dst, const BIT_ARRAY* src)
   
    # utility
    char* bit_array_to_str(const BIT_ARRAY* bitarr, char* str)
    bit_index_t bit_array_length(const BIT_ARRAY* bit_arr)
    bit_index_t bit_array_num_bits_set(const BIT_ARRAY* bitarr)

    # bit juggling
    void bit_array_set_bit(BIT_ARRAY* bitarr, bit_index_t b) nogil
    void bit_array_toggle_bit(BIT_ARRAY* bitarr, bit_index_t b) nogil
    void bit_array_set_bits(BIT_ARRAY* bitarr, size_t n, ...) nogil
    char bit_array_get_bit(const BIT_ARRAY* bitarr, bit_index_t b) nogil
    bit_index_t bit_array_get_bits(const BIT_ARRAY* bitarr, bit_index_t end, bit_index_t* dst, bit_index_t n)
    void bit_array_clear_bit(BIT_ARRAY* bitarr, bit_index_t b)
    void bit_array_clear_all(BIT_ARRAY* bitarr) nogil

    # logical operations
    void bit_array_and(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2) nogil
    void bit_array_or(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2) nogil
    void bit_array_xor(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2) nogil
    void bit_array_not(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2) nogil

    # cyclic shifting
    void bit_array_cycle_right(BIT_ARRAY* bitarr, bit_index_t dist)
    void bit_array_cycle_left (BIT_ARRAY* bitarr, bit_index_t dist)

    # comparisons
    # 
    # (from bit_array.h) comparison functions return
    #   1 iff bitarr1 > bitarr2
    #   0 iff bitarr1 == bitarr2
    #  -1 iff bitarr1 < bitarr2
    int bit_array_cmp(const BIT_ARRAY* bitarr1, const BIT_ARRAY* bitarr2)
