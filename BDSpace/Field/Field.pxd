from BDSpace.Space cimport Space

cdef class Field(Space):
    cdef:
        str __type
