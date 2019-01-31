from BDSpace.Space cimport Space

cdef class Field(Space):
    cdef:
        str __type

    cpdef double[:] points_scalar(self, double[:, :] xyz, double value)
    cpdef double[:, :] points_vector(self, double[:, :] xyz, double[:] value)
