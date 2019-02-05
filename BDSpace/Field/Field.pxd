from BDSpace.Space cimport Space

cdef class Field(Space):
    cdef:
        str __type

    cdef double[:] __points_scalar(self, double[:, :] xyz, double value)
    cdef double[:, :] __points_vector(self, double[:, :] xyz, double[:] value)
    cpdef double[:] scalar_field(self, double[:, :] xyz)
    cpdef double[:, :] vector_field(self, double[:, :] xyz)


cdef class ConstantScalarConservativeField(Field):
    cdef:
        double __potential


cdef class ConstantVectorConservativeField(Field):
    cdef:
        double[:] __potential
