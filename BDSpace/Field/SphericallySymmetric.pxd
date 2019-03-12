from .Field cimport Field


cdef class SphericallySymmetric(Field):
    cdef:
        double __r

    cdef inline double __get_r(self, double[:] xyz) nogil
    cdef double scalar_field_r_law(self, double r) nogil
    cpdef double scalar_field_r_point(self, double r)
    cpdef double[:] scalar_field_r(self, double[:] r)
    cdef double vector_field_r_law(self, double r) nogil
    cpdef double vector_field_r_point(self, double r)
    cpdef double[:] vector_field_r(self, double[:] r)

cdef class HyperbolicPotentialSphericalConservativeField(SphericallySymmetric):
    cdef:
        double __a
