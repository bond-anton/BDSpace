from BDSpace.Space cimport Space

cdef class ParametricCurve(Space):
    cdef:
        double __start
        double __stop
    cpdef double[:] x(self, double[:] t)
    cpdef double[:] y(self, double[:] t)
    cpdef double[:] z(self, double[:] t)
    cpdef double[:, :] generate_points(self, double[:] t)
    cpdef double[:, :] tangent(self, double[:] t)
    cpdef double length(self, double a=*, double b=*, double precision=*, bint print_details=*)
