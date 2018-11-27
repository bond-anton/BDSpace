cdef double __reduce_angle(double angle, bint center=*, bint positive=*)
cdef double[:] __reduce_angles(double[:] euler_angles, bint center=*, bint positive=*)
cpdef double vector_norm(double[:] v)
cpdef double[:] unit_vector(double[:] v)
