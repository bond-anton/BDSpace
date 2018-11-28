from cpython.array cimport array
cdef double __reduce_angle(double angle, bint center=*, bint positive=*)
cdef double[:] __reduce_angles(double[:] euler_angles, bint center=*, bint positive=*)
cpdef double vector_norm(double[:] v)
cpdef double[:] unit_vector(double[:] v)
cdef double[:] __extend_vector_dimensions(double[:] v, Py_ssize_t s)
cpdef double angles_between_vectors(double[:] v1, double[:] v2)
cdef double[:] __cartesian_to_spherical_point(double[:] xyz)
cdef array* __cartesian_to_spherical_points(double[:, :] xyz)
