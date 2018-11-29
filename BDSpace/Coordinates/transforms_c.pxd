cdef double __reduce_angle(double angle, bint center=*, bint positive=*)
cdef double[:] __reduce_angles(double[:] euler_angles, bint center=*, bint positive=*)
cpdef double vector_norm(double[:] v)
cpdef double[:] unit_vector(double[:] v)
cdef double[:] __extend_vector_dimensions(double[:] v, Py_ssize_t s)
cpdef double angles_between_vectors(double[:] v1, double[:] v2)
cdef double[:] __cartesian_to_spherical_point(double[:] xyz)
cdef double[:] __cartesian_to_spherical_points(double[:] xyz)
cdef double[:] __spherical_to_cartesian_point(double[:] r_theta_phi)
cdef double[:] __spherical_to_cartesian_points(double[:] r_theta_phi)
