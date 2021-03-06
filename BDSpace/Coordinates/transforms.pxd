cdef double __reduce_angle(double angle, bint center=*, bint positive=*) nogil
cpdef reduce_angle(double angle, bint keep_sign=*)
cpdef double[:] reduce_angles(double[:] euler_angles, bint keep_sign=*)
cpdef double vector_norm(double[:] v)
cpdef double[:] unit_vector(double[:] v)
cdef double[:] __extend_vector_dimensions(double[:] v, Py_ssize_t s)
cpdef double angles_between_vectors(double[:] v1, double[:] v2)
cpdef double[:] cartesian_to_spherical_point(double[:] xyz)
cpdef double[:, :] cartesian_to_spherical(double[:, :] xyz)
cpdef double[:] spherical_to_cartesian_point(double[:] r_theta_phi)
cpdef double[:, :] spherical_to_cartesian(double[:, :] r_theta_phi)
cpdef double[:] invert_spherical_point(double[:] r_theta_phi)
cpdef double[:, :] invert_spherical(double[:, :] r_theta_phi)
cpdef double[:] cartesian_to_cylindrical_point(double[:] xyz)
cpdef double[:, :] cartesian_to_cylindrical(double[:, :] xyz)
cpdef double[:] cylindrical_to_cartesian_point(double[:] rho_phi_z)
cpdef double[:, :] cylindrical_to_cartesian(double[:, :] rho_phi_z)
cpdef double[:] spherical_to_cylindrical_point(double[:] r_theta_phi)
cpdef double[:, :] spherical_to_cylindrical(double[:, :] r_theta_phi)
cpdef double[:] cylindrical_to_spherical_point(double[:] rho_phi_z)
cpdef double[:, :] cylindrical_to_spherical(double[:, :] rho_phi_z)
