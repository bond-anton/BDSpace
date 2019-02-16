import numbers
import numpy as np

from cython import boundscheck, wraparound
from cython.parallel import prange

from cpython.array cimport array, clone
from libc.math cimport fabs, sin, cos, atan2, acos, sqrt, M_PI


cdef double __reduce_angle(double angle, bint center=True, bint positive=False) nogil:
    """
    Adjusts angle to be in the range [-2*pi; 2*pi]
    :param angle: angle or array-like of input angle
    :param center: if True (default) adjust angle to be within [-pi; pi]
    :param half: if True (default False) adjust angle to be within [0; pi]
    :return: reduced angle or array of reduced angles
    """
    cdef:
        double reduced_angle
    if fabs(angle) > 2 * M_PI:
        reduced_angle = angle - 2 * M_PI * (angle // (2 * M_PI))
    else:
        reduced_angle = angle
    if center:
        if reduced_angle < -M_PI:
            reduced_angle += 2 * M_PI
        elif reduced_angle > M_PI:
            reduced_angle -= 2 * M_PI
    if positive:
        if reduced_angle < 0:
            reduced_angle += 2 * M_PI
    return reduced_angle


cpdef reduce_angle(double angle, bint keep_sign=False):
    return __reduce_angle(angle, center=False, positive=not keep_sign)


@boundscheck(False)
@wraparound(False)
cpdef double[:] reduce_angles(double[:] angles, bint keep_sign=False):
    cdef:
        int i
        Py_ssize_t s = angles.size
        array[double] reduced_angles, template = array('d')
    reduced_angles = clone(template, s, zero=False)
    with nogil:
        for i in prange(s):
            reduced_angles[i] = __reduce_angle(angles[i], center=False, positive=not keep_sign)
    return reduced_angles

@boundscheck(False)
@wraparound(False)
cpdef double vector_norm(double[:] v):
    cdef:
        int i
        double result = 0.0
    for i in range(v.size):
        result += v[i] * v[i]
    return sqrt(result)


@boundscheck(False)
@wraparound(False)
cpdef double[:] unit_vector(double[:] v):
    """
    returns unit vector for given vector v
    :param v: input vector
    :return: unit vector u parallel v of length 1
    """
    cdef:
        int i
        Py_ssize_t s = v.size
        double length = vector_norm(v)
        array[double] result, template = array('d')
    result = clone(template, s, zero=False)

    if length == 0:
        for i in range(s):
            result[i] = 0.0
    else:
        for i in range(s):
            result[i] = v[i] / length
    return result


@boundscheck(False)
@wraparound(False)
cdef double[:] __extend_vector_dimensions(double[:] v, Py_ssize_t s):
    cdef:
        Py_ssize_t s1 = v.size
        array[double] result, template = array('d')
        int i
    result = clone(template, s, False)
    for i in range(s):
        if i < s1:
            result[i] = v[i]
        else:
            result[i] = 0.0
    return result


@boundscheck(False)
@wraparound(False)
cpdef double angles_between_vectors(double[:] v1, double[:] v2):
    cdef:
        int i
        Py_ssize_t s1 = v1.size, s2 = v2.size
        Py_ssize_t s = s1
        double [:] v1_u = unit_vector(v1)
        double [:] v2_u = unit_vector(v2)
        double cos_angle = 0
    if s1 > s2:
        v2_u = __extend_vector_dimensions(v2_u, s)
    elif s1 < s2:
        s = s2
        v1_u = __extend_vector_dimensions(v1_u, s)
    for i in range(s):
        cos_angle += v1_u[i] * v2_u[i]
    return acos(cos_angle)


@boundscheck(False)
@wraparound(False)
cpdef double[:] cartesian_to_spherical_point(double[:] xyz):
    cdef:
        Py_ssize_t s = xyz.shape[0]
        array[double] r_theta_phi, template = array('d')
        double xy
    if s < 3:
        xyz = __extend_vector_dimensions(xyz, 3)
    r_theta_phi = clone(template, 3, False)
    xy = xyz[0]*xyz[0] + xyz[1]*xyz[1]
    r_theta_phi[0] = sqrt(xy + xyz[2]*xyz[2])
    r_theta_phi[1] = atan2(sqrt(xy), xyz[2])
    r_theta_phi[2] = __reduce_angle(atan2(xyz[1], xyz[0]), center=False, positive=True)
    return r_theta_phi


@boundscheck(False)
@wraparound(False)
cpdef double[:, :] cartesian_to_spherical(double[:, :] xyz):
    cdef:
        unsigned int i, s = xyz.shape[0]
        double[:, :] r_theta_phi = np.empty((s, 3), dtype=np.double)
        # double xy
    for i in range(s):
        r_theta_phi[i] = cartesian_to_spherical_point(xyz[i])
        # xy = xyz[i]*xyz[i] + xyz[i + 1]*xyz[i + 1]
        # r_theta_phi[i] = sqrt(xy + xyz[i + 2]*xyz[i + 2])
        # r_theta_phi[i + 1] = atan2(sqrt(xy), xyz[i + 2])
        # r_theta_phi[i + 2] = __reduce_angle(atan2(xyz[i + 1], xyz[i]), center=False, positive=True)
    return r_theta_phi


@boundscheck(False)
@wraparound(False)
cpdef double[:] spherical_to_cartesian_point(double[:] r_theta_phi):
    cdef:
        Py_ssize_t s = r_theta_phi.shape[0]
        array[double] xyz, template = array('d')
        double xy
    if s < 3:
        r_theta_phi = __extend_vector_dimensions(r_theta_phi, 3)
    xyz = clone(template, 3, False)
    r_theta_phi[2] = __reduce_angle(r_theta_phi[2], center=False, positive=True)
    xy = r_theta_phi[0] * sin(r_theta_phi[1])
    xyz[0] = xy * cos(r_theta_phi[2])
    xyz[1] = xy * sin(r_theta_phi[2])
    xyz[2] = r_theta_phi[0] * cos(r_theta_phi[1])
    return xyz


@boundscheck(False)
@wraparound(False)
cpdef double[:, :] spherical_to_cartesian(double[:, :] r_theta_phi):
    cdef:
        unsigned int i, s = r_theta_phi.shape[0]
        double[:, :] xyz = np.empty((s, 3), dtype=np.double)
    for i in range(s):
        xyz[i] = spherical_to_cartesian_point(r_theta_phi[i])
    return xyz


@boundscheck(False)
@wraparound(False)
cpdef double[:] cartesian_to_cylindrical_point(double[:] xyz):
    cdef:
        Py_ssize_t s = xyz.shape[0]
        array[double] rho_phi_z, template = array('d')
    if s < 3:
        xyz = __extend_vector_dimensions(xyz, 3)
    rho_phi_z = clone(template, 3, False)
    rho_phi_z[0] = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1])
    rho_phi_z[1] = atan2(xyz[1], xyz[0])
    rho_phi_z[2] = xyz[2]
    return rho_phi_z


@boundscheck(False)
@wraparound(False)
cpdef double[:, :] cartesian_to_cylindrical(double[:, :] xyz):
    cdef:
        unsigned int i, s = xyz.shape[0]
        double[:, :] rho_phi_z = np.empty((s, 3), dtype=np.double)
    for i in range(s):
        rho_phi_z[i] = cartesian_to_cylindrical_point(xyz[i])
    return rho_phi_z


@boundscheck(False)
@wraparound(False)
cpdef double[:] cylindrical_to_cartesian_point(double[:] rho_phi_z):
    cdef:
        Py_ssize_t s = rho_phi_z.shape[0]
        array[double] xyz, template = array('d')
    if s < 3:
        rho_phi_z = __extend_vector_dimensions(rho_phi_z, 3)
    xyz = clone(template, 3, False)
    rho_phi_z[1] = __reduce_angle(rho_phi_z[1], center=False, positive=True)
    xyz[0] = rho_phi_z[0] * cos(rho_phi_z[1])
    xyz[1] = rho_phi_z[0] * sin(rho_phi_z[1])
    xyz[2] = rho_phi_z[2]
    return xyz


@boundscheck(False)
@wraparound(False)
cpdef double[:, :] cylindrical_to_cartesian(double[:, :] rho_phi_z):
    cdef:
        unsigned int i, s = rho_phi_z.shape[0]
        double[:, :] xyz = np.empty((s, 3), dtype=np.double)
    for i in range(s):
        xyz[i] = cylindrical_to_cartesian_point(rho_phi_z[i])
    return xyz


@boundscheck(False)
@wraparound(False)
cpdef double[:] spherical_to_cylindrical_point(double[:] r_theta_phi):
    cdef:
        Py_ssize_t s = r_theta_phi.shape[0]
        array[double] rho_phi_z, template = array('d')
    if s < 3:
        r_theta_phi = __extend_vector_dimensions(r_theta_phi, 3)
    rho_phi_z = clone(template, 3, False)
    rho_phi_z[0] = r_theta_phi[0] * sin(r_theta_phi[1])
    rho_phi_z[1] = r_theta_phi[2]
    rho_phi_z[2] = r_theta_phi[0] * cos(r_theta_phi[1])
    return rho_phi_z


@boundscheck(False)
@wraparound(False)
cpdef double[:, :] spherical_to_cylindrical(double[:, :] r_theta_phi):
    cdef:
        unsigned int i, s = r_theta_phi.shape[0]
        double[:, :] rho_phi_z = np.empty((s, 3), dtype=np.double)
    for i in range(s):
        rho_phi_z[i] = spherical_to_cylindrical_point(r_theta_phi[i])
    return rho_phi_z


@boundscheck(False)
@wraparound(False)
cpdef double[:] cylindrical_to_spherical_point(double[:] rho_phi_z):
    cdef:
        Py_ssize_t s = rho_phi_z.shape[0]
        array[double] r_theta_phi, template = array('d')
    if s < 3:
        rho_phi_z = __extend_vector_dimensions(rho_phi_z, 3)
    r_theta_phi = clone(template, 3, False)
    r_theta_phi[0] = sqrt(rho_phi_z[0] * rho_phi_z[0] + rho_phi_z[2] * rho_phi_z[2])
    r_theta_phi[1] = atan2(rho_phi_z[0], rho_phi_z[2])
    r_theta_phi[2] = rho_phi_z[1]
    return r_theta_phi


@boundscheck(False)
@wraparound(False)
cpdef double[:, :] cylindrical_to_spherical(double[:, :] rho_phi_z):
    cdef:
        unsigned int i, s = rho_phi_z.shape[0]
        r_theta_phi = np.empty((s, 3), dtype=np.double)
    for i in range(s):
        r_theta_phi[i] = cylindrical_to_spherical_point(rho_phi_z[i])
    return r_theta_phi
