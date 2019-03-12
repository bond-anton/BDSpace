import numpy as np

from cython import boundscheck, wraparound
from cython.parallel import prange

from libc.math cimport fabs, fmod, sqrt, M_PI
from cpython.array cimport array, clone

from .Field cimport Field
from BDSpace.Coordinates.transforms cimport cartesian_to_spherical_point, cartesian_to_spherical
from BDSpace.Coordinates.transforms cimport spherical_to_cartesian_point, spherical_to_cartesian


cdef class SphericallySymmetric(Field):

    def __init__(self, str name, str field_type, double r):
        self.__r = r
        super(SphericallySymmetric, self).__init__(name, field_type)

    @property
    def r(self):
        return self.__r

    @r.setter
    def r(self, double r):
        self.__r = r

    @boundscheck(False)
    @wraparound(False)
    cdef inline double __get_r(self, double[:] xyz) nogil:
        return sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2])

    cdef double scalar_field_r_law(self, double r) nogil:
        return 0.0

    cpdef double scalar_field_r_point(self, double r):
        return self.scalar_field_r_law(r)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field_r(self, double[:] r):
        cdef:
            int i, s = r.shape[0]
            array[double] result, template = array('d')
        result = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                result[i] = self.scalar_field_r_law(r[i])
        return result

    cdef double vector_field_r_law(self, double r) nogil:
        return 0.0

    cpdef double vector_field_r_point(self, double r):
        return self.vector_field_r_law(r)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] vector_field_r(self, double[:] r):
        cdef:
            int i, s = r.shape[0]
            array[double] result, template = array('d')
        result = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                result[i] = self.vector_field_r_law(r[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field(self, double[:, :] xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
        cdef:
            int i, s = xyz.shape[0]
            array[double] result, template = array('d')
        result = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                result[i] = self.scalar_field_r_law(self.__get_r(xyz[i]))
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double scalar_field_polar_point(self, double[:] rtp):
        """
        Calculates scalar field value at point rtp in polar coordinates
        :param rtp: array of polar coordinates of the point
        :return: scalar field value
        """
        return self.scalar_field_r_law(rtp[0])

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field_polar(self, double[:, :] rtp):
        """
        Calculates scalar field value at points rtp in polar coordinates
        :param rtp: array of N points with shape (N, 3)
        :return: scalar values array
        """
        cdef:
            int i, s = rtp.shape[0]
            array[double] result, template = array('d')
        result = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                result[i] = self.scalar_field_r_law(rtp[i, 0])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] vector_field_polar_point(self, double[:] rtp):
        """
        Calculates vector field value at point rtp in polar coordinates
        :param rtp: array of polar coordinates of the point
        :return: vector field value
        """
        cdef:
            array[double] result = clone(array('d'), 3, zero=False)
            double mag = self.vector_field_r_law(rtp[0])
        if mag < 0:
            result[0] = fabs(mag)
            result[1] = M_PI - rtp[1]
            result[2] = fmod((rtp[2] + M_PI), (2 * M_PI))
        else:
            result[0] = mag
            result[1] = rtp[1]
            result[2] = rtp[2]
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field_polar(self, double[:, :] rtp):
        """
        Calculates vector field value at points rtp in polar coordinates
        :param rtp: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            int i, s = rtp.shape[0]
            double[:, :] result = np.empty((s, 3), dtype=np.double)
            double mag
        with nogil:
            for i in prange(s):
                mag = self.vector_field_r_law(rtp[i, 0])
                if mag < 0:
                    result[i, 0] = fabs(mag)
                    result[i, 1] = M_PI - rtp[i, 1]
                    result[i, 2] = fmod((rtp[i, 2] + M_PI), (2 * M_PI))
                else:
                    result[i, 0] = mag
                    result[i, 1] = rtp[i, 1]
                    result[i, 2] = rtp[i, 2]
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] vector_field_point(self, double[:] xyz):
        """
        Calculates vector field value at point xyz
        :param xyz: array of cartesian coordinates of the point
        :return: vector field value
        """
        return spherical_to_cartesian_point(self.vector_field_polar_point(cartesian_to_spherical_point(xyz)))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: vector field values array
        """
        return spherical_to_cartesian(self.vector_field_polar(cartesian_to_spherical(xyz)))


cdef class HyperbolicPotentialSphericalConservativeField(SphericallySymmetric):

    def __init__(self, str name, str field_type, double r, double a):
        self.__r = r
        self.__a = a
        super(HyperbolicPotentialSphericalConservativeField, self).__init__(name, field_type, r)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, double a):
        self.__a = a

    cdef double scalar_field_r_law(self, double r) nogil:
        if r < self.__r:
            return self.__a / self.__r
        return self.__a / r

    cdef double vector_field_r_law(self, double r) nogil:
        if r < self.__r:
            return 0.0
        return self.__a / (r * r)
