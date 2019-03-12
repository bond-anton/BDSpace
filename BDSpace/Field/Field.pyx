import numpy as np

from cython import boundscheck, wraparound
from cython.parallel import prange

from libc.math cimport sqrt
from cpython.array cimport array, clone

from BDSpace.Space cimport Space
from BDSpace.Coordinates.transforms cimport spherical_to_cartesian_point, spherical_to_cartesian


cdef class Field(Space):

    def __init__(self, str name, str field_type):
        self.__type = field_type
        super(Field, self).__init__(name, coordinate_system=None)

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, str field_type):
        self.__type = field_type

    def __str__(self):
        description = 'Field: %s (%s)\n' % (self.name, self.type)
        if self.parent is not None:
            description += 'Parent entity:\n'
            description += str(self.parent)
        else:
            description += str(self.coordinate_system)
        return description

    cpdef bint add_element(self, Space element):
        return False

    cpdef bint remove_element(self, Space element):
        return False

    @boundscheck(False)
    @wraparound(False)
    cdef double[:] __points_scalar(self, double[:, :] xyz, double value):
        cdef:
            int i, s = xyz.shape[0]
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                values[i] = value
        return values

    @boundscheck(False)
    @wraparound(False)
    cdef double[:, :] __points_vector(self, double[:, :] xyz, double[:] value):
        cdef:
            int i, s = xyz.shape[0]
            double[:, :] values = np.empty((s, 3), dtype=np.double)
        with nogil:
            for i in prange(s):
                values[i, 0] = value[0]
                values[i, 1] = value[1]
                values[i, 2] = value[2]
        return values

    cpdef double scalar_field_point(self, double[:] xyz):
        """
        Calculates scalar field value at point xyz
        :param xyz: array of cartesian coordinates of the point
        :return: scalar field value
        """
        return 0.0

    cpdef double[:] scalar_field(self, double[:, :] xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
        return self.__points_scalar(xyz, 0.0)

    cpdef double scalar_field_polar_point(self, double[:] rtp):
        """
        Calculates scalar field value at point rtp in polar coordinates
        :param rtp: array of polar coordinates of the point
        :return: scalar field value
        """
        cdef:
            double[:] xyz = spherical_to_cartesian_point(rtp)
        return self.scalar_field_point(xyz)

    cpdef double[:] scalar_field_polar(self, double[:, :] rtp):
        """
        Calculates scalar field value at points rtp in polar coordinates
        :param rtp: array of N points with shape (N, 3)
        :return: scalar values array
        """
        cdef:
            double[:, :] xyz = spherical_to_cartesian(rtp)
        return self.scalar_field(xyz)

    cpdef double[:] vector_field_point(self, double[:] xyz):
        """
        Calculates vector field value at point xyz
        :param xyz: array of cartesian coordinates of the point
        :return: vector field value
        """
        cdef:
            array[double] template = array('d')
        return clone(template, 3, zero=True)

    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            array[double] template = array('d')
        return self.__points_vector(xyz, clone(template, 3, zero=True))

    cpdef double[:] vector_field_polar_point(self, double[:] rtp):
        """
        Calculates vector field value at point rtp in polar coordinates
        :param rtp: array of polar coordinates of the point
        :return: vector field value
        """
        cdef:
            double[:] xyz = spherical_to_cartesian_point(rtp)
        return self.vector_field_point(xyz)

    cpdef double[:, :] vector_field_polar(self, double[:, :] rtp):
        """
        Calculates vector field value at points rtp in polar coordinates
        :param rtp: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            double[:, :] xyz = spherical_to_cartesian(rtp)
        return self.vector_field(xyz)


cdef class ConstantScalarConservativeField(Field):

    def __init__(self, str name, str field_type, double potential):
        self.__potential = potential
        super(ConstantScalarConservativeField, self).__init__(name, field_type)

    @property
    def potential(self):
        return self.__potential

    @potential.setter
    def potential(self, double potential):
        self.__potential = potential

    cpdef double scalar_field_point(self, double[:] xyz):
        return self.__potential

    cpdef double[:] scalar_field(self, double[:, :] xyz):
        return self.__points_scalar(xyz, self.__potential)


cdef class ConstantVectorConservativeField(Field):

    def __init__(self, str name, str field_type, double[:] potential):
        self.__potential = potential
        super(ConstantVectorConservativeField, self).__init__(name, field_type)

    @property
    def potential(self):
        return self.__potential

    @potential.setter
    def potential(self, double[:] potential):
        self.__potential = potential

    cpdef double scalar_field_point(self, double[:] xyz):
        return xyz[0] * self.__potential[0] + xyz[1] * self.__potential[1] + xyz[2] * self.__potential[2]

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field(self, double[:, :] xyz):
        cdef:
            int i, s = xyz.shape[0]
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                values[i] = xyz[i, 0] * self.__potential[0] + xyz[i, 1] * self.__potential[1]\
                            + xyz[i, 2] * self.__potential[2]
        return values

    cpdef double[:] vector_field_point(self, double[:] xyz):
        return self.__potential

    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        return self.__points_vector(xyz, self.__potential)
