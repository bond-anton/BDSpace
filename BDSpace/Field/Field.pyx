import numpy as np

from cython import boundscheck, wraparound
from cython.parallel import prange

from libc.math cimport sqrt
from cpython.array cimport array, clone

from BDSpace.Space cimport Space


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

    cpdef double[:] scalar_field(self, double[:, :] xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
        return self.__points_scalar(xyz, 0.0)

    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            array[double] template = array('d')
        return self.__points_vector(xyz, clone(template, 3, zero=True))


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

    cpdef double[:] scalar_field(self, double[:, :] xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
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
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                values[i] = xyz[i, 0] * self.__potential[0] + xyz[i, 1] * self.__potential[1]\
                            + xyz[i, 2] * self.__potential[2]
        return values

    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: vector field values array
        """
        return self.__points_vector(xyz, self.__potential)


cdef class HyperbolicPotentialSphericalConservativeField(Field):

    def __init__(self, str name, str field_type, double a, double r):
        self.__r = r
        self.__a = a
        super(HyperbolicPotentialSphericalConservativeField, self).__init__(name, field_type)

    @property
    def r(self):
        return self.__r

    @r.setter
    def r(self, double r):
        self.__r = r

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, double a):
        self.__a = a

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field(self, double[:, :] xyz):
        cdef:
            int i, s = xyz.shape[0]
            double r
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                r = sqrt(xyz[i, 0] * xyz[i, 0] + xyz[i, 1] * xyz[i, 1] + xyz[i, 2] * xyz[i, 2])
                if r < self.__r:
                    r = self.__r
                values[i] = self.__a / r
        return values

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        cdef:
            int i, j, s = xyz.shape[0]
            double r, r2, r2_min = self.__r * self.__r
            double[:, :] values = np.empty((s, 3), dtype=np.double)
        with nogil:
            for i in prange(s):
                r2 = xyz[i, 0] * xyz[i, 0] + xyz[i, 1] * xyz[i, 1] + xyz[i, 2] * xyz[i, 2]
                if r2 >= r2_min:
                    r = sqrt(r2)
                    values[i, 0] = self.__a * xyz[i, 0] / r2 / r
                    values[i, 1] = self.__a * xyz[i, 1] / r2 / r
                    values[i, 2] = self.__a * xyz[i, 2] / r2 / r
                else:
                    values[i, 0] = 0.0
                    values[i, 0] = 0.0
                    values[i, 0] = 0.0
        return values
