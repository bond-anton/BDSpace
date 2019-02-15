import numpy as np

from cython import boundscheck, wraparound

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
            Py_ssize_t i
            Py_ssize_t s = xyz.shape[0]
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        for i in range(s):
            values[i] = value
        return values

    @boundscheck(False)
    @wraparound(False)
    cdef double[:, :] __points_vector(self, double[:, :] xyz, double[:] value):
        cdef:
            Py_ssize_t i, j
            Py_ssize_t s = xyz.shape[0]
            Py_ssize_t c = xyz.shape[1]
            double[:, :] values = np.empty((s, c), dtype=np.double)
        for i in range(s):
            for j in range(c):
                values[i, j] = value[j]
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
        return self.__points_vector(xyz, clone(template, xyz.shape[1], zero=True))


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
            Py_ssize_t i, j
            Py_ssize_t s = xyz.shape[0]
            Py_ssize_t c = xyz.shape[1]
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        for i in range(s):
            values[i] = 0
            for j in range(c):
                values[i] += xyz[i, j] * self.__potential[j]
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

    cpdef double[:] scalar_field(self, double[:, :] xyz):
        cdef:
            Py_ssize_t i
            Py_ssize_t s = xyz.shape[0]
            Py_ssize_t c = xyz.shape[1]
            double r
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        for i in range(s):
            r = sqrt(xyz[i, 0]*xyz[i, 0] + xyz[i, 1]*xyz[i, 1] + xyz[i, 2]*xyz[i, 2])
            if r < self.__r:
                r = self.__r
            values[i] = self.__a / r
        return values

    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        cdef:
            Py_ssize_t i, j
            Py_ssize_t s = xyz.shape[0]
            Py_ssize_t c = xyz.shape[1]
            double r2, r2_min = self.__r * self.__r
            double[:, :] values = np.empty((s, c), dtype=np.double)
        for i in range(s):
            r2 = xyz[i, 0]*xyz[i, 0] + xyz[i, 1]*xyz[i, 1] + xyz[i, 2]*xyz[i, 2]
            if r2 < r2_min:
                r2 = r2_min
            for j in range(c):
                values[i, j] = self.__a * xyz[i, j] / r2
        return values
