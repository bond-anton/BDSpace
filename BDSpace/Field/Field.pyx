from __future__ import division, print_function
import numpy as np

from cython import boundscheck, wraparound

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

    def add_element(self, element):
        raise NotImplementedError('No element could be added to or removed from Field')

    def remove_element(self, element):
        raise NotImplementedError('No element could be added to or removed from Field')

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
