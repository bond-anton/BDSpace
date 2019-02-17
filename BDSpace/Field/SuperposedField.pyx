import numpy as np

from cython import boundscheck, wraparound
from cython.parallel import prange

from cpython.array cimport array, clone

from .Field cimport Field


cdef class SuperposedField(Field):

    def __init__(self, str name, list fields):
        self.__fields = []
        self.type = None
        self.fields = fields
        super(SuperposedField, self).__init__(name, self.type)

    @property
    def fields(self):
        return self.__fields

    @fields.setter
    def fields(self, list fields):
        self.__fields = []
        for field in fields:
            if not isinstance(field, Field):
                raise ValueError('Fields must be iterable of Field class instances')
            if self.type is None:
                self.type = field.type
            if self.type != field.type:
                raise ValueError('All fields must be iterable of Field class instances')
            self.__fields.append(field)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field(self, double[:, :] xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
        cdef:
            int i, j, s = xyz.shape[0], n_fields = len(self.__fields)
            double[:, :] global_xyz = self.to_global_coordinate_system(xyz)
            double[:] total_field = self.__points_scalar(global_xyz, 0.0)
            double[:] field_contribution
        for j in range(n_fields):
            field_contribution = self.__fields[j].scalar_field(self.__fields[j].to_local_coordinate_system(global_xyz))
            with nogil:
                for i in prange(s):
                    total_field[i] += field_contribution[i]
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            int i, j, k, s = xyz.shape[0], n_fields = len(self.__fields)
            array[double] template = array('d')
            double[:, :] field_contribution
            double[:, :] total_field = self.__points_vector(xyz, clone(template, xyz.shape[1], zero=True))
            double[:, :] global_xyz = self.to_global_coordinate_system(xyz)
        for k in range(n_fields):
            field_contribution = self.__fields[k].vector_field(self.__fields[k].to_local_coordinate_system(global_xyz))
            with nogil:
                for i in prange(s):
                    for j in prange(3):
                        total_field[i][j] += field_contribution[i][j]
        return total_field
