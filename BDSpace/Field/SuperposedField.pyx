from __future__ import division, print_function
import numpy as np

from .Field cimport Field


cdef class SuperposedField(Field):

    def __init__(self, name, fields):
        self.__fields = None
        self.type = None
        self.fields = fields
        super(SuperposedField, self).__init__(name, self.type)

    @property
    def fields(self):
        return self.__fields

    @fields.setter
    def fields(self, fields):
        for field in fields:
            if not isinstance(field, Field):
                raise ValueError('Fields must be iterable of Field class instances')
            if self.type is None:
                self.type = field.type
            if self.type != field.type:
                raise ValueError('All fields must be iterable of Field class instances')

        self.__fields = fields

    def scalar_field(self, xyz):
        """
        Calculates superposed scalar field value at points xyz
        :param xyz: array of one or more points in global coordinate system
        :return: scalar values array
        """
        xyz = np.array(xyz, dtype=np.float)
        if xyz.size == 3:
            total_field = 0
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                total_field = np.zeros((xyz.shape[0],))
            else:
                raise ValueError('N-points array shape must be (N, 3)')
        else:
            raise ValueError('at least 3 coordinates are needed for point')
        for field in self.fields:
            total_field += field.scalar_field(field.to_local_coordinate_system(xyz))
        return total_field

    def vector_field(self, xyz):
        """
        Calculates superposed vector field value at points xyz
        :param xyz: array of one or more points in global coordinate system
        :return: vector values array
        """
        xyz = np.array(xyz, dtype=np.float)
        if xyz.size == 3:
            total_field = np.array([0, 0, 0])
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                total_field = np.zeros_like(xyz)
            else:
                raise ValueError('N-points array shape must be (N, 3)')
        else:
            raise ValueError('at least 3 coordinates are needed for point')
        for field in self.fields:
            vector_field_local = field.vector_field(field.to_local_coordinate_system(xyz))
            total_field += field.to_global_coordinate_system(vector_field_local)
        return total_field
