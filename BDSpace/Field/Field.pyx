from __future__ import division, print_function
import numpy as np

from BDSpace.Space cimport Space


cdef class Field(Space):

    def __init__(self, str name, str field_type):
        self.__type = None
        self.type = field_type
        super(Field, self).__init__(name, coordinate_system=None)

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, field_type):
        self.__type = str(field_type)

    def __str__(self):
        description = 'Field: %s (%s)\n' % (self.name, self.type)
        if self.parent is not None:
            description += 'Parent entity:\n'
            description += str(self.parent)
        else:
            description += str(self.coordinate_system)
        return description

    def add_element(self, element):
        print('No element could be added to Field')
        pass

    def remove_element(self, element):
        print('No element could be added to or removed from Field')
        pass

    def scalar_field(self, xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of one or more points
        :return: scalar values array
        """
        xyz = np.array(xyz, dtype=np.float)
        if xyz.size == 3:
            field = 0
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                field = np.zeros_like(xyz.shape[0])
            else:
                raise ValueError('N-points array shape must be (N, 3)')
        else:
            raise ValueError('at least 3 coordinates are needed for point')
        return field

    def vector_field(self, xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of one or more points
        :return: vector field values array
        """
        xyz = np.array(xyz, dtype=np.float)
        if xyz.size == 3:
            field = np.array([0, 0, 0])
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                field = np.zeros_like(xyz)
            else:
                raise ValueError('N-points array shape must be (N, 3)')
        else:
            raise ValueError('at least 3 coordinates are needed for point')
        return field
