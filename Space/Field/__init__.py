import numpy as np
from Space import Space


class Field(Space):

    def __init__(self, name, parent, field_type):
        if not isinstance(parent, Space):
            raise ValueError('parent must be an instance of Space class')
        self.type = str(field_type)
        super(Field, self).__init__(name, coordinate_system=None)
        parent.add_element(self)

    def __str__(self):
        return 'Field: %s (%s)' % (self.name, self.type)

    def add_element(self, element):
        print 'No element could be added to Field'
        pass

    def remove_element(self, element):
        print 'No element could be added to or removed from Field'
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
                field = np.zeros_like(len(xyz))
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
