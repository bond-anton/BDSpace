import numpy as np

from cython import boundscheck, wraparound

from cpython.array cimport array, clone

from BDSpace.Coordinates.Cartesian cimport Cartesian

from ._version import __version__


cdef class Space(object):

    def __init__(self, str name, Cartesian coordinate_system=None):
        self.__name = name
        if coordinate_system is None:
            self.__coordinate_system = Cartesian()
        else:
            self.__coordinate_system = coordinate_system
        self.__parent = None
        self.__elements = {}

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, str name):
        self.__name = name

    @property
    def coordinate_system(self):
        return self.__coordinate_system

    @coordinate_system.setter
    def coordinate_system(self, Cartesian coordinate_system):
        self.__coordinate_system = coordinate_system

    @property
    def parent(self):
        return self.__parent

    @parent.setter
    def parent(self, parent):
        if isinstance(parent, Space) or parent is None:
            self.__parent = parent
        else:
            raise ValueError('Only Space object or None are accepted for parent.')

    @property
    def elements(self):
        return self.__elements

    def __str__(self):
        description = 'BDSpace: %s\n' % self.name
        description += str(self.coordinate_system)
        return description

    cpdef double[:] to_global_coordinate_system_vector(self, double[:] xyz):
        """
        convert local points coordinates xyz to global coordinate system coordinates
        :param xyz: 3D vector
        :return: 3D vector in global coordinates system
        """
        cdef:
            double[:] parent_xyz = self.coordinate_system.to_parent_vector(xyz)
        if self.__parent is None:
            return parent_xyz
        else:
            return self.parent.to_global_coordinate_system_vector(parent_xyz)

    cpdef double[:, :] to_global_coordinate_system(self, double[:, :] xyz):
        """
        convert local points coordinates xyz to global coordinate system coordinates
        :param xyz: array of points shaped Nx3
        :return: array of points in global coordinates system
        """
        parent_xyz = self.coordinate_system.to_parent(xyz)
        if self.parent is None:
            return parent_xyz
        else:
            return self.parent.to_global_coordinate_system(parent_xyz)

    @boundscheck(False)
    cpdef Cartesian basis_in_global_coordinate_system(self):
        """
        returns local coordinate system basis in global coordinate system as Cartesian class object
        :return: local Cartesian coordinate system in global coordinate system
        """
        cdef:
            double[:, :] basis = np.empty((3, 3), dtype=np.double)
            double[:] origin
            array[double] template = array('d')
        origin = clone(template, 3, zero=False)
        origin[0] = self.coordinate_system.origin[0]
        origin[1] = self.coordinate_system.origin[1]
        origin[2] = self.coordinate_system.origin[2]
        basis[0, 0] = self.coordinate_system.basis[0, 0]
        basis[0, 1] = self.coordinate_system.basis[0, 1]
        basis[0, 2] = self.coordinate_system.basis[0, 2]
        basis[1, 0] = self.coordinate_system.basis[1, 0]
        basis[1, 1] = self.coordinate_system.basis[1, 1]
        basis[1, 2] = self.coordinate_system.basis[1, 2]
        basis[2, 0] = self.coordinate_system.basis[2, 0]
        basis[2, 1] = self.coordinate_system.basis[2, 1]
        basis[2, 2] = self.coordinate_system.basis[2, 2]
        if self.parent is not None:
            basis[0, 0] += origin[0]
            basis[0, 1] += origin[1]
            basis[0, 2] += origin[2]
            basis[1, 0] += origin[0]
            basis[1, 1] += origin[1]
            basis[1, 2] += origin[2]
            basis[2, 0] += origin[0]
            basis[2, 1] += origin[1]
            basis[2, 2] += origin[2]
            basis = self.parent.to_global_coordinate_system(basis)
            origin = self.parent.to_global_coordinate_system_vector(origin)
            basis[0, 0] -= origin[0]
            basis[0, 1] -= origin[1]
            basis[0, 2] -= origin[2]
            basis[1, 0] -= origin[0]
            basis[1, 1] -= origin[1]
            basis[1, 2] -= origin[2]
            basis[2, 0] -= origin[0]
            basis[2, 1] -= origin[1]
            basis[2, 2] -= origin[2]
        name = self.coordinate_system.name
        labels = self.coordinate_system.labels
        coordinate_system = Cartesian(basis=basis, origin=origin, name=name, labels=labels)
        return coordinate_system

    cpdef double[:] to_local_coordinate_system_vector(self, double[:] xyz):
        """
        convert global points coordinates xyz to local coordinate system coordinates
        :param xyz: array of points shaped Nx3
        :return: array of points in local coordinates system
        """
        cdef:
            Cartesian basis_global = self.basis_in_global_coordinate_system()
        return basis_global.to_local_vector(xyz)

    cpdef double[:, :] to_local_coordinate_system(self, double[:, :] xyz):
        """
        convert global points coordinates xyz to local coordinate system coordinates
        :param xyz: array of points shaped Nx3
        :return: array of points in local coordinates system
        """
        cdef:
            Cartesian basis_global = self.basis_in_global_coordinate_system()
        return basis_global.to_local(xyz)

    cpdef bint add_element(self, Space element):
        if element == self:
            return False
        else:
            if element.parent is None:
                element_name = element.name
                name_counter = 1
                name_format = ' %d'
                while element_name in self.elements.keys():
                    element_name = element.name + name_format % name_counter
                    name_counter += 1
                element.name = element_name
                self.elements[element.name] = element
                element.parent = self
                return True
            elif element.parent == self:
                return False
            else:
                return False

    @boundscheck(False)
    @wraparound(False)
    cpdef bint remove_element(self, Space element):
        if element.parent == self:
            for key in self.elements.keys():
                if self.elements[key] == element:
                    del self.elements[key]
            element.parent = None
            return True
        return False

    cpdef void detach_from_parent(self):
        if self.__parent is not None:
            self.__parent.remove_element(self)

    @boundscheck(False)
    @wraparound(False)
    cpdef void print_tree(self, int level=0):
        print('-' * level + ' ' * (level > 0) + self.name)
        for key in self.elements.keys():
            self.elements[key].print_tree(level=level+1)
