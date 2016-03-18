import numpy as np
from Coordinates import Cartesian
from _version import __version__


class Space(object):

    def __init__(self, name, coordinate_system=None):
        if coordinate_system is None:
            self.coordinate_system = Cartesian()
        elif isinstance(coordinate_system, Cartesian):
            self.coordinate_system = coordinate_system
        else:
            raise ValueError('coordinates system must be instance of Cartesian class')
        self.name = str(name)
        self.parent = None
        self.elements = {}
        self.points = None

    def __str__(self):
        return 'Space: ' + self.name

    def to_global_coordinate_system(self, xyz):
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

    def basis_in_global_coordinate_system(self):
        """
        returns local coordinate system basis in global coordinate system as Cartesian class object
        :return: local Cartesian coordinate system in global coordinate system
        """
        origin = np.copy(self.coordinate_system.origin)
        basis = np.copy(self.coordinate_system.basis)
        if self.parent is not None:
            basis = self.parent.to_global_coordinate_system(basis + origin)
            origin = self.parent.to_global_coordinate_system(origin)
            basis = basis - origin
        name = self.coordinate_system.name
        labels = self.coordinate_system.labels
        coordinate_system = Cartesian(basis=basis, origin=origin, name=name, labels=labels)
        return coordinate_system

    def add_element(self, element):
        if isinstance(element, Space):
            if element == self:
                raise ValueError('Space can not be its own subspace')
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
                elif element.parent == self:
                    print 'Space ' + element.name + 'is already a subspace of ' + self.name
                else:
                    raise ValueError(element.name + ' is subspace of another space: ' + element.parent.name +
                                     '. Please delete first.')
        else:
            raise ValueError('Only another space could be included as a subspace')

    def remove_element(self, element):
        if isinstance(element, Space):
            if element.parent == self:
                for key in self.elements.keys():
                    if self.elements[key] == element:
                        del self.elements[key]
                element.parent = None
            else:
                print element.name + ' is not a subspace of ' + self.name
        else:
            raise ValueError('Only another space could be detached')

    def detach_from_parent(self):
        if self.parent is not None:
            self.parent.remove_element(self)

    def print_tree(self, level=0):
        print '-' * level + ' ' * (level > 0) + self.name
        for key in self.elements.keys():
            self.elements[key].print_tree(level=level+1)
