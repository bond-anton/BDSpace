from Coordinates import Cartesian


class Space(object):

    def __init__(self, name, coordinates_system=None):
        if coordinates_system is None:
            self.coordinates_system = Cartesian()
        elif isinstance(coordinates_system, Cartesian):
            self.coordinates_system = coordinates_system
        else:
            raise ValueError('coordinates system must be instance of Cartesian class')
        self.name = str(name)
        self.parent = None
        self.elements = {}

    def __str__(self):
        return 'Space: ' + self.name

    def add_element(self, element):
        if isinstance(element, Space):
            if element == self:
                raise ValueError('Space can not be its own subspace')
            else:
                if element.parent is None:
                    element_name = element.name
                    name_counter = 1
                    format = ' %d'
                    while element_name in self.elements.keys():
                        element_name = element.name + format % name_counter
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
