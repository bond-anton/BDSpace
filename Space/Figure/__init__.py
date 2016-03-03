from Space import Space


class Figure(Space):

    def __init__(self, name, coordinate_system=None):
        super(Figure, self).__init__(name, coordinate_system=coordinate_system)
