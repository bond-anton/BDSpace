from Space import Space


class Curve(Space):

    def __init__(self, name, coordinate_system=None):
        super(Curve, self).__init__(name, coordinate_system=coordinate_system)

    def __str__(self):
        return 'Curve: ' + self.name
