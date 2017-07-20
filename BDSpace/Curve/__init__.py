from BDSpace import Space


class Curve(Space):

    def __init__(self, name, coordinate_system=None):
        super(Curve, self).__init__(name, coordinate_system=coordinate_system)

    def __str__(self):
        description = 'Curve: %s\n' % self.name
        description += str(self.coordinate_system)
        return description
