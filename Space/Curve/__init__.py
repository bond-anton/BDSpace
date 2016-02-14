from Space.Coordinates.coordinates import Cartesian

class Curve(object):

    def __init__(self, name='Unknown curve', CS=None):
        self.name = str(name)
        if CS is None:
            self.CS = Cartesian()
        else:
            self.CS = CS

    def __str__(self):
        return 'Curve: ' + self.name
