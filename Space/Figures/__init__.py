import numpy as np

from Space import Space
# from Space.Coordinates import Cartesian


class Figure(Space):

    def __init__(self, name, coordinate_system=None, points=None):
        super(Figure, self).__init__(name, coordinate_system=coordinate_system)
        self.points = None
        self.set_points(points)

    def set_points(self, points=None):
        if points is None:
            self.points = None
        else:
            self.points = np.array(points, dtype=np.float)
