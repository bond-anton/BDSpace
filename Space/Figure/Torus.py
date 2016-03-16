from __future__ import division
import numpy as np

from Space.Coordinates.transforms import reduce_angle
from Space.Figure import Figure


class ToricWedge(Figure):

    def __init__(self, name='Toric wedge', coordinate_system=None,
                 phi=np.pi/2, theta=np.array([0.0, np.pi/2]),
                 r_torus=1.0, r_tube=np.array([0, 0.25])):
        self.r_torus = float(r_torus)
        self.r_tube = np.array(r_tube, dtype=np.float)
        self.r_tube = self.r_tube.reshape(self.r_tube.size,)
        reduced_theta = reduce_angle(theta, keep_sign=True)
        theta_min = min(reduced_theta)
        theta_max = max(reduced_theta)
        if theta_max - theta_min >= 2*np.pi:
            theta_min = 0.0
            theta_max = 2 * np.pi
        self.theta = np.array([theta_min, theta_max])
        self.phi = reduce_angle(phi)

        super(ToricWedge, self).__init__(name, coordinate_system=coordinate_system)

    def inner_volume(self):
        if np.allclose(self.theta[1] - self.theta[0], 2*np.pi) and np.allclose(self.phi, 2*np.pi):
            return 2 * np.pi**2 * self.r_torus * self.r_tube**2
        else:
            return 0



class ToricSector(ToricWedge):

    def __init__(self, name='Toric sector', coordinate_system=None,
                 phi=np.pi/2,
                 r_torus=1.0, r_tube=np.array([0, 0.25])):

        super(ToricSector, self).__init__(name, coordinate_system=coordinate_system,
                                          phi=phi, theta=np.array([0, 2*np.pi]),
                                          r_torus=r_torus, r_tube=r_tube)


class Torus(ToricSector):

    def __init__(self, name='Torus', coordinate_system=None,
                 r_torus=1.0, r_tube=np.array([0, 0.25])):

        super(ToricSector, self).__init__(name, coordinate_system=coordinate_system,
                                          phi=2*np.pi, r_torus=r_torus, r_tube=r_tube)
