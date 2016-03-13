from __future__ import division
import numpy as np

from Space.Coordinates.transforms import reduce_angle
from Space.Figure import Figure


class ConicalWedge(Figure):

    def __init__(self, name='Conical wedge', coordinate_system=None,
                 phi=np.pi/2, theta=np.pi/6, z=np.array([0.0, 1.0]), z_offset=0.0, r_min=0.0):
        reduced_angle = reduce_angle(float(theta))
        if reduced_angle > np.pi/2:
            raise ValueError('Cone half angle should be between 0 and 2*pi radians')
        self.theta = reduced_angle
        self.phi = reduce_angle(float(phi))
        self.r_min = abs(float(r_min))
        z_min = r_min / np.tan(theta)
        self.z_offset = abs(float(z_offset))
        z_points = np.array(z, dtype=float)
        if np.sign(min(z_points) * max(z_points)) <= 0:
            z_points = np.union1d(z_points, np.array([0.0]))
        if z_offset - z_min > 0:
            if max(z_points) > self.z_offset - z_min > min(z_points):
                z_points = np.union1d(z_points, np.array([self.z_offset - z_min]))
            if min(z_points) < -self.z_offset + z_min < max(z_points):
                z_points = np.union1d(z_points, np.array([-self.z_offset + z_min]))
        self.z = z_points

        super(ConicalWedge, self).__init__(name, coordinate_system=coordinate_system)
