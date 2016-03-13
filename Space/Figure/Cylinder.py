from __future__ import division
import numpy as np

from Space.Coordinates.transforms import reduce_angle
from Space.Figure import Figure


class CylindricalWedge(Figure):

    def __init__(self, name='Cylindrical wedge', coordinate_system=None,
                 r_inner=0.0, r_outer=1.0, phi=np.pi/2, z=np.array([0.0, 1.0])):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.phi = reduce_angle(float(phi))
        self.z = np.array(z, dtype=float)
        super(CylindricalWedge, self).__init__(name, coordinate_system=coordinate_system)

    def inner_volume(self):
        return 0.0

    def external_volume(self):
        h = abs(max(self.z) - min(self.z))
        v_inner = h * self.r_inner**2 * self.phi / 2
        v_outer = h * self.r_outer**2 * self.phi / 2
        return v_outer - v_inner

    def inner_surface_area(self):
        return 0.0

    def external_surface_area(self):
        h = abs(max(self.z) - min(self.z))
        s_outer = h * self.r_outer * self.phi
        s_inner = h * self.r_inner * self.phi
        s_bases = (self.r_outer**2 - self.r_inner**2) * self.phi / 2
        if self.phi == 2 * np.pi:
            s_cut = 0
        else:
            s_cut = 2 * h * (self.r_outer - self.r_inner)
        return s_outer + s_inner + s_bases + s_cut


class Cylinder(CylindricalWedge):
    def __init__(self, name='Cylinder', coordinate_system=None,
                 r_inner=0.0, r_outer=1.0, z=np.array([0.0, 1.0])):
        super(Cylinder, self).__init__(name, coordinate_system=coordinate_system,
                                       r_inner=r_inner, r_outer=r_outer, phi=np.pi*2, z=z)
