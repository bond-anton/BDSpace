from __future__ import division
import numpy as np

from Space.Figure import Figure


class CylindricalWedge(Figure):

    def __init__(self, name='Cylindrical wedge', coordinate_system=None,
                 r_inner=0.0, r_outer=1.0, phi=np.pi/2, z=[0.0, 1.0]):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.phi = phi
        self.z = np.array(z, dtype=float)
        super(CylindricalWedge, self).__init__(name, coordinate_system=coordinate_system)


class Cylinder(CylindricalWedge):
    def __init__(self, name='Cylinder', coordinate_system=None,
                 r_inner=0.0, r_outer=1.0, z=[0.0, 1.0]):
        super(Cylinder, self).__init__(name, coordinate_system=coordinate_system,
                                       r_inner=r_inner, r_outer=r_outer, phi=np.pi*2, z=z)
