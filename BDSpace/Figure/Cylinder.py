import numpy as np

from BDSpace.Coordinates.transforms import reduce_angle
from BDSpace.Figure import Figure


class CylindricalWedge(Figure):

    def __init__(self, name='Cylindrical wedge', coordinate_system=None,
                 r_inner=0.0, r_outer=1.0, phi=np.pi/2, z=np.array([0.0, 1.0])):
        self.__r_inner = None
        self.r_inner = r_inner
        self.__r_outer = None
        self.r_outer = r_outer
        self.__phi = None
        self.phi = phi
        self.__z = None
        self.z = z
        super(CylindricalWedge, self).__init__(name, coordinate_system=coordinate_system)

    @property
    def r_inner(self):
        return self.__r_inner

    @r_inner.setter
    def r_inner(self, r_inner):
        self.__r_inner = np.float64(r_inner)

    @property
    def r_outer(self):
        return self.__r_outer

    @r_outer.setter
    def r_outer(self, r_outer):
        self.__r_outer = np.float64(r_outer)

    @property
    def phi(self):
        return self.__phi

    @phi.setter
    def phi(self, phi):
        self.__phi = reduce_angle(np.float64(phi))

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, z):
        self.__z = np.array(z, dtype=np.float64)

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
