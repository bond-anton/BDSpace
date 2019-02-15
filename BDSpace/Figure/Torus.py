import numpy as np

from BDSpace.Coordinates.transforms import reduce_angle
from BDSpace.Figure import Figure


class ToricWedge(Figure):

    def __init__(self, name='Toric wedge', coordinate_system=None,
                 phi=np.pi/2, theta=np.array([0.0, np.pi/2]),
                 r_torus=1.0, r_tube=np.array([0, 0.25])):
        self.__r_torus = None
        self.r_torus = r_torus
        self.__r_tube = None
        self.r_tube = r_tube
        self.__theta = None
        self.theta = theta
        self.__phi = None
        self.phi = reduce_angle(phi)
        super(ToricWedge, self).__init__(name, coordinate_system=coordinate_system)

    @property
    def r_torus(self):
        return self.__r_torus

    @r_torus.setter
    def r_torus(self, r_torus):
        self.r_torus = abs(np.float64(r_torus))

    @property
    def r_tube(self):
        return self.__r_tube

    @r_tube.setter
    def r_tube(self, r_tube):
        self.__r_tube = abs(np.array(r_tube, dtype=np.float64))
        self.__r_tube = self.__r_tube.reshape(self.r_tube.size, )
        # Torus is not allowed to be thicker then its radius
        self.__r_tube[np.where(self.__r_tube > self.r_torus)] = self.r_torus

    @property
    def theta(self):
        return self.__theta

    @theta.setter
    def theta(self, theta):
        reduced_theta = reduce_angle(theta, keep_sign=True)
        theta_min = min(reduced_theta)
        theta_max = max(reduced_theta)
        if theta_max - theta_min >= 2 * np.pi:
            theta_min = 0.0
            theta_max = 2 * np.pi
        self.__theta = np.array([theta_min, theta_max])

    @property
    def phi(self):
        return self.__phi

    @phi.setter
    def phi(self, phi):
        self.__phi = reduce_angle(np.float64(phi))

    def inner_volume(self):
        if np.allclose(self.theta[1] - self.theta[0], 2*np.pi) and np.allclose(self.phi, 2*np.pi):
            return 2 * np.pi**2 * self.r_torus * min(self.r_tube)**2
        else:
            return 0

    def external_volume(self):
        external_volume = self.r_torus * max(self.r_tube)**2 * self.phi * (self.theta[1] - self.theta[0]) / 2
        if np.allclose(self.theta[1] - self.theta[0], 2*np.pi) and np.allclose(self.phi, 2*np.pi):
            return external_volume
        else:
            internal_volume = self.r_torus * min(self.r_tube)**2 * self.phi * (self.theta[1] - self.theta[0]) / 2
            return external_volume - internal_volume

    def inner_surface_area(self):
        if np.allclose(self.theta[1] - self.theta[0], 2*np.pi) and np.allclose(self.phi, 2*np.pi):
            return 4 * np.pi**2 * self.r_torus * min(self.r_tube)
        else:
            return 0

    def external_surface_area(self):
        if np.allclose(self.theta[1] - self.theta[0], 2*np.pi) and np.allclose(self.phi, 2*np.pi):
            return 4 * np.pi**2 * self.r_torus * max(self.r_tube)
        else:
            external_toric = self.r_torus * max(self.r_tube) * self.phi * (self.theta[1] - self.theta[0]) / 2
            internal_toric = self.r_torus * min(self.r_tube) * self.phi * (self.theta[1] - self.theta[0]) / 2
            if np.allclose(self.theta[1] - self.theta[0], 2*np.pi):
                tangent_sides = 0
            else:
                tangent_sides = (self.r_torus + max(self.r_tube))**2 * self.phi / 2
                tangent_sides -= (self.r_torus - max(self.r_tube))**2 * self.phi / 2
                tangent_sides -= (self.r_torus + min(self.r_tube))**2 * self.phi / 2
                tangent_sides += (self.r_torus - min(self.r_tube))**2 * self.phi / 2
            if np.allclose(self.phi, 2*np.pi):
                normal_sides = 0
            else:
                normal_sides = (max(self.r_tube) - min(self.r_tube))**2 * (self.theta[1] - self.theta[0])
            return external_toric + internal_toric + tangent_sides + normal_sides


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
