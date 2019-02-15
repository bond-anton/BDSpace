import numpy as np

from BDSpace.Coordinates.transforms import reduce_angle
from BDSpace.Figure import Figure


class ConicalWedge(Figure):

    def __init__(self, name='Conical wedge', coordinate_system=None,
                 phi=np.pi/2, theta=np.pi/6, z=np.array([0.0, 1.0]), z_offset=0.0, r_min=0.0):
        self.__theta = None
        self.theta = theta
        self.__phi = None
        self.phi = phi
        self.__r_min = None
        self.r_min = r_min
        self.__z_offset = None
        self.z_offset = z_offset
        self.__z = None
        self.z = z
        super(ConicalWedge, self).__init__(name, coordinate_system=coordinate_system)

    @property
    def theta(self):
        return self.__theta

    @theta.setter
    def theta(self, theta):
        reduced_angle = reduce_angle(np.float64(theta))
        if reduced_angle > np.pi / 2:
            raise ValueError('Cone half angle should be between 0 and 2*pi radians')
        self.__theta = reduced_angle

    @property
    def phi(self):
        return self.__phi

    @phi.setter
    def phi(self, phi):
        self.__phi = reduce_angle(np.float64(phi))

    @property
    def r_min(self):
        return self.__r_min

    @r_min.setter
    def r_min(self, r_min):
        self.__r_min = abs(np.float64(r_min))

    @property
    def z_offset(self):
        return self.__z_offset

    @z_offset.setter
    def z_offset(self, z_offset):
        self.__z_offset = abs(np.float64(z_offset))

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, z):
        z_min = self.r_min / np.tan(self.theta)
        z_points = np.array(z, dtype=np.float64)
        if np.sign(min(z_points) * max(z_points)) <= 0:
            z_points = np.union1d(z_points, np.array([0.0]))
        if self.z_offset - z_min > 0:
            if max(z_points) > self.z_offset - z_min > min(z_points):
                z_points = np.union1d(z_points, np.array([self.z_offset - z_min]))
            if min(z_points) < -self.z_offset + z_min < max(z_points):
                z_points = np.union1d(z_points, np.array([-self.z_offset + z_min]))
        self.__z = z_points

    def inner_volume(self):
        return 0.0

    def external_volume(self):
        z_cut = self.r_min / np.tan(self.theta)
        z_min = min(self.z)
        z_max = max(self.z)
        if np.sign(z_min * z_max) <= 0:
            external_volume = _coaxial_cone_volume(z_max + z_cut, self.theta, self.phi, self.z_offset)
            external_volume += _coaxial_cone_volume(-z_min + z_cut, self.theta, self.phi, self.z_offset)
            external_volume -= 2 * _coaxial_cone_volume(z_cut, self.theta, self.phi, self.z_offset)
        else:
            z_cut += min(abs(z_max), abs(z_min))
            external_volume = _coaxial_cone_volume(z_max - z_min + z_cut, self.theta, self.phi, self.z_offset)
            external_volume -= _coaxial_cone_volume(z_cut, self.theta, self.phi, self.z_offset)
        return external_volume

    def inner_surface_area(self):
        return 0.0

    def external_surface_area(self):
        z_cut = self.r_min / np.tan(self.theta)
        z_min = min(self.z)
        z_max = max(self.z)
        if np.sign(z_min * z_max) <= 0:
            surface_area = _coaxial_cone_surface_area(z_max + z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_surface_area(-z_min + z_cut, self.theta, self.phi, self.z_offset)
            surface_area -= 2 * _coaxial_cone_surface_area(z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_base_surface_area(z_max + z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_base_surface_area(-z_min + z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_side_cut_area(z_max + z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_side_cut_area(-z_min + z_cut, self.theta, self.phi, self.z_offset)
            surface_area -= 2 * _coaxial_cone_side_cut_area(z_cut, self.theta, self.phi, self.z_offset)
        else:
            z_cut += min(abs(z_max), abs(z_min))
            surface_area = _coaxial_cone_surface_area(z_max - z_min + z_cut, self.theta, self.phi, self.z_offset)
            surface_area -= _coaxial_cone_surface_area(z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_base_surface_area(z_max - z_min + z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_base_surface_area(z_cut, self.theta, self.phi, self.z_offset)
            surface_area += _coaxial_cone_side_cut_area(z_max - z_min + z_cut, self.theta, self.phi, self.z_offset)
            surface_area -= _coaxial_cone_side_cut_area(z_cut, self.theta, self.phi, self.z_offset)
        return surface_area


def _cone_volume(h, theta, phi):
    if h > 0:
        r = h * np.tan(theta)
        return 1/6 * r**2 * h * phi
    else:
        return 0


def _cone_surface_area(h, theta, phi):
    if h > 0:
        r = h * np.tan(theta)
        return phi * np.sqrt(r**2 + h**2) / 2
    else:
        return 0


def _cone_base_surface_area(h, theta, phi):
    if h > 0:
        r = h * np.tan(theta)
        return phi / 2 * r**2
    else:
        return 0


def _cone_side_cut_area(h, theta, phi):
    if h > 0 and phi < 2*np.pi:
        r = h * np.tan(theta)
        return r * h
    else:
        return 0


def _coaxial_cone_volume(h, theta, phi, z_offset):
    return _cone_volume(h, theta, phi) - _cone_volume(h - z_offset, theta, phi)


def _coaxial_cone_surface_area(h, theta, phi, z_offset):
    return _cone_surface_area(h, theta, phi) + _cone_surface_area(h - z_offset, theta, phi)


def _coaxial_cone_base_surface_area(h, theta, phi, z_offset):
    return _cone_base_surface_area(h, theta, phi) - _cone_base_surface_area(h - z_offset, theta, phi)


def _coaxial_cone_side_cut_area(h, theta, phi, z_offset):
    return _cone_side_cut_area(h, theta, phi) - _cone_side_cut_area(h - z_offset, theta, phi)
