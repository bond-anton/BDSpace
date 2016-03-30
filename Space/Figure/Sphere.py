from __future__ import division, print_function
import numpy as np

from Space.Coordinates.transforms import reduce_angle
from Space.Figure import Figure


class SphericalShape(Figure):

    def __init__(self, name='Spherical shape', coordinate_system=None,
                 r_inner=0, r_outer=1.0, phi=np.pi/2):
        self.r_inner = max(min(r_inner, r_outer), 0.0)
        self.r_outer = max(max(r_inner, r_outer), 0.0)
        self.phi = reduce_angle(float(phi))
        super(SphericalShape, self).__init__(name, coordinate_system=coordinate_system)


class SphericalWedge(SphericalShape):

    def __init__(self, name='Spherical wedge', coordinate_system=None,
                 r_inner=0, r_outer=1.0, phi=np.pi/2, theta=np.array([0, np.pi/2])):
        max_theta = reduce_angle(max(theta))
        if max_theta > np.pi:
            max_theta = 2*np.pi - max_theta
        min_theta = reduce_angle(min(theta))
        if min_theta > np.pi:
            min_theta = 2*np.pi - min_theta
        theta_range = [min_theta, max_theta]
        self.theta = np.array(theta_range, dtype=np.float)
        super(SphericalWedge, self).__init__(name, coordinate_system=coordinate_system,
                                             r_inner=r_inner, r_outer=r_outer, phi=phi)

    def inner_volume(self):
        if self.phi == 2 * np.pi and (self.theta[1] - self.theta[0]) == np.pi:
            return 4/3 * np.pi * self.r_inner**3
        else:
            return 0

    def external_volume(self):
        if self.phi == 2 * np.pi and (self.theta[1] - self.theta[0]) == np.pi:
            return 4/3 * np.pi * self.r_outer**3
        else:
            v_inner = 1/3 * self.phi * self.r_inner**3 * (np.cos(self.theta[0]) - np.cos(self.theta[1]))
            v_outer = 1/3 * self.phi * self.r_outer**3 * (np.cos(self.theta[0]) - np.cos(self.theta[1]))
            return v_outer - v_inner

    def inner_surface_area(self):
        if self.phi == 2 * np.pi and (self.theta[1] - self.theta[0]) == np.pi:
            return 4 * np.pi * self.r_inner**2
        else:
            return 0

    def external_surface_area(self):
        if self.phi == 2 * np.pi and (self.theta[1] - self.theta[0]) == np.pi:
            return 4 * np.pi * self.r_outer**2
        else:
            s_inner = self.phi * self.r_inner**2 * (np.cos(self.theta[0] - np.cos(self.theta[1])))
            s_outer = self.phi * self.r_outer**2 * (np.cos(self.theta[0] - np.cos(self.theta[1])))
            if self.phi < 2 * np.pi:
                s_sides = (self.theta[1] - self.theta[0]) * (self.r_outer**2 - self.r_inner**2)
            else:
                s_sides = 0
            s_sides += (self.phi / 2) * (self.r_outer**2 - self.r_inner**2) * np.sin(self.theta[0])
            s_sides += (self.phi / 2) * (self.r_outer**2 - self.r_inner**2) * np.sin(self.theta[1])
            return s_inner + s_outer + s_sides


class SphericalCone(SphericalWedge):

    def __init__(self, name='Spherical cone', coordinate_system=None,
                 r_inner=0, r_outer=1.0, theta=np.pi/4):
        super(SphericalCone, self).__init__(name, coordinate_system=coordinate_system,
                                            r_inner=r_inner, r_outer=r_outer, phi=2*np.pi, theta=np.array([0, theta]))


class Sphere(SphericalCone):

    def __init__(self, name='Sphere', coordinate_system=None, r_inner=0, r_outer=1.0):
        super(Sphere, self).__init__(name, coordinate_system=coordinate_system,
                                     r_inner=r_inner, r_outer=r_outer, theta=np.pi)


class SphericalSegmentWedge(SphericalShape):

    def __init__(self, name='Spherical section', coordinate_system=None, r_inner=0, r_outer=1.0,
                 h1=0, h2=1.0, phi=np.pi/2):
        self.r_inner = max(min(r_inner, r_outer), 0.0)
        self.r_outer = max(max(r_inner, r_outer), 0.0)
        self.h1 = max(min(h1, h2), -self.r_outer)
        self.h2 = min(max(h1, h2), self.r_outer)
        super(SphericalSegmentWedge, self).__init__(name, coordinate_system=coordinate_system,
                                                    r_inner=r_inner, r_outer=r_outer, phi=phi)

    def inner_volume(self):
        if self.phi == 2 * np.pi and self.h1 <= - self.r_inner and self.h2 >= self.r_inner:
            return 4/3 * np.pi * self.r_inner**3
        else:
            return 0

    def external_volume(self):
        h_outer = self.h2 - self.h1
        v_outer = h_outer * (self.r_outer**2 - self.h1**2 - self.h1 * h_outer - 1/3 * h_outer**2) * self.phi / 2
        h1_inner = np.sign(self.h1) * min(self.r_inner, abs(self.h1))
        h2_inner = np.sign(self.h2) * min(self.r_inner, abs(self.h2))
        h_inner = h2_inner - h1_inner
        v_inner = h_inner * (self.r_inner**2 - h1_inner**2 - h1_inner * h_inner - 1/3 * h_inner**2) * self.phi / 2
        if self.phi == 2 * np.pi and self.h1 <= - self.r_inner and self.h2 >= self.r_inner:
            return v_outer
        else:
            return v_outer - v_inner

    def inner_surface_area(self):
        if self.phi == 2 * np.pi and self.h1 <= - self.r_inner and self.h2 >= self.r_inner:
            return 4 * np.pi * self.r_inner**2
        else:
            return 0

    def external_surface_area(self):
        theta1 = np.arccos(self.h1 / self.r_outer)
        theta2 = np.arccos(self.h2 / self.r_outer)
        h_outer = self.h2 - self.h1
        s_outer = self.phi * self.r_outer * h_outer
        h1_inner = np.sign(self.h1) * min(self.r_inner, abs(self.h1))
        h2_inner = np.sign(self.h2) * min(self.r_inner, abs(self.h2))
        theta1_inner = np.arccos(h1_inner / self.r_inner)
        theta2_inner = np.arccos(h2_inner / self.r_inner)
        h_inner = h2_inner - h1_inner
        s_inner = self.phi * self.r_inner * h_inner
        if np.sign(self.h1 * self.h2) >= 0:
            delta_theta = theta1 - theta2
            s_cut = self.r_outer**2 * (delta_theta - np.sin(delta_theta))
            s_cut += 2 * h_outer * self.r_outer * min(np.sin(theta1), np.sin(theta2))
            s_cut += h_outer * self.r_outer * abs(np.sin(theta2) - np.sin(theta1))
        else:
            delta_theta = theta1 - np.pi/2
            s_cut = self.r_outer**2 * (delta_theta - np.sin(delta_theta))
            s_cut += 2 * abs(self.h1) * self.r_outer * min(np.sin(theta1), np.sin(theta2))
            s_cut += abs(self.h1) * self.r_outer * abs(np.sin(theta2) - np.sin(theta1))
            delta_theta = np.pi/2 - theta2
            s_cut = self.r_outer**2 * (delta_theta - np.sin(delta_theta))
            s_cut += 2 * self.h2 * self.r_outer * min(np.sin(theta1), np.sin(theta2))
            s_cut += self.h2 * self.r_outer * abs(np.sin(theta2) - np.sin(theta1))
        if h_inner > 0:
            if np.sign(h1_inner * h2_inner) >= 0:
                delta_theta = theta1_inner - theta2_inner
                s_cut -= self.r_inner**2 * (delta_theta - np.sin(delta_theta))
                s_cut -= 2 * h_inner * self.r_inner * min(np.sin(theta1_inner), np.sin(theta2_inner))
                s_cut -= h_inner * self.r_inner * abs(np.sin(theta2_inner) - np.sin(theta1_inner))
            else:
                delta_theta = theta1_inner - np.pi/2
                s_cut -= self.r_inner**2 * (delta_theta - np.sin(delta_theta))
                s_cut -= 2 * abs(h1_inner) * self.r_inner * min(np.sin(theta1_inner), np.sin(theta2_inner))
                s_cut -= abs(h1_inner) * self.r_inner * abs(np.sin(theta2_inner) - np.sin(theta1_inner))
                delta_theta = np.pi/2 - theta2_inner
                s_cut = self.r_inner**2 * (delta_theta - np.sin(delta_theta))
                s_cut -= 2 * h2_inner * self.r_inner * min(np.sin(theta1_inner), np.sin(theta2_inner))
                s_cut -= h2_inner * self.r_inner * abs(np.sin(theta2_inner) - np.sin(theta1_inner))
        if self.phi == 2 * np.pi:
            s_cut = 0
        s_base = self.phi / 2 * (self.r_outer * np.sin(theta1))**2
        s_base += self.phi / 2 * (self.r_outer * np.sin(theta2))**2
        if h_inner > 0:
            s_base -= self.phi / 2 * (self.r_inner * np.sin(theta1_inner))**2
            s_base -= self.phi / 2 * (self.r_inner * np.sin(theta2_inner))**2
        if self.phi == 2 * np.pi and self.h1 <= - self.r_inner and self.h2 >= self.r_inner:
            return s_outer + s_base
        else:
            return s_outer + s_inner + s_cut + s_base


class SphericalSegment(SphericalSegmentWedge):

    def __init__(self, name='Spherical section', coordinate_system=None, r_inner=0, r_outer=1.0, h1=0, h2=1.0):
        super(SphericalSegment, self).__init__(name, coordinate_system=coordinate_system,
                                               r_inner=r_inner, r_outer=r_outer, h1=h1, h2=h2, phi=2*np.pi)


class SphericalCap(SphericalSegment):

    def __init__(self, name='Spherical section', coordinate_system=None, r_inner=0, r_outer=1.0, h1=0):
        super(SphericalCap, self).__init__(name, coordinate_system=coordinate_system,
                                           r_inner=r_inner, r_outer=r_outer, h1=h1, h2=r_outer)
