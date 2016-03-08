from __future__ import division
import numpy as np

from Space.Figure import Figure


class SphericalShape(Figure):

    def __init__(self, name='Spherical shape', coordinate_system=None,
                 r_inner=0, r_outer=1.0, phi=np.pi/2):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.phi = phi
        super(SphericalShape, self).__init__(name, coordinate_system=coordinate_system)


class SphericalWedge(SphericalShape):
    def __init__(self, name='Spherical wedge', coordinate_system=None,
                 r_inner=0, r_outer=1.0, phi=np.pi/2, theta=np.pi/2):
        self.theta = theta
        super(SphericalWedge, self).__init__(name, coordinate_system=coordinate_system,
                                             r_inner=r_inner, r_outer=r_outer, phi=phi)


class SphericalCone(SphericalWedge):

    def __init__(self, name='Spherical cone', coordinate_system=None,
                 r_inner=0, r_outer=1.0, theta=np.pi/4):
        super(SphericalCone, self).__init__(name, coordinate_system=coordinate_system,
                                            r_inner=r_inner, r_outer=r_outer, phi=2*np.pi, theta=theta)


class Sphere(SphericalCone):

    def __init__(self, name='Sphere', coordinate_system=None, r_inner=0, r_outer=1.0):
        super(Sphere, self).__init__(name, coordinate_system=coordinate_system,
                                     r_inner=r_inner, r_outer=r_outer, theta=np.pi)


class SphericalSectionWedge(SphericalShape):

    def __init__(self, name='Spherical section', coordinate_system=None, r_inner=0, r_outer=1.0,
                 h1=0, h2=1.0, phi=np.pi/2):
        self.h1 = h1
        self.h2 = h2
        super(SphericalSectionWedge, self).__init__(name, coordinate_system=coordinate_system,
                                                    r_inner=r_inner, r_outer=r_outer, phi=phi)


class SphericalSection(SphericalSectionWedge):

    def __init__(self, name='Spherical section', coordinate_system=None, r_inner=0, r_outer=1.0, h1=0, h2=1.0):
        super(SphericalSection, self).__init__(name, coordinate_system=coordinate_system,
                                               r_inner=r_inner, r_outer=r_outer, h1=h1, h2=h2, phi=2*np.pi)


class SphericalCap(SphericalSection):

    def __init__(self, name='Spherical section', coordinate_system=None, r_inner=0, r_outer=1.0, h1=0):
        super(SphericalCap, self).__init__(name, coordinate_system=coordinate_system,
                                           r_inner=r_inner, r_outer=r_outer, h1=h1, h2=r_outer)
