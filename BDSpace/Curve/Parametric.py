from __future__ import division, print_function
import numpy as np

from BDSpace.Curve import Curve
from ._helpers import check_equation


class ParametricCurve(Curve):

    def __init__(self, name='Parametric curve', coordinate_system=None,
                 x=None, y=None, z=None, start=0.0, stop=0.0):
        super(ParametricCurve, self).__init__(name, coordinate_system=coordinate_system)
        self.__x = None
        self.__y = None
        self.__z = None
        self.x = x
        self.y = y
        self.z = z
        self.__start = None
        self.start = start
        self.__stop = None
        self.stop = stop

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, x):
        if x is None:
            self.__x = lambda t: None
        elif check_equation(x):
            self.__x = x
        else:
            raise ValueError('Parametric equation must be single argument callable compatible with 1D numpy arrays')

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, y):
        if y is None:
            self.__y = lambda t: None
        elif check_equation(y):
            self.__y = y
        else:
            raise ValueError('Parametric equation must be single argument callable compatible with 1D numpy arrays')

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, z):
        if z is None:
            self.__z = lambda t: None
        elif check_equation(z):
            self.__z = z
        else:
            raise ValueError('Parametric equation must be single argument callable compatible with 1D numpy arrays')

    @property
    def start(self):
        return self.__start

    @start.setter
    def start(self, start):
        self.__start = np.float64(start)

    @property
    def stop(self):
        return self.__stop

    @stop.setter
    def stop(self, stop):
        self.__stop = np.float64(stop)

    def generate_points(self, t):
        xyz = np.zeros((len(t), 3), dtype=np.float64)
        xyz[:, 0] = self.x(t)
        xyz[:, 1] = self.y(t)
        xyz[:, 2] = self.z(t)
        return xyz

    def tangent(self, t):
        xyz = self.generate_points(t)
        return np.diff(xyz, axis=0) / np.diff(t).reshape(len(t) - 1, 1)

    def length(self, a=None, b=None, precision=1e-6, return_details=False):
        if a is None:
            a = self.start
        if b is None:
            b = self.stop
        num_points = int(100 * abs(b-a) / (2 * np.pi)) + 1
        iteration = 0
        while True:
            iteration += 1
            t = np.linspace(a, b, num=num_points, endpoint=True, dtype=np.float)
            xyz = self.tangent(t)
            dl = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
            length = np.trapz(dl, t[1:])
            xyz = np.diff(self.generate_points(t), axis=0)
            length_polygonal = np.sum(np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2))
            if length == length_polygonal == 0.0:
                error = 0.0
                break
            error = abs(length - length_polygonal) / max(length, length_polygonal)
            if error <= abs(precision):
                break
            num_points *= 2
        if return_details:
            return max(length, length_polygonal), error, iteration
        return max(length, length_polygonal)


class Line(ParametricCurve):

    def __init__(self, name='Line', coordinate_system=None, origin=(0, 0, 0), a=1, b=1, c=1, start=0, stop=1):
        self.origin = np.array(origin, dtype=np.float)
        self.__a = None
        self.__b = None
        self.__c = None
        self.a = a
        self.b = b
        self.c = c
        super(Line, self).__init__(name=name, coordinate_system=coordinate_system,
                                   x=lambda t: self.origin[0] + self.a * t,
                                   y=lambda t: self.origin[1] + self.b * t,
                                   z=lambda t: self.origin[2] + self.c * t,
                                   start=start, stop=stop)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, a):
        self.__a = np.float64(a)

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, b):
        self.__b = np.float64(b)

    @property
    def c(self):
        return self.__c

    @c.setter
    def c(self, c):
        self.__c = np.float64(c)


class Arc(ParametricCurve):

    def __init__(self, name='Arc', coordinate_system=None, a=1, b=1, start=0, stop=np.pi * 2, right=True):
        self.__a = None
        self.__b = None
        self.a = max(a, b)
        self.b = min(a, b)
        direction = 1 if right else -1
        super(Arc, self).__init__(name=name, coordinate_system=coordinate_system,
                                  x=lambda t: self.a * np.cos(t),
                                  y=lambda t: direction * self.b * np.sin(t),
                                  z=lambda t: np.zeros_like(t),
                                  start=start, stop=stop)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, a):
        self.__a = np.float64(a)

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, b):
        self.__b = np.float64(b)

    def get_eccentricity(self):
        return np.sqrt((self.a**2 - self.b**2) / self.a**2)

    def get_focus(self):
        return self.a * self.get_eccentricity()


class Helix(ParametricCurve):

    def __init__(self, name='Helix', coordinate_system=None, radius=1, pitch=1, start=0, stop=10, right=True):
        self.__radius = None
        self.__pitch = None
        self.radius = radius
        self.pitch = pitch
        direction = -1 if right else 1
        super(Helix, self).__init__(name=name, coordinate_system=coordinate_system,
                                    x=lambda t: self.radius - self.radius * np.cos(t),
                                    y=lambda t: direction * self.radius * np.sin(t),
                                    z=lambda t: self.pitch / (2 * np.pi) * t,
                                    start=start, stop=stop)

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self, radius):
        self.__radius = np.float64(radius)

    @property
    def pitch(self):
        return self.__pitch

    @pitch.setter
    def pitch(self, pitch):
        self.__pitch = np.float64(pitch)
