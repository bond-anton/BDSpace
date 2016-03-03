from __future__ import division
import numpy as np

from Space.Curve import Curve


class ParametricCurve(Curve):

    def __init__(self, name='Parametric curve', coordinate_system=None,
                 x=None, y=None, z=None, start=0.0, stop=0.0):
        super(ParametricCurve, self).__init__(name, coordinate_system=coordinate_system)
        self.x = None
        self.y = None
        self.z = None
        self.set_equations(x, y, z)
        self.start = start
        self.stop = stop

    def set_equations(self, x, y, z):
        if x is None or y is None or z is None:
            self.x = lambda t: None
            self.y = lambda t: None
            self.z = lambda t: None
        elif not check_equations([x, y, z]):
            raise ValueError('Parametric equations must be single argument callable compatible with 1D numpy arrays')
        else:
            self.x = x
            self.y = y
            self.z = z

    def generate_points(self, t):
        xyz = np.zeros((len(t), 3), dtype=np.float)
        xyz[:, 0] = self.x(t)
        xyz[:, 1] = self.y(t)
        xyz[:, 2] = self.z(t)
        return xyz

    def tangent(self, t):
        xyz = self.generate_points(t)
        return np.diff(xyz, axis=0)

    def length(self, a=None, b=None):
        if a is None:
            a = self.start
        if b is None:
            b = self.stop
        t = np.linspace(a, b, num=500 * abs(b-a) / (2 * np.pi), endpoint=True, dtype=np.float)
        xyz = self.tangent(t)
        dl = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
        length = np.trapz(dl, t[1:])
        return length


class Helix(ParametricCurve):

    def __init__(self, name='Helix', coordinate_system=None, r=1, h=1, start=0, stop=10, right=True):
        self.r = r
        self.h = h
        direction = 1 if right else -1
        super(Helix, self).__init__(name=name, coordinate_system=coordinate_system,
                                    x=lambda t: self.r * np.cos(t) - self.r,
                                    y=lambda t: direction * self.r * np.sin(t),
                                    z=lambda t: self.h / (2 * np.pi) * t,
                                    start=start, stop=stop)


class Arc(ParametricCurve):

    def __init__(self, name='Arc', coordinate_system=None, r=1, start=0, stop=np.pi * 2, right=True):
        self.r = r
        direction = 1 if right else -1
        super(Arc, self).__init__(name=name, coordinate_system=coordinate_system,
                                  x=lambda t: direction * self.r * np.sin(t),
                                  y=lambda t: np.zeros_like(t),
                                  z=lambda t: -self.r * np.cos(t) + self.r,
                                  start=start, stop=stop)


def check_equations(equations):
        result = []
        for eq in equations:
            if hasattr(eq, '__call__'):
                parameter = np.arange(5, dtype=np.float)
                answer = eq(parameter)
                if isinstance(answer, np.ndarray):
                    if answer is None:
                        result.append(True)
                        continue
                    if answer.size == parameter.size and answer.shape == parameter.shape:
                        result.append(True)
                        continue
                result.append(False)
            else:
                result.append(False)
        return np.array(result).all()
