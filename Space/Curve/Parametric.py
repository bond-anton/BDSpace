import numpy as np
import sympy as sym

from Space.Curve import Curve


class ParametricCurve(Curve):

    def __init__(self, name='Parametric curve', CS=None, arg=sym.symbols('t'),
                 x=None, y=None, z=None, start=None, stop=None):
        super(ParametricCurve, self).__init__(name, CS)
        self.arg = arg
        if x is None:
            self.x = self.arg * 0
        else:
            self.x = x
        if y is None:
            self.y = self.arg * 0
        else:
            self.y = y
        if z is None:
            self.z = self.arg * 0
        else:
            self.z = z
        if start is None:
            self.start = 0.0
        else:
            self.start = start
        if stop is None:
            self.stop = 0.0
        else:
            self.stop = stop

    def points(self, t):
        xyz = np.zeros((t.size, 3), dtype=np.float)
        fx = sym.lambdify(self.arg, self.x, modules='numpy')
        fy = sym.lambdify(self.arg, self.y, modules='numpy')
        fz = sym.lambdify(self.arg, self.z, modules='numpy')
        xyz[:, 0] = fx(t)
        xyz[:, 1] = fy(t)
        xyz[:, 2] = fz(t)
        return xyz

    def tangent(self, t):
        xyz = np.zeros((t.size, 3), dtype=np.float)
        fx = sym.lambdify(self.arg, sym.diff(self.x, self.arg), modules='numpy')
        fy = sym.lambdify(self.arg, sym.diff(self.y, self.arg), modules='numpy')
        fz = sym.lambdify(self.arg, sym.diff(self.z, self.arg), modules='numpy')
        xyz[:, 0] = fx(t)
        xyz[:, 1] = fy(t)
        xyz[:, 2] = fz(t)
        return xyz

    def length(self, symbolic=False, a=None, b=None):
        if a is None:
            a = self.start
        if b is None:
            b = self.stop
        if symbolic:
            fx = sym.diff(self.x, self.arg)
            fy = sym.diff(self.y, self.arg)
            fz = sym.diff(self.z, self.arg)
            L = sym.N(sym.integrate(sym.sqrt(fx**2 + fy**2 + fz**2), (self.arg, a, b)))
        else:
            t = np.linspace(a, b, num=101 * abs(b-a) / (2 * np.pi), endpoint=True)
            xyz = self.tangent(t)
            dL = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
            L = np.trapz(dL, t)
        return L

class Helix(ParametricCurve):

    def __init__(self, name='Helix', CS=None, r=1, h=1, start=0, stop=10, right=True):
        self.r = r
        self.h = h
        direction = 1 if right else -1
        arg = sym.symbols('t')
        x = self.r * sym.cos(arg) - self.r
        y = direction * self.r * sym.sin(arg)
        z = self.h / (2 * sym.pi) * arg
        super(Helix, self).__init__(name, CS, arg, x, y, z, start, stop)

class Arc(ParametricCurve):

    def __init__(self, name='Arc', CS=None, r=1, start=0, stop=np.pi * 2, right=True):
        self.r = r
        direction = 1 if right else -1
        arg = sym.symbols('t')
        x = direction * self.r * sym.sin(arg)
        y = 0.0 # direction * self.r * sym.sin(arg)
        z = -self.r * sym.cos(arg) + self.r
        super(Arc, self).__init__(name, CS, arg, x, y, z, start, stop)