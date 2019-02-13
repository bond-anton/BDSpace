import numpy as np

from ._helpers import check_equation
from cython import boundscheck, wraparound

from cpython.array cimport array, clone

from BDSpace.Space cimport Space
from BDSpace.Coordinates cimport Cartesian


cdef class ParametricCurve(Space):

    def __init__(self, str name='Parametric curve', Cartesian coordinate_system=None,
                 double start=0.0, double stop=0.0):
        super(ParametricCurve, self).__init__(name, coordinate_system=coordinate_system)
        self.__start = None
        self.start = start
        self.__stop = None
        self.stop = stop

    cpdef double[:] x(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        return result

    cpdef double[:] y(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        return result

    cpdef double[:] z(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        return result

    @property
    def start(self):
        return self.__start

    @start.setter
    def start(self, double start):
        self.__start = start

    @property
    def stop(self):
        return self.__stop

    @stop.setter
    def stop(self, double stop):
        self.__stop = stop

    cpdef double[:, :] generate_points(self, double[:] t):
        cdef:
            double[:, :] xyz = np.empty((t.shape[0], 3), dtype=np.double)
        xyz[:, 0] = self.x(t)
        xyz[:, 1] = self.y(t)
        xyz[:, 2] = self.z(t)
        return xyz

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] tangent(self, double[:] t):
        cdef:
            double[:, :] xyz = self.generate_points(t)
            unsigned int i, s = t.shape[0]
            double[:, :] result = np.empty((s - 1, 3), dtype=np.double)
        for i in range(s - 1):
            result[i, 0] = (xyz[i + 1, 0] - xyz[i, 0]) / (t[i + 1] - t[i])
        return result

    cpdef double length(self, double a=None, double b=None, double precision=1e-6, bint print_details=False):
        cdef:
            unsigned int num_points, iteration = 0
        if a is None:
            a = self.__start
        if b is None:
            b = self.__stop
        num_points = int(100 * abs(b-a) / (2 * np.pi)) + 1
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
