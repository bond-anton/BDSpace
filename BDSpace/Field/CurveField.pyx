import numpy as np

from cython import boundscheck, wraparound

from cpython.array cimport array, clone
from libc.math cimport sqrt

from .Field cimport Field
from BDSpace.Curve.Parametric cimport ParametricCurve
from BDMesh.Mesh1D cimport Mesh1D


cdef class CurveField(Field):

    def __init__(self, str name, str field_type, ParametricCurve curve):
        self.__curve = curve
        self.__tree_mesh = self.__curve.mesh_tree()
        self.__curve.add_element(self)
        self.__a = 0.0
        super(CurveField, self).__init__(name, field_type)

    @property
    def curve(self):
        return self.__curve

    @curve.setter
    def curve(self, ParametricCurve curve):
        self.__curve.remove_element(self)
        self.__curve = curve
        self.__tree_mesh = self.__curve.mesh_tree()
        self.__curve.add_element(self)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, double a):
        self.__a = a

    cpdef double linear_density_point(self, double t):
        return self.__a

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] linear_density(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.linear_density_point(t[i])
        return result


cdef class HyperbolicPotentialCurveConservativeField(CurveField):

    def __init__(self, str name, str field_type, ParametricCurve curve, double r):
        self.__r = r
        super(HyperbolicPotentialCurveConservativeField, self).__init__(name, field_type, curve)

    @property
    def r(self):
        return self.__r

    @r.setter
    def r(self, double r):
        self.__r = r

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field(self, double[:, :] xyz):
        cdef:
            Mesh1D flat_mesh = self.__tree_mesh.flatten()
            double[:] t = flat_mesh.physical_nodes
            double[:, :] curve_points = self.__curve.generate_points(t)
            double[:] nl = self.linear_density(t)
            unsigned int i, j, s = xyz.shape[0], ms = flat_mesh.num
            double d
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        for i in range(s):
            values[i] = 0.0
            for j in range(ms):
                x = xyz[i, 0] - curve_points[j, 0]
                y = xyz[i, 1] - curve_points[j, 1]
                z = xyz[i, 2] - curve_points[j, 2]
                d = sqrt(x * x + y * y + z * z)
                if d < self.__r:
                    d = self.__r
                values[i] += nl[j] * flat_mesh.solution[j] / d
        return values

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        cdef:
            Mesh1D flat_mesh = self.__tree_mesh.flatten()
            double[:] t = flat_mesh.physical_nodes
            double[:, :] curve_points = self.__curve.generate_points(t)
            double[:] nl = self.linear_density(t)
            unsigned int i, j, s = xyz.shape[0], ms = flat_mesh.num
            double d2, d2_min = self.__r * self.__r
            double x, y, z
            double[:, :] values = np.empty((s, 3), dtype=np.double)
        for i in range(s):
            values[i, 0] = 0.0
            values[i, 1] = 0.0
            values[i, 2] = 0.0
            for j in range(ms):
                x = xyz[i, 0] - curve_points[j, 0]
                y = xyz[i, 1] - curve_points[j, 1]
                z = xyz[i, 2] - curve_points[j, 2]
                d2 = x * x + y * y + z * z
                if d2 < d2_min:
                    d2 = d2_min
                values[i, 0] += nl[j] * flat_mesh.solution[j] * x / d2
                values[i, 1] += nl[j] * flat_mesh.solution[j] * y / d2
                values[i, 2] += nl[j] * flat_mesh.solution[j] * z / d2
        return values
