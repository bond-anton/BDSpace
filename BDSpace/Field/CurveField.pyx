import numpy as np

from cython import boundscheck, wraparound
from cython.parallel import prange

from cpython.array cimport array, clone
from libc.math cimport sqrt

from .Field cimport Field
from BDSpace.Curve.Parametric cimport ParametricCurve


cdef class CurveField(Field):

    def __init__(self, str name, str field_type, ParametricCurve curve):
        super(CurveField, self).__init__(name, field_type)
        self.__curve = curve
        self.__tree_mesh = self.__curve.mesh_tree()
        self.__flat_mesh = self.__tree_mesh.flatten()
        self.__curve.add_element(self)
        self.__a = 0.0

    @property
    def curve(self):
        return self.__curve

    @curve.setter
    def curve(self, ParametricCurve curve):
        self.__curve.remove_element(self)
        self.__curve = curve
        self.__tree_mesh = self.__curve.mesh_tree()
        self.__flat_mesh = self.__tree_mesh.flatten()
        self.__curve.add_element(self)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, double a):
        self.__a = a

    cdef double __linear_density_point(self, double t) nogil:
        return self.__a

    cpdef double linear_density_point(self, double t):
        return self.__linear_density_point(t)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] linear_density(self, double[:] t):
        cdef:
            int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                result[i] = self.__linear_density_point(t[i])
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
            double[:, :] global_xyz = self.to_global_coordinate_system(xyz)
            double[:, :] curve_xyz = self.__curve.to_local_coordinate_system(global_xyz)
            double[:] t = self.__flat_mesh.physical_nodes
            double[:, :] curve_points = self.__curve.generate_points(t)
            double[:] nl = self.linear_density(t)
            double[:] dl = self.__flat_mesh.solution
            int i, j, s = xyz.shape[0], ms = self.__flat_mesh.num
            double d, x, y, z
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        with nogil:
            for i in prange(s):
                values[i] = 0.0
                for j in prange(ms):
                    x = curve_xyz[i, 0] - curve_points[j, 0]
                    y = curve_xyz[i, 1] - curve_points[j, 1]
                    z = curve_xyz[i, 2] - curve_points[j, 2]
                    d = sqrt(x * x + y * y + z * z)
                    if d < self.__r:
                        d = self.__r
                    values[i] += nl[j] * dl[j] / d
        return values

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        cdef:
            double[:, :] curve_xyz = self.__coordinate_system.to_parent(xyz)
            double[:] t = self.__flat_mesh.physical_nodes
            double[:, :] curve_points = self.__curve.generate_points(t)
            double[:] nl = self.linear_density(t)
            double[:] dl = self.__flat_mesh.solution
            int i, j, s = xyz.shape[0], ms = self.__flat_mesh.num
            double d, d2, d2_min = self.__r * self.__r
            double x, y, z
            double[:, :] values = np.empty((s, 3), dtype=np.double)
        with nogil:
            for i in prange(s):
                values[i, 0] = 0.0
                values[i, 1] = 0.0
                values[i, 2] = 0.0
                for j in prange(ms):
                    x = curve_xyz[i, 0] - curve_points[j, 0]
                    y = curve_xyz[i, 1] - curve_points[j, 1]
                    z = curve_xyz[i, 2] - curve_points[j, 2]
                    d2 = x * x + y * y + z * z
                    if d2 >= d2_min:
                        d = sqrt(d2)
                        values[i, 0] += nl[j] * dl[j] * x / d2 / d
                        values[i, 1] += nl[j] * dl[j] * y / d2 / d
                        values[i, 2] += nl[j] * dl[j] * z / d2 / d
        return values
