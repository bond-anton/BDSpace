import numpy as np
from cpython.array cimport array, clone

from .Field cimport Field
from BDSpace.Curve.Parametric cimport ParametricCurve
from BDMesh.Mesh1D cimport Mesh1D


cdef class CurveField(Field):

    def __init__(self, str name, str field_type, ParametricCurve curve):
        self.__curve = curve
        self.__tree_mesh = self.__curve.mesh_tree()
        super(CurveField, self).__init__(name, field_type)

    @property
    def curve(self):
        return self.__curve

    @curve.setter
    def curve(self, ParametricCurve curve):
        self.__curve = curve
        self.__tree_mesh = self.__curve.mesh_tree()

    cpdef double linear_density(self, double t):
        return 0.0


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

    cpdef double[:] scalar_field(self, double[:, :] xyz):
        cdef:
            Mesh1D flat_mesh = self.__tree_mesh.flatten()
            unsigned int i, j, s = xyz.shape[0]
            double t, d
            array[double] values, template = array('d')
        values = clone(template, s, zero=False)
        for i in range(s):
            values[i] = 0.0
            for j in range(flat_mesh.num):
                t = flat_mesh.__physical_boundary_1 + flat_mesh.j() * flat_mesh.__local_nodes[j]
                d = self.__curve.distance_to_point(t, xyz[i])
                if d < self.__r:
                    d = self.__r
                values[i] += self.linear_density(t) * flat_mesh.solution[j] / d
        return values

    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        cdef:
            Mesh1D flat_mesh = self.__tree_mesh.flatten()
            unsigned int i, j, s = xyz.shape[0]
            double t, d, d2, d2_min = self.__r * self.__r
            double[:, :] values = np.empty((s, 3), dtype=np.double)
        for i in range(s):
            values[i, 0] = 0.0
            values[i, 1] = 0.0
            values[i, 2] = 0.0
            for j in range(flat_mesh.num):
                t = flat_mesh.__physical_boundary_1 + flat_mesh.j() * flat_mesh.__local_nodes[j]
                d = self.__curve.distance_to_point(t, xyz[i])
                d2 = d * d
                if d2 < d2_min:
                    d2 = d2_min
                # values[i, 0] += self.linear_density(t) * flat_mesh.solution[j] / d2
                # values[i, 1] += self.linear_density(t) * flat_mesh.solution[j] / d2
                # values[i, 2] += self.linear_density(t) * flat_mesh.solution[j] / d2
                # HINT: v = xyz -L(t)
        return values