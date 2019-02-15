from .Field cimport Field
from BDSpace.Curve.Parametric cimport ParametricCurve
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

cdef class CurveField(Field):
    cdef:
        ParametricCurve __curve
        TreeMesh1DUniform __tree_mesh
        double __a

    cpdef double linear_density_point(self, double t)
    cpdef double[:] linear_density(self, double[:] t)


cdef class HyperbolicPotentialCurveConservativeField(CurveField):
    cdef:
        double __r
