from .Field cimport Field
from BDSpace.Curve.Parametric cimport ParametricCurve
from BDMesh.Mesh1D cimport Mesh1D
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

cdef class CurveField(Field):
    cdef:
        ParametricCurve __curve
        TreeMesh1DUniform __tree_mesh
        Mesh1D __flat_mesh
        double __a

    cdef double __linear_density_point(self, double t) nogil
    cpdef double linear_density_point(self, double t)
    cpdef double[:] linear_density(self, double[:] t)


cdef class HyperbolicPotentialCurveConservativeField(CurveField):
    cdef:
        double __r
