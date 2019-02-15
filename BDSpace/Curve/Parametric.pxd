from BDMesh.Mesh1DUniform cimport Mesh1DUniform
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform
from BDSpace.Space cimport Space


cdef class ParametricCurve(Space):
    cdef:
        double __start
        double __stop
    cpdef double x_point(self, double t)
    cpdef double y_point(self, double t)
    cpdef double z_point(self, double t)
    cpdef double[:] x(self, double[:] t)
    cpdef double[:] y(self, double[:] t)
    cpdef double[:] z(self, double[:] t)
    cpdef double[:, :] generate_points(self, double[:] t)
    cpdef double tangent_x_point(self, double t, double step=*, bint left=*, bint right=*)
    cpdef double tangent_y_point(self, double t, double step=*, bint left=*, bint right=*)
    cpdef double tangent_z_point(self, double t, double step=*, bint left=*, bint right=*)
    cpdef double[:] tangent_x(self, double[:] t, double step=*)
    cpdef double[:] tangent_y(self, double[:] t, double step=*)
    cpdef double[:] tangent_z(self, double[:] t, double step=*)
    cpdef double[:, :] tangent(self, double[:] t, double step=*)
    cdef double __length_tangent_array(self, double[:] t, double tangent_step=*)
    cdef double __length_poly_array(self, double[:] t)
    cdef double __length_tangent_mesh(self, Mesh1DUniform mesh, double tangent_step=*)
    cpdef double length(self, double precision=*, unsigned int max_iterations=*, double tangent_step=*)
    cpdef TreeMesh1DUniform mesh_tree(self, double precision=*, unsigned int max_iterations=*, double tangent_step=*)
    cpdef double distance_to_point(self, double t, double[:] xyz)


cdef class Line(ParametricCurve):
    cdef:
        double[:] __origin
        double __a
        double __b
        double __c


cdef class Arc(ParametricCurve):
    cdef:
        double __a
        double __b
        short __direction
    cpdef double eccentricity(self)
    cpdef double focus(self)


cdef class Helix(ParametricCurve):
    cdef:
        double __radius
        double __pitch
        short __direction
