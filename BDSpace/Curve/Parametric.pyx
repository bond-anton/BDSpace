import numpy as np

from cython import boundscheck, wraparound

from cpython.array cimport array, clone
from libc.math cimport sin, cos, sqrt, M_PI
from BDMesh.Mesh1D cimport Mesh1D
from BDMesh.Mesh1DUniform cimport Mesh1DUniform
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

from BDSpace.Space cimport Space
from BDSpace.Coordinates.Cartesian cimport Cartesian
from ._helpers cimport trapz_1d, refinement_points


cdef class ParametricCurve(Space):

    def __init__(self, str name='Parametric curve', Cartesian coordinate_system=None,
                 double start=0.0, double stop=0.0):
        super(ParametricCurve, self).__init__(name, coordinate_system=coordinate_system)
        self.__start = start
        self.__stop = stop
        self.__dt = 1.0e-10
        self.__precision = 1.0e-6

    cpdef double x_point(self, double t):
        return 0.0

    cpdef double y_point(self, double t):
        return 0.0

    cpdef double z_point(self, double t):
        return 0.0

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] x(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        for i in range(s):
            result[i] = self.x_point(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] y(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        for i in range(s):
            result[i] = self.y_point(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] z(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        for i in range(s):
            result[i] = self.z_point(t[i])
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

    @property
    def dt(self):
        return self.__dt

    @dt.setter
    def dt(self, double dt):
        self.__dt = dt

    @property
    def precision(self):
        return self.__precision

    @precision.setter
    def precision(self, double precision):
        self.__precision = precision

    cpdef double[:, :] generate_points(self, double[:] t):
        cdef:
            double[:, :] xyz = np.empty((t.shape[0], 3), dtype=np.double)
        xyz[:, 0] = self.x(t)
        xyz[:, 1] = self.y(t)
        xyz[:, 2] = self.z(t)
        return xyz

    cpdef double tangent_x_point(self, double t, bint left=True, bint right=True):
        cdef:
            double step2 = self.__dt / 2
        if left and right:
            return (self.x_point(t + step2) - self.x_point(t - step2)) / self.__dt
        elif left:
            return (self.x_point(t) - self.x_point(t - self.__dt)) / self.__dt
        else:
            return (self.x_point(t + self.__dt) - self.x_point(t)) / self.__dt

    cpdef double tangent_y_point(self, double t, bint left=True, bint right=True):
        cdef:
            double step2 = self.__dt / 2
        if left and right:
            return (self.y_point(t + step2) - self.y_point(t - step2)) / self.__dt
        elif left:
            return (self.y_point(t) - self.y_point(t - self.__dt)) / self.__dt
        else:
            return (self.y_point(t + self.__dt) - self.y_point(t)) / self.__dt

    cpdef double tangent_z_point(self, double t, bint left=True, bint right=True):
        cdef:
            double step2 = self.__dt / 2
        if left and right:
            return (self.z_point(t + step2) - self.z_point(t - step2)) / self.__dt
        elif left:
            return (self.z_point(t) - self.z_point(t - self.__dt)) / self.__dt
        else:
            return (self.z_point(t + self.__dt) - self.z_point(t)) / self.__dt

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] tangent_x(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0] - 1
            array[double] result, template = array('d')
        result = clone(template, s + 1, zero=False)
        result[0] = self.tangent_x_point(t[0], left=False)
        result[s] = self.tangent_x_point(t[s], right=False)
        for i in range(1, s):
            result[i] = self.tangent_x_point(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] tangent_y(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0] - 1
            array[double] result, template = array('d')
        result = clone(template, s + 1, zero=False)
        result[0] = self.tangent_y_point(t[0], left=False)
        result[s] = self.tangent_y_point(t[s], right=False)
        for i in range(1, s):
            result[i] = self.tangent_y_point(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] tangent_z(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0] - 1
            array[double] result, template = array('d')
        result = clone(template, s + 1, zero=False)
        result[0] = self.tangent_z_point(t[0], left=False)
        result[s] = self.tangent_z_point(t[s], right=False)
        for i in range(1, s):
            result[i] = self.tangent_z_point(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] tangent(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0] - 1
            double[:, :] result = np.empty((s + 1, 3), dtype=np.double)
        result[0, 0] = self.tangent_x_point(t[0], left=False)
        result[0, 1] = self.tangent_y_point(t[0], left=False)
        result[0, 2] = self.tangent_z_point(t[0], left=False)
        result[s, 0] = self.tangent_x_point(t[s], right=False)
        result[s, 1] = self.tangent_y_point(t[s], right=False)
        result[s, 2] = self.tangent_z_point(t[s], right=False)
        for i in range(1, s):
            result[i, 0] = self.tangent_x_point(t[s])
            result[i, 1] = self.tangent_y_point(t[s])
            result[i, 2] = self.tangent_z_point(t[s])
        return result

    @boundscheck(False)
    @wraparound(False)
    cdef double __length_tangent_array(self, double[:] t):
        cdef:
            unsigned int i, num_points = t.shape[0]
            double[:, :] xyz = self.tangent(t)
            array[double] dl, template = array('d')
        dl = clone(template, num_points, zero=False)
        for i in range(num_points):
            dl[i] = sqrt(xyz[i, 0] * xyz[i, 0] + xyz[i, 1] * xyz[i, 1] + xyz[i, 2] * xyz[i, 2])
        return trapz_1d(dl, t)

    @boundscheck(False)
    @wraparound(False)
    cdef double __length_poly_array(self, double[:] t):
        cdef:
            unsigned int i, num_points = t.shape[0] - 1
            double[:, :] xyz = self.generate_points(t)
            double result = 0.0, dx, dy, dz
        for i in range(num_points):
            dx = xyz[i + 1, 0] - xyz[i, 0]
            dy = xyz[i + 1, 1] - xyz[i, 1]
            dz = xyz[i + 1, 2] - xyz[i, 2]
            result += sqrt(dx * dx + dy * dy + dz * dz)
        return result

    @boundscheck(False)
    @wraparound(False)
    cdef double __length_tangent_mesh(self, Mesh1DUniform mesh):
        cdef:
            unsigned int i, num_points = mesh.num
            double[:, :] xyz = self.generate_points(mesh.physical_nodes)
            double[:, :] xyz_t = self.tangent(mesh.physical_nodes)
            double result_t, result_p, result_t_acc = 0.0, result_p_acc = 0.0, dx, dy, dz, dl1, dl2
            array[double] solution, error, template = array('d')
        solution = clone(template, num_points, zero=False)
        error = clone(template, num_points, zero=False)
        solution[0] = 0.0
        error[0] = 0.0
        for i in range(num_points - 1):
            dx = xyz[i + 1, 0] - xyz[i, 0]
            dy = xyz[i + 1, 1] - xyz[i, 1]
            dz = xyz[i + 1, 2] - xyz[i, 2]
            result_p = sqrt(dx * dx + dy * dy + dz * dz)
            result_p_acc += result_p
            dl1 = sqrt(xyz_t[i, 0] * xyz_t[i, 0] + xyz_t[i, 1] * xyz_t[i, 1] + xyz_t[i, 2] * xyz_t[i, 2])
            dl2 = sqrt(xyz_t[i + 1, 0] * xyz_t[i + 1, 0] + xyz_t[i + 1, 1] * xyz_t[i + 1, 1]\
                       + xyz_t[i + 1, 2] * xyz_t[i + 1, 2])
            result_t = (mesh.physical_nodes[i + 1] - mesh.physical_nodes[i]) * (dl1 + dl2) / 2
            result_t_acc += result_t
            solution[i + 1] = result_t
            error[i + 1] = abs(result_t - result_p)
        mesh.solution = solution
        mesh.residual = error
        return result_t_acc

    @boundscheck(False)
    @wraparound(False)
    cpdef double length(self, unsigned int max_iterations=100):
        cdef:
            TreeMesh1DUniform meshes_tree = self.mesh_tree(max_iterations)
            Mesh1D flat_mesh = meshes_tree.flatten()
            unsigned int i
            double length_tangent = 0
        for i in range(flat_mesh.num):
            length_tangent += flat_mesh.solution[i]
        return length_tangent

    @boundscheck(False)
    @wraparound(False)
    cpdef TreeMesh1DUniform mesh_tree(self, unsigned int max_iterations=100):
        cdef:
            Mesh1DUniform root_mesh, mesh, refinement_mesh
            TreeMesh1DUniform meshes_tree
            unsigned int num_points = 3, iteration = 0, level, i, to_refine
            double[:] t
            double length_tangent = 0.0
            long[:, :] refinements
        root_mesh = Mesh1DUniform(self.__start, self.__stop,
                                  boundary_condition_1=0.0,
                                  boundary_condition_2=0.0,
                                  physical_step=(self.__stop - self.__start) / 2.0)
        meshes_tree = TreeMesh1DUniform(root_mesh, refinement_coefficient=2, aligned=True)
        while iteration < max_iterations:
            iteration += 1
            level = max(meshes_tree.levels)
            to_refine = 0
            for mesh in meshes_tree.__tree[level]:
                length_tangent = self.__length_tangent_mesh(mesh)
                refinements = refinement_points(mesh, self.__precision)
                for i in range(refinements.shape[0]):
                    to_refine += 1
                    refinement_mesh = Mesh1DUniform(
                        mesh.to_physical(np.array([mesh.__local_nodes[refinements[i][0]]]))[0],
                        mesh.to_physical(np.array([mesh.__local_nodes[refinements[i][1]]]))[0],
                        boundary_condition_1=0.0,
                        boundary_condition_2=0.0,
                        physical_step=mesh.physical_step/meshes_tree.refinement_coefficient)
                    meshes_tree.add_mesh(refinement_mesh)
                meshes_tree.remove_coarse_duplicates()
            if to_refine == 0:
                break
        return meshes_tree

    cpdef double distance_to_point_square(self, double t, double[:] xyz):
        cdef:
            double x, y, z
        x = self.x_point(t) - xyz[0]
        y = self.y_point(t) - xyz[1]
        z = self.z_point(t) - xyz[2]
        return x*x + y*y + z*z

    cpdef double distance_to_point(self, double t, double[:] xyz):
        return sqrt(self.distance_to_point_square(t, xyz))


cdef class Line(ParametricCurve):

    def __init__(self, str name='Line', Cartesian coordinate_system=None,
                 double[:] origin=np.zeros(3, dtype=np.double),
                 double a=1.0, double b=1.0, double c=1.0,
                 double start=0.0, double stop=1.0):
        self.__origin = origin
        self.__a = a
        self.__b = b
        self.__c = c
        super(Line, self).__init__(name=name, coordinate_system=coordinate_system,
                                   start=start, stop=stop)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, double a):
        self.__a = a

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, double b):
        self.__b = b

    @property
    def c(self):
        return self.__c

    @c.setter
    def c(self, double c):
        self.__c = c

    cpdef double x_point(self, double t):
        return self.__origin[0] + self.__a * t

    cpdef double y_point(self, double t):
        return self.__origin[1] + self.__b * t

    cpdef double z_point(self, double t):
        return self.__origin[2] + self.__c * t

    cpdef double tangent_x_point(self, double t, bint left=True, bint right=True):
        return self.__a

    cpdef double tangent_y_point(self, double t, bint left=True, bint right=True):
        return self.__b

    cpdef double tangent_z_point(self, double t, bint left=True, bint right=True):
        return self.__c


cdef class Arc(ParametricCurve):

    def __init__(self, str name='Arc', Cartesian coordinate_system=None,
                 double a=1.0, double b=1.0,
                 double start=0.0, double stop=M_PI * 2, bint right=True):
        self.__a = max(a, b)
        self.__b = min(a, b)
        if right:
            self.__direction = 1
        else:
            self.__direction = -1
        super(Arc, self).__init__(name=name, coordinate_system=coordinate_system,
                                  start=start, stop=stop)

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, double a):
        self.__a = a

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, double b):
        self.__b = b

    @property
    def direction(self):
        return self.__direction

    @direction.setter
    def direction(self, short direction):
        self.__direction = direction

    @property
    def right(self):
        if self.__direction > 0:
            return True
        return False

    @right.setter
    def right(self, bint right):
        if right:
            self.__direction = 1
        else:
            self.__direction = -1

    @property
    def left(self):
        if self.__direction > 0:
            return False
        return True

    @left.setter
    def left(self, bint left):
        if left:
            self.__direction = -1
        else:
            self.__direction = 1

    cpdef double x_point(self, double t):
        return self.__a * cos(t)

    cpdef double y_point(self, double t):
        return self.__direction * self.__b * sin(t)

    cpdef double z_point(self, double t):
        return 0.0

    cpdef double tangent_x_point(self, double t, bint left=True, bint right=True):
        return -self.__a * sin(t)

    cpdef double tangent_y_point(self, double t, bint left=True, bint right=True):
        return self.__direction * self.__b * cos(t)

    cpdef double tangent_z_point(self, double t, bint left=True, bint right=True):
        return 0.0

    cpdef double eccentricity(self):
        return sqrt((self.__a * self.__a - self.__b * self.__b) / (self.__a * self.__a))

    cpdef double focus(self):
        return self.__a * self.eccentricity()


cdef class Helix(ParametricCurve):

    def __init__(self, str name='Helix', Cartesian coordinate_system=None,
                 double radius=1.0, double pitch=1.0,
                 double start=0.0, double stop=10.0, double right=True):
        self.__radius = radius
        self.__pitch = pitch
        if right:
            self.__direction = 1
        else:
            self.__direction = -1
        super(Helix, self).__init__(name=name, coordinate_system=coordinate_system,
                                    start=start, stop=stop)

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self, double radius):
        self.__radius = radius

    @property
    def pitch(self):
        return self.__pitch

    @pitch.setter
    def pitch(self, double pitch):
        self.__pitch = pitch

    @property
    def direction(self):
        return self.__direction

    @direction.setter
    def direction(self, short direction):
        self.__direction = direction

    @property
    def right(self):
        if self.__direction > 0:
            return True
        return False

    @right.setter
    def right(self, bint right):
        if right:
            self.__direction = 1
        else:
            self.__direction = -1

    @property
    def left(self):
        if self.__direction > 0:
            return False
        return True

    @left.setter
    def left(self, bint left):
        if left:
            self.__direction = -1
        else:
            self.__direction = 1

    cpdef double x_point(self, double t):
        return self.__radius - self.__radius * cos(t)

    cpdef double y_point(self, double t):
        return self.__direction * self.__radius * sin(t)

    cpdef double z_point(self, double t):
        return self.__pitch / (2 * M_PI) * t

    cpdef double tangent_x_point(self, double t, bint left=True, bint right=True):
        return self.__radius * sin(t)

    cpdef double tangent_y_point(self, double t, bint left=True, bint right=True):
        return self.__direction * self.__radius * cos(t)

    cpdef double tangent_z_point(self, double t, bint left=True, bint right=True):
        return self.__pitch / (2 * M_PI)
