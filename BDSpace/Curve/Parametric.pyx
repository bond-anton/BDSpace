import numpy as np

from cython import boundscheck, wraparound

from cpython.array cimport array, clone
from libc.math cimport sin, cos, sqrt, M_PI
from libc.float cimport DBL_MIN
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

    cpdef double x_point(self, double t):
        return 0.0

    cpdef double y_point(self, double t):
        return 0.0

    cpdef double z_point(self, double t):
        return 0.0

    cpdef double[:] x(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        for i in range(s):
            result[i] = self.x_point(t[i])
        return result

    cpdef double[:] y(self, double[:] t):
        cdef:
            unsigned int i, s = t.shape[0]
            array[double] result, template = array('d')
        result = clone(template, t.shape[0], zero=False)
        for i in range(s):
            result[i] = self.y_point(t[i])
        return result

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

    cpdef double[:, :] generate_points(self, double[:] t):
        cdef:
            double[:, :] xyz = np.empty((t.shape[0], 3), dtype=np.double)
        xyz[:, 0] = self.x(t)
        xyz[:, 1] = self.y(t)
        xyz[:, 2] = self.z(t)
        return xyz

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] tangent(self, double[:] t, double step=10*DBL_MIN):
        cdef:
            unsigned int i, s = t.shape[0] - 1
            double[:, :] result = np.empty((s + 1, 3), dtype=np.double)
            double step2 = step / 2
        result[0, 0] = (self.x(t[0] + step)) - self.x(t[0]) / step
        result[0, 1] = (self.y(t[0] + step)) - self.y(t[0]) / step
        result[0, 2] = (self.z(t[0] + step)) - self.z(t[0]) / step
        result[s, 0] = (self.x(t[s]) - self.x(t[s] - step)) / step
        result[s, 1] = (self.y(t[s]) - self.y(t[s] - step)) / step
        result[s, 2] = (self.z(t[s]) - self.z(t[s] - step)) / step
        for i in range(1, s):
            result[i, 0] = (self.x(t[i] + step2) - self.x(t[i] - step2)) / step
            result[i, 1] = (self.y(t[i] + step2) - self.y(t[i] - step2)) / step
            result[i, 2] = (self.z(t[i] + step2) - self.z(t[i] - step2)) / step
        return result

    @boundscheck(False)
    @wraparound(False)
    cdef double __length_tangent_array(self, double[:] t, double tangent_step=10*DBL_MIN):
        cdef:
            unsigned int i, num_points = t.shape[0] - 1
            double[:, :] xyz = self.tangent(t, tangent_step)
            array[double] dl, template = array('d')
        dl = clone(template, num_points, zero=False)
        for i in range(num_points):
            dl[i] = sqrt(xyz[i, 0] * xyz[i, 0] + xyz[i, 1] * xyz[i, 1] + xyz[i, 2] * xyz[i, 2])
        return trapz_1d(t[:num_points], dl)

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

    cdef double __length_tangent_mesh(self, Mesh1DUniform mesh, double tangent_step=10*DBL_MIN):
        cdef:
            unsigned int i, num_points = mesh.num  -1
            double[:, :] xyz = self.generate_points(mesh.physical_nodes)
            double[:, :] xyz_t = self.tangent(mesh.physical_nodes, tangent_step)
            double result_t, result_p, result_t_acc = 0.0, result_p_acc = 0.0, dx, dy, dz, dl1, dl2
            array[double] solution, error, template = array('d')
        solution = clone(template, num_points + 1, zero=False)
        error = clone(template, num_points + 1, zero=False)
        solution[0] = 0.0
        error[0] = 0.0
        for i in range(num_points):
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

    cpdef double length(self, double precision=1e-6, bint print_details=False):
        cdef:
            unsigned int i, num_points, iteration = 0
            double[:] dl ,t
        num_points = int(100 * abs(self.__stop - self.__start) / (2 * M_PI)) + 1
        while True:
            iteration += 1
            t = np.linspace(self.__start, self.__stop, num=num_points, endpoint=True, dtype=np.double)
            length_tangent = self.__length_tangent_array(t)
            length_polygonal = self.__length_poly_array(t)
            if length_tangent == length_polygonal == 0.0:
                error = 0.0
            else:
                error = abs(length_tangent - length_polygonal) / max(length_tangent, length_polygonal)
            if error <= abs(precision):
                break
            num_points *= 2
        if print_details:
            print(max(length_tangent, length_polygonal), error, iteration)
        return max(length_tangent, length_polygonal)


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

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] x(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__origin[0] + self.__a * t[i]
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] y(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__origin[0] + self.__b * t[i]
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] z(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__origin[0] + self.__c * t[i]
        return result


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

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] x(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__a * cos(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] y(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__direction * self.__b * sin(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] z(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int s = t.shape[0]
        result = clone(template, s, zero=True)
        return result

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

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] x(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__radius - self.__radius * cos(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] y(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__direction * self.__radius * sin(t[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] z(self, double[:] t):
        cdef:
            array[double] result, template = array('d')
            unsigned int i, s = t.shape[0]
        result = clone(template, s, zero=False)
        for i in range(s):
            result[i] = self.__pitch / (2 * M_PI) * t[i]
        return result
