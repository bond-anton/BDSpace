from __future__ import division, print_function
import numpy as np

from cython import boundscheck, wraparound

from cpython.array cimport array, clone
from cpython.object cimport Py_EQ, Py_NE

from BDQuaternions cimport Rotation, EulerAngles

from .transforms cimport unit_vector
from ._utils cimport check_points_array


cdef class Cartesian(object):
    """
    3D cartesian coordinate system
    """

    def __init__(self, basis=None, origin=None, name='Cartesian CS', labels=None,
                 euler_angles_convention=None):
        # The basis rotation is kept as Rotation quaternion
        self.__rotation = Rotation()
        self.euler_angles_convention = euler_angles_convention
        self.__name = str(name)
        self.labels = labels
        # Basis in parent CS
        self.basis = basis
        # Origin in parent CS
        self.origin = origin

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = str(name)

    @property
    def labels(self):
        return self.__labels

    @labels.setter
    def labels(self, labels):
        if labels is None:
            self.__labels = self.__rotation.euler_angles_convention.axes_labels
        elif isinstance(labels, (list, tuple, np.ndarray)):
            if len(labels) == 3:
                self.__labels = [str(labels[i]) for i in range(3)]
            else:
                raise ValueError('Labels must be iterable of size 3')
        else:
            raise ValueError('Labels must be iterable of size 3')

    @property
    def euler_angles_convention(self):
        return self.__rotation.euler_angles_convention

    @euler_angles_convention.setter
    def euler_angles_convention(self, euler_angles_convention):
        self.__rotation.euler_angles_convention = euler_angles_convention

    @property
    def euler_angles(self):
        return self.__rotation.euler_angles

    @euler_angles.setter
    def euler_angles(self, euler_angles):
        self.__rotation.euler_angles = euler_angles

    @property
    def basis(self):
        return self.__rotation.rotation_matrix.T

    @basis.setter
    @boundscheck(False)
    @wraparound(False)
    def basis(self, basis):
        """
        Sets basis for the cartesian coordinate system. Basis will be converted to 3x3 array.
        Each basis vector will be normed and placed in separate row like this:
          x  y  z
        1 x1 y1 z1
        2 x2 y2 z2
        3 x3 y3 z3
        :param basis: 3x3 numpy array.
        """
        if basis is None:
            basis = np.identity(3, dtype=np.double)
        if isinstance(basis, np.ndarray):
            try:
                basis = basis.astype(np.double)
            except:
                raise ValueError('only numeric basis coordinates are supported')
            if basis.shape == (3, 3) or basis.size == 9:
                basis = basis.reshape((3, 3))
                for i in range(3):
                    basis[i] = unit_vector(basis[i])
                if not np.allclose(np.dot(basis[0], basis[1]), [0]):
                    raise ValueError('only orthogonal vectors accepted')
                if not np.allclose(np.cross(basis[0], basis[1]), basis[2]):
                    raise ValueError('only right-hand basis accepted')
                self.__rotation.rotation_matrix = basis.T
            else:
                raise ValueError('complete 3D basis is needed')

    @property
    def origin(self):
        return self.__origin

    @origin.setter
    def origin(self, origin):
        if origin is None:
            self.__origin = np.zeros(3, dtype=np.double)
        else:
            origin = np.array(origin, dtype=np.double).ravel()
            if origin.size != 3:
                raise ValueError('Origin must be 3 numeric coordinates')
            self.__origin = origin

    def __richcmp__(x, y, int op):
        if op == Py_EQ:
            if isinstance(x, Cartesian) and isinstance(y, Cartesian):
                return np.allclose(x.basis, y.basis) and np.allclose(x.origin, y.origin)
            return False
        elif op == Py_NE:
            if isinstance(x, Cartesian) and isinstance(y, Cartesian):
                return not np.allclose(x.basis, y.basis) or not np.allclose(x.origin, y.origin)
            return True
        else:
            return NotImplemented

    def __str__(self):
        information = 'Cartesian coordinate system: %s\n' % self.name
        information += 'Origin: ' + str(np.asarray(self.origin)) + '\n'
        information += 'Basis:\n'
        information += self.labels[0] + ': ' + str(np.asarray(self.basis[0])) + '\n'
        information += self.labels[1] + ': ' + str(np.asarray(self.basis[1])) + '\n'
        information += self.labels[2] + ': ' + str(np.asarray(self.basis[2])) + '\n'
        information += 'Orientation: %s:\n' % self.euler_angles_convention.description
        information += str(self.euler_angles) + '\n'
        return information

    @boundscheck(False)
    @wraparound(False)
    cpdef rotate(self, Rotation rotation, double[:] rot_center=None):
        """
        Apply rotation specified by Rotation quaternion instance
        :param rotation: Rotation quaternion instance
        :param rot_center: center of rotation, if None the origin of the CS is used
        """
        cdef:
            int i
            Py_ssize_t s = self.__origin.size
            array[double] origin_shift, template = array('d')
        origin_shift = clone(template, s, zero=False)
        self.__rotation *= rotation
        if rot_center is not None:
            for i in range(s):
                origin_shift[i] = self.__origin[i] - rot_center[i]
            self.__origin = rot_center + rotation.rotate(origin_shift)[0]

    cpdef rotate_axis_angle(self, double[:] axis, double theta, double[:] rot_center=None):
        """
        rotate basis around axis in parent CS
        :param axis: axis of rotation
        :param theta: angle of rotation
        :param rot_center: center of rotation, if None the origin of the CS is used
        """
        rotation = Rotation()
        rotation.euler_angles_convention = self.__rotation.euler_angles_convention
        rotation.axis_angle = (axis, theta)
        self.rotate(rotation, rot_center=rot_center)
        
    cpdef rotate_euler_angles(self, double[:] euler_angles, double[:] rot_center=None):
        """
        rotate basis in parent CS using three Euler's angles and selected center of rotation
        :param euler_angles: Euler's angles
        :param rot_center: rotation center, if None the origin of the CS is used
        """
        rotation = Rotation()
        rotation.euler_angles_convention = self.__rotation.euler_angles_convention
        rotation.euler_angles = EulerAngles(euler_angles, rotation.euler_angles_convention)
        self.rotate(rotation, rot_center=rot_center)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] to_parent(self, xyz):
        """
        calculates coordinates of given points in parent (global) CS
        :param xyz: local coordinates array
        """
        cdef:
            int i
            double[:, :] xyz_parent = self.__rotation.rotate(xyz)
            Py_ssize_t[:] s = xyz_parent.shape
        for i in range(0, s[0]):
            xyz_parent[i, 0] += self.__origin[0]
            xyz_parent[i, 1] += self.__origin[1]
            xyz_parent[i, 2] += self.__origin[2]
        return xyz_parent

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] to_local(self, xyz):
        """
        calculates local coordinates for points in parent CS/
        :param xyz: coordinates in parent (global) coordinate system.
        """
        cdef:
            int i
            double[:, :] xyz_local = check_points_array(xyz)
            Py_ssize_t[:] s = xyz_local.shape
        for i in range(0, s[0]):
            xyz_local[i, 0] -= self.__origin[0]
            xyz_local[i, 1] -= self.__origin[1]
            xyz_local[i, 2] -= self.__origin[2]
        xyz_local = self.__rotation.reciprocal().rotate(xyz_local)
        return xyz_local
