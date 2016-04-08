from __future__ import division, print_function
import numpy as np

from Quaternions import Rotation
from _transforms import unit_vector


class Cartesian(object):
    """
    3D cartesian coordinate system
    """

    def __init__(self, basis=None, origin=None, name='Cartesian CS', labels=None,
                 euler_angles_convention=None):
        # The basis rotation is kept as Rotation quaternion
        self._rotation = Rotation(euler_angles_convention=euler_angles_convention)
        # Basis in parent CS
        self._set_basis(basis)
        # Origin in parent CS
        self._origin = None
        self._set_origin(origin)
        self.labels = None
        if labels is not None:
            self.labels = labels
        else:
            self.labels = ['x', 'y', 'z']
        self.name = str(name)

    def _set_euler_angles_convention(self, euler_angles_convention):
        self._rotation.euler_angles_convention = euler_angles_convention

    def _get_euler_angles_convention(self):
        return self._rotation.euler_angles_convention

    euler_angles_convention = property(_get_euler_angles_convention, _set_euler_angles_convention)

    def _set_euler_angles(self, euler_angles):
        self._rotation.euler_angles = euler_angles

    def _get_euler_angles(self):
        return self._rotation.euler_angles

    euler_angles = property(_get_euler_angles, _set_euler_angles)

    def _set_basis(self, basis):
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
            basis = np.identity(3, dtype=np.float)
        if isinstance(basis, np.ndarray):
            try:
                basis = basis.astype(np.float)
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
                self._rotation.rotation_matrix = basis.T
            else:
                raise ValueError('complete 3D basis is needed')

    def _get_basis(self):
        return self._rotation.rotation_matrix.T

    basis = property(_get_basis, _set_basis)

    def _set_origin(self, origin):
        if origin is None:
            origin = [0, 0, 0]
        origin = np.array(origin, dtype=np.float).ravel()
        if origin.size != 3:
            raise ValueError('Origin must be 3 numeric coordinates')
        self._origin = origin

    def _get_origin(self):
        return self._origin

    origin = property(_get_origin, _set_origin)

    def __eq__(self, other):
        result = isinstance(other, self.__class__)
        result = result and np.allclose(self.basis, other.basis)
        result = result and np.allclose(self.origin, other.origin)
        return result

    def __str__(self):
        information = 'Cartesian coordinate system: %s\n' % self.name
        information += 'Origin: ' + str(self.origin) + '\n'
        information += 'Basis:\n'
        information += self.labels[0] + ': ' + str(self.basis[0]) + '\n'
        information += self.labels[1] + ': ' + str(self.basis[1]) + '\n'
        information += self.labels[2] + ': ' + str(self.basis[2]) + '\n'
        information += 'Orientation: %s:\n' % self.euler_angles_convention['description']
        information += str(self.euler_angles) + '\n'
        return information

    def rotate_axis_angle(self, axis, theta, rot_center=None):
        """
        rotate basis around axis in parent CS
        :param axis: axis of rotation
        :param theta: angle of rotation
        :param rot_center: center of rotation, if None the origin of the CS is used
        """
        if rot_center is None:
            rot_center = self.origin
        rotation = Rotation(euler_angles_convention=self.euler_angles_convention)
        rotation.axis_angle = (axis, theta)
        rot_matrix = rotation.rotation_matrix
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(rot_matrix, shift)
        self.basis = np.dot(rot_matrix, self.basis)
        
    def rotate_euler_angles(self, euler_angles, rot_center=None):
        """
        rotate basis in parent CS using three Euler's angles and selected center of rotation
        :param euler_angles: Euler's angles
        :param rot_center: rotation center, if None the origin of the CS is used
        """
        if rot_center is None:
            rot_center = self.origin
        rotation = Rotation(euler_angles_convention=self.euler_angles_convention)
        rotation.euler_angles = euler_angles
        rot_matrix = rotation.rotation_matrix
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(rot_matrix, shift)
        self.basis = np.dot(rot_matrix, self.basis)
    
    def to_parent(self, xyz):
        """
        calculates coordinates of given points in parent (global) CS
        :param xyz: local coordinates array
        """
        xyz = np.array(xyz, dtype=np.float)
        if len(xyz.shape) == 2:
            if xyz.shape[1] != 3 and xyz.shape[0] == 3:
                xyz = xyz.T
            elif 3 not in xyz.shape:
                raise ValueError('Input must be a single point or an array of points coordinates with shape Nx3')
        elif len(xyz.shape) == 1 and xyz.size == 3:
            xyz = xyz.reshape(1, 3)
        else:
            raise ValueError('Input must be a single point or an array of points coordinates with shape Nx3')
        if xyz.size > 3:
            return np.dot(xyz, self.basis) + self.origin
        else:
            return (np.dot(xyz, self.basis) + self.origin)[0]

    def to_local(self, xyz):
        """
        calculates local coordinates for points in parent CS/
        :param xyz: coordinates in parent (global) coordinate system.
        """
        xyz = np.array(xyz, dtype=np.float)
        if len(xyz.shape) == 2:
            if xyz.shape[1] != 3 and xyz.shape[0] == 3:
                xyz = xyz.T
            elif 3 not in xyz.shape:
                raise ValueError('Input must be a single point or an array of points coordinates with shape Nx3')
        elif len(xyz.shape) == 1 and xyz.size == 3:
            xyz = xyz.reshape(1, 3)
        else:
            raise ValueError('Input must be a single point or an array of points coordinates with shape Nx3')
        if xyz.size > 3:
            return np.dot(xyz - self.origin, self.basis.T)
        else:
            return np.dot(xyz - self.origin, self.basis.T)[0]
