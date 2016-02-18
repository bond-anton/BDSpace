import numpy as np
import math as m

from transforms import unit_vector, reduce_angle, rotation_matrix, rotation_matrix_euler_angles


class Cartesian(object):
    """
    3D cartesian coordinate system
    """

    def __init__(self, basis=None, origin=None, name='Cartesian CS', labels=None):
        # Euler angles Z-X'-Z" notation
        self.euler_angles = np.array([0, 0, 0], dtype=np.float)
        # Basis in parent CS
        self.basis = None
        if basis is not None:
            self.set_basis(basis)
        else:
            self.set_basis(np.identity(3, dtype=np.float))
        # Origin in parent CS
        self.origin = None
        self.set_origin(origin)
        self.labels = None
        if labels is not None:
            self.labels = labels
        else:
            self.labels = ['x', 'y', 'z']
        self.name = str(name)

    def __eq__(self, other):
        result = isinstance(other, self.__class__)
        result = result and np.allclose(self.basis, other.basis)
        result = result and np.allclose(self.origin, other.origin)
        return result

    def __ne__(self, other):
        return not self.__eq__(other)

    def set_origin(self, origin=None):
        if origin is None:
            origin = [0, 0, 0]
        origin = np.array(origin).astype(np.float).ravel()
        if origin.size != 3:
            raise ValueError('Origin must be 3 numeric coordinates')
        self.origin = origin
    
    def set_basis(self, basis):
        if isinstance(basis, np.ndarray):
            try:
                basis = basis.astype(np.float)
            except:
                raise ValueError('only numeric basis coordinates are supported')
            if basis.shape == (3, 3) or basis.size == 9:
                basis = basis.reshape((3, 3))
                for i in range(3):
                    basis[:, i] = unit_vector(basis[:, i])
                if not np.allclose(np.dot(basis[:, 0], basis[:, 1]), [0]):
                    raise ValueError('only orthogonal vectors accepted')
                if not np.allclose(np.cross(basis[:, 0], basis[:, 1]), basis[:, 2]):
                        raise ValueError('only right-hand basis accepted')
                self.basis = basis
            else:
                raise ValueError('complete 3D basis is needed')
            self.calc_euler_angles()
    
    def calc_euler_angles(self):
        if np.allclose(self.basis[2, 2], [1.0]):
            self.euler_angles[1] = m.acos(1.0)
            self.euler_angles[2] = 0.0
            self.euler_angles[0] = m.atan2(self.basis[1, 0], self.basis[0, 0])
        elif np.allclose(self.basis[2, 2], [-1.0]):
            self.euler_angles[1] = m.acos(-1.0)
            self.euler_angles[2] = 0.0
            self.euler_angles[0] = m.atan2(self.basis[1, 0], self.basis[0, 0])
        else:
            self.euler_angles[1] = m.acos(self.basis[2, 2])
            self.euler_angles[0] = m.atan2(self.basis[0, 2], -self.basis[1, 2])
            self.euler_angles[2] = m.atan2(self.basis[2, 0], self.basis[2, 1])
        self.euler_angles = reduce_angle(self.euler_angles)
    
    def set_euler_angles(self, euler_angles):
        self.euler_angles = reduce_angle(euler_angles)
        self.basis = rotation_matrix_euler_angles(euler_angles)
        self.calc_euler_angles()
        
    def rotate_axis_angle(self, axis, theta, rot_center=None):
        """
        rotate basis around axis in parent CS
        :param axis:
        :param theta:
        :param rot_center:
        """
        if rot_center is None:
            rot_center = self.origin
        R = rotation_matrix(axis, theta)
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(R, shift)
        self.basis = np.dot(R, self.basis)
        self.calc_euler_angles()
        
    def rotate_euler_angles(self, euler_angles, rot_center=None):
        """
        rotate basis around axis in parent CS
        :param euler_angles:
        :param rot_center:
        """
        if rot_center is None:
            rot_center = self.origin
        R = rotation_matrix_euler_angles(euler_angles)
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(R, shift)
        self.basis = np.dot(R, self.basis)
        self.calc_euler_angles()
    
    def to_global(self, xyz):
        """
        calculate coordinate of points in parent CS
        :param xyz:
        """
        if xyz.size == 3:
            return np.dot(self.basis, xyz) + self.origin
        else:
            coords = np.zeros_like(xyz)
            for i in range(xyz.shape[0]):
                coords[i, :] = self.to_global(xyz[i, :])
            return coords

    def to_local(self, xyz, xyz_coordinate_system=None):
        """
        calculates local coords for points in CS_xyz
        if CS_xyz is None global CS is assumed
        :param xyz:
        :param xyz_coordinate_system:
        """
        if xyz_coordinate_system is not None:
            coordinates = xyz_coordinate_system.to_global(xyz)
        if xyz.size == 3:
            return np.dot(xyz - self.origin, self.basis)
        else:
            coordinates = np.zeros_like(xyz)
            for i in range(xyz.shape[0]):
                coordinates[i, :] = self.to_local(xyz[i, :])
            return coordinates
