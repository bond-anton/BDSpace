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
        self.set_basis(basis)
        # Origin in parent CS
        self.origin = None
        self.set_origin(origin)
        self.labels = None
        if labels is not None:
            self.labels = labels
        else:
            self.labels = ['x', 'y', 'z']
        self.name = str(name)
        self.calculate_euler_angles()

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
                    basis[:, i] = unit_vector(basis[:, i])
                if not np.allclose(np.dot(basis[:, 0], basis[:, 1]), [0]):
                    raise ValueError('only orthogonal vectors accepted')
                if not np.allclose(np.cross(basis[:, 0], basis[:, 1]), basis[:, 2]):
                        raise ValueError('only right-hand basis accepted')
                self.basis = basis
            else:
                raise ValueError('complete 3D basis is needed')
            self.calculate_euler_angles()
    
    def calculate_euler_angles(self):
        """
        method calculates Euler's angles for the coordinate system in Z-X'-Z" notation
        The limits for the angles are [0; 2*pi] for Z and Z", and [0: pi] for X'
        """
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
        """
        Sets coordinate system orientation using given 3 Euler's angles in Z-X'-Z" notation
        The input angles could be of any value. After the rotation the correct Euler's angles will be calculated.
        :param euler_angles: input Euler's angles
        """
        self.euler_angles = reduce_angle(euler_angles)
        self.basis = rotation_matrix_euler_angles(euler_angles)
        self.calculate_euler_angles()
        
    def rotate_axis_angle(self, axis, theta, rot_center=None):
        """
        rotate basis around axis in parent CS
        :param axis: axis of rotation
        :param theta: angle of rotation
        :param rot_center: center of rotation, if None the origin of the CS is used
        """
        if rot_center is None:
            rot_center = self.origin
        rot_matrix = rotation_matrix(axis, theta)
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(rot_matrix, shift)
        self.basis = np.dot(rot_matrix, self.basis)
        self.calculate_euler_angles()
        
    def rotate_euler_angles(self, euler_angles, rot_center=None):
        """
        rotate basis in parent CS using three Euler's angles and selected center of rotation
        :param euler_angles: Euler's angles
        :param rot_center: rotation center, if None the origin of the CS is used
        """
        if rot_center is None:
            rot_center = self.origin
        rot_matrix = rotation_matrix_euler_angles(euler_angles)
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(rot_matrix, shift)
        self.basis = np.dot(rot_matrix, self.basis)
        self.calculate_euler_angles()
    
    def to_global(self, xyz):
        """
        calculate coordinate of points in parent CS
        :param xyz: local coordinates array
        """
        if xyz.size == 3:
            return np.dot(self.basis, xyz) + self.origin
        else:
            coordinates = np.zeros_like(xyz)
            for i in range(xyz.shape[0]):
                coordinates[i, :] = self.to_global(xyz[i, :])
            return coordinates

    def to_local(self, xyz, xyz_coordinate_system=None):
        """
        calculates local coordinates for points in CS_xyz
        if CS_xyz is None global CS is assumed
        :param xyz: coordinates in external coordinate system if specified or in global coordinate system.
        :param xyz_coordinate_system: external coordinate system, if None global coordinate system is assumed.
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
