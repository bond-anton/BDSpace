'''
Created on Mar 25, 2014

@author: Anton Bondarenko
'''

import numpy as np
import math as m

from transforms import unit_vector, adjust_euler_angles, rotation_matrix, rotation_matrix_euler_angles

class Cartesian(object):
    '''
    3D cartesian coordinate system
    '''

    
    def __init__(self, basis=None, origin=None, name='Cartesian CS', labels=None):
        # Euler angles Z-X1-Z2 notation
        self.euler_angles = np.array([0, 0, 0], dtype=np.float)
        # Basis in parent CS
        if basis is not None:
            self.set_basis(basis)
        else:
            self.set_basis(np.identity(3, dtype=np.float))
        # Origin in parent CS
        if origin is not None:
            self.origin = origin
        else:
            self.origin = np.array([0, 0, 0], dtype=np.float)
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

    def set_origin(self, origin):
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
                if not np.allclose(np.dot(basis[:, 0], basis[:, 1]), 0):
                    raise ValueError('only orthogonal vectors accepted')
                if not np.allclose(np.cross(basis[:, 0], basis[:, 1]), basis[:, 2]):
                        raise ValueError('only right-hand basis accepted')
                self.basis = basis
            else:
                raise ValueError('complete 3D basis is needed')
            self.calc_eulers()
    
    def calc_eulers(self):
        if np.allclose(self.basis[2, 2], 1.0):
            self.euler_angles[1] = m.acos(1.0)
            self.euler_angles[2] = 0.0
            self.euler_angles[0] = m.atan2(self.basis[1, 0], self.basis[0, 0])
        elif np.allclose(self.basis[2, 2], -1.0):
            self.euler_angles[1] = m.acos(-1.0)
            self.euler_angles[2] = 0.0
            self.euler_angles[0] = m.atan2(self.basis[1, 0], self.basis[0, 0])
        else:
            self.euler_angles[1] = m.acos(self.basis[2, 2])
            self.euler_angles[0] = m.atan2(self.basis[0, 2], -self.basis[1, 2])
            self.euler_angles[2] = m.atan2(self.basis[2, 0], self.basis[2, 1])
        self.euler_angles = adjust_euler_angles(self.euler_angles)
    
    def set_eulers(self, euler_angles):
        self.euler_angles = adjust_euler_angles(euler_angles)
        self.basis = rotation_matrix_euler_angles(euler_angles)
        self.calc_eulers()
        
    def rotate_axis_angle(self, axis, theta, rot_center=None):
        '''
        rotate basis around axis in parent CS
        '''
        if rot_center is None:
            rot_center = self.origin
        R = rotation_matrix(axis, theta)
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(R, shift)
        self.basis = np.dot(R, self.basis)
        self.calc_eulers()
        
    def rotate_eulers_angles(self, euler_angles, rot_center=None):
        '''
        rotate basis around axis in parent CS
        '''
        if rot_center is None:
            rot_center = self.origin
        R = rotation_matrix_euler_angles(euler_angles)
        shift = self.origin - rot_center
        self.origin = rot_center + np.dot(R, shift)
        self.basis = np.dot(R, self.basis)
        self.calc_eulers()
    
    def to_global(self, xyz):
        '''
        calculate coordinate of points in parent CS
        '''
        if xyz.size == 3:
            return np.dot(self.basis, xyz) + self.origin
        else:
            coords = np.zeros_like(xyz)
            for i in range(xyz.shape[0]):
                coords[i, :] = self.to_global(xyz[i, :])
            return coords

    def to_local(self, xyz, CS_xyz=None):
        '''
        calculates local coords for points in CS_xyz
        if CS_xyz is None global CS is assumed
        '''
        #print xyz
        if CS_xyz is not None:
            coords = CS_xyz.to_global(xyz)        
        if xyz.size == 3:
            #print xyz - self.origin
            #print self.basis
            return np.dot(xyz - self.origin, self.basis)
        else:
            coords = np.zeros_like(xyz)
            for i in range(xyz.shape[0]):
                coords[i, :] = self.to_local(xyz[i, :])
            return coords

    def to_polar(self, xyz):
        '''
        convert cartesian to polar coordinates
        '''
        r_theta_phi = np.zeros(xyz.shape)
        if xyz.size == 3:
            xy = xyz[0]**2 + xyz[1]**2
            r_theta_phi[0] = np.sqrt(xy + xyz[2]**2)
            r_theta_phi[1] = np.arctan2(np.sqrt(xy), xyz[2])
            r_theta_phi[2] = np.arctan2(xyz[1], xyz[0])
        else:
            xy = xyz[:, 0]**2 + xyz[:, 1]**2
            r_theta_phi[:, 0] = np.sqrt(xy + xyz[:, 2]**2)
            r_theta_phi[:, 1] = np.arctan2(np.sqrt(xy), xyz[:, 2])
            r_theta_phi[:, 2] = np.arctan2(xyz[:, 1], xyz[:, 0])
        return r_theta_phi
    
    def to_cartesian(self, r_theta_phi):
        '''
        convert polar to cartesian coordinates
        '''
        xyz = np.zeros(r_theta_phi.shape)
        if r_theta_phi.size == 3:
            xy = r_theta_phi[0] * np.sin(r_theta_phi[1])
            xyz[0] = xy * np.cos(r_theta_phi[2])
            xyz[1] = xy * np.sin(r_theta_phi[2])
            xyz[2] = r_theta_phi[0] * np.cos(r_theta_phi[1])
        else:
            xy = r_theta_phi[:, 0] * np.sin(r_theta_phi[:, 1])
            xyz[:, 0] = xy * np.cos(r_theta_phi[:, 2])
            xyz[:, 1] = xy * np.sin(r_theta_phi[:, 2])
            xyz[:, 2] = r_theta_phi[:, 0] * np.cos(r_theta_phi[:, 1])
        return xyz
