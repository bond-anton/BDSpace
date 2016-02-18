from __future__ import division
import numbers
import numpy as np
import math as m


def reduce_angle(angle, keep_sign=False):
    """
    Adjusts rotation angle to be in the range [-2*pi; 2*pi]
    :param angle: angle or array-like of input angle
    :param keep_sign: if False (default) adjust angle to be within [0; 2*pi]
    :return: reduced angle or array of reduced angles
    """
    if isinstance(angle, numbers.Number):
        if angle > 2 * m.pi:
            reduced_angle = angle - 2 * m.pi * (angle // (2 * m.pi))
        elif angle < -2 * m.pi:
            reduced_angle = angle + 2 * m.pi * (abs(angle) // (2 * m.pi))
        else:
            reduced_angle = angle
        if not keep_sign and reduced_angle < 0:
            reduced_angle += 2 * m.pi
    elif isinstance(angle, np.ndarray):
        reduced_angle = np.copy(angle)
        big_pos_idx = np.where(reduced_angle > 2 * np.pi)
        big_neg_idx = np.where(reduced_angle < -2 * np.pi)
        reduced_angle[big_pos_idx] = angle[big_pos_idx] - 2 * np.pi * (angle[big_pos_idx] // (2 * np.pi))
        reduced_angle[big_neg_idx] = angle[big_neg_idx] + 2 * np.pi * (abs(angle[big_neg_idx]) // (2 * np.pi))
        if not keep_sign:
            reduced_angle[np.where(reduced_angle < 0)] += 2 * np.pi
    elif isinstance(angle, (tuple, list)):
        reduced_angle = reduce_angle(np.array(angle), keep_sign=keep_sign)
    else:
        raise ValueError('Input angle must be either number or iterable of numbers.')
    return reduced_angle


def unit_vector(v):
    """
    returns unit vector for given vector v
    :param v: input vector
    :return: unit vector u parallel v of length 1
    """
    length = np.sqrt(np.dot(v, v))
    if length == 0:
        raise ValueError('Can not calculate unit vector for null-vector')
    else:
        return v / length


def rotation_matrix(axis, theta):
    """
    Calculates 3D rotation matrix for given axis and angle
    :param axis: rotation axis (3D vector)
    :param theta: rotation angle
    :return: rotation matrix
    """
    axis_e = unit_vector(axis)
    a = m.cos(theta/2)
    b, c, d = -axis_e * m.sin(theta / 2)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])


def rotation_matrix_euler_angles(euler_angles):
    """
    Calculates 3D rotation matrix for cartesian coordinates given 3 Euler's angles.
    Z-X'-Z" notation is assumed.
    :param euler_angles: 3 Euler's angles
    :return: rotation matrix
    """
    euler_angles = reduce_angle(euler_angles)
    c = np.cos(euler_angles)
    s = np.sin(euler_angles)
    return np.array([[c[0]*c[2]-c[1]*s[0]*s[2], -c[0]*s[2]-c[1]*c[2]*s[0],  s[0]*s[1]],
                    [c[2]*s[0]+c[0]*c[1]*s[2],   c[0]*c[1]*c[2]-s[0]*s[2], -c[0]*s[1]],
                    [s[1]*s[2],                  c[2]*s[1],                 c[1]]])


def rotate_vector(vec, axis, theta):
    """
    Rotates vector around given axis by given angle
    :param vec: input vector to rotate
    :param axis: rotation axis
    :param theta: rotation angle
    :return: rotated vector
    """
    return np.dot(rotation_matrix(axis, theta), vec)


def to_polar(xyz):
    """
    convert cartesian to polar coordinates
    :param xyz: Cartesian coordinates (at least 3)
    """
    xyz = np.array(xyz)
    r_theta_phi = np.zeros(xyz.shape)
    if xyz.size == 3:
        xy = xyz[0]**2 + xyz[1]**2
        r_theta_phi[0] = np.sqrt(xy + xyz[2]**2)
        r_theta_phi[1] = np.arctan2(np.sqrt(xy), xyz[2])
        r_theta_phi[2] = reduce_angle(np.arctan2(xyz[1], xyz[0]), keep_sign=False)
    elif xyz.size > 3:
        if len(xyz.shape) == 2 and xyz.shape[1] == 3:
            xy = xyz[:, 0]**2 + xyz[:, 1]**2
            r_theta_phi[:, 0] = np.sqrt(xy + xyz[:, 2]**2)
            r_theta_phi[:, 1] = np.arctan2(np.sqrt(xy), xyz[:, 2])
            r_theta_phi[:, 2] = reduce_angle(np.arctan2(xyz[:, 1], xyz[:, 0]), keep_sign=False)
        else:
            raise ValueError('N-points array shape must be (N, 3)')
    else:
        raise ValueError('at least 3 coordinates are needed for conversion')
    return r_theta_phi


def to_cartesian(r_theta_phi):
    """
    convert polar to cartesian coordinates
    :param r_theta_phi: polar coordinates (at least 3)
    """
    r_theta_phi = np.array(r_theta_phi)
    xyz = np.zeros(r_theta_phi.shape)
    if r_theta_phi.size == 3:
        r_theta_phi[2] = reduce_angle(r_theta_phi[2], keep_sign=False)
        if (r_theta_phi < 0).any() or r_theta_phi[1] > np.pi or r_theta_phi[2] > 2 * np.pi:
            raise ValueError('r must be >= 0, theta must be within [0; pi], phi must be within [0; 2*pi]')
        xy = r_theta_phi[0] * np.sin(r_theta_phi[1])
        xyz[0] = xy * np.cos(r_theta_phi[2])
        xyz[1] = xy * np.sin(r_theta_phi[2])
        xyz[2] = r_theta_phi[0] * np.cos(r_theta_phi[1])
    elif r_theta_phi.size > 3:
        if len(r_theta_phi.shape) == 2 and r_theta_phi.shape[1] == 3:
            r_theta_phi[:, 2] = reduce_angle(r_theta_phi[:, 2], keep_sign=False)
            if (r_theta_phi < 0).any() or (r_theta_phi[:, 1] > np.pi).any() or (r_theta_phi[:, 2] > 2 * np.pi).any():
                print r_theta_phi
                raise ValueError('r must be >= 0, theta must be within [0; pi], phi must be within [0; 2*pi]')
            xy = r_theta_phi[:, 0] * np.sin(r_theta_phi[:, 1])
            xyz[:, 0] = xy * np.cos(r_theta_phi[:, 2])
            xyz[:, 1] = xy * np.sin(r_theta_phi[:, 2])
            xyz[:, 2] = r_theta_phi[:, 0] * np.cos(r_theta_phi[:, 1])
        else:
            raise ValueError('N-points array shape must be (N, 3)')
    else:
        raise ValueError('at least 3 coordinates are needed for conversion')
    return xyz
