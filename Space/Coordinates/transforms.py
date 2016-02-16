from __future__ import division
import numpy as np
import math as m


def unit_vector(v):
    return v / np.sqrt(np.dot(v, v))


def rotation_matrix(axis, theta):
    axis_e = unit_vector(axis)
    a = m.cos(theta/2)
    b, c, d = -axis_e * m.sin(theta / 2)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])


def rotation_matrix_euler_angles(euler_angles):
    euler_angles = adjust_rotation_angles(euler_angles)
    c = np.cos(euler_angles)
    s = np.sin(euler_angles)
    return np.array([[c[0]*c[2]-c[1]*s[0]*s[2], -c[0]*s[2]-c[1]*c[2]*s[0],  s[0]*s[1]],
                    [c[2]*s[0]+c[0]*c[1]*s[2],   c[0]*c[1]*c[2]-s[0]*s[2], -c[0]*s[1]],
                    [s[1]*s[2],                  c[2]*s[1],                 c[1]]])


def adjust_rotation_angles(angles, keep_sign=False):
    """
    Adjusts rotation angles to be in the range [-2*pi; 2*pi]
    :param angles: array of input angles
    :param keep_sign: if False (default) adjust angle to be within [0; 2*pi]
    :return: adjusted array of angles
    """
    adjusted_angles = []
    for angle in angles:
        if angle > 2 * m.pi:
            angle -= 2*m.pi*(angle // (2*m.pi))
        elif angle < -2 * m.pi:
            angle += 2*m.pi*(abs(angle) // (2*m.pi))
        if not keep_sign and angle < 0:
            angle += 2*m.pi
        adjusted_angles.append(angle)
    return adjusted_angles


def rotate_vector(vec, axis, theta):
    return np.dot(rotation_matrix(axis, theta), vec)
