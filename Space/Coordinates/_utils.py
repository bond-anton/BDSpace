from __future__ import division
import numpy as np


def check_points_array(xyz):
    """
    Checks if input is a single 3d point or an array of 3d points
    :param xyz: single 3d point or Nx3 array of 3d points
    :return: Nx3 numpy array of 3d points
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
    return xyz
