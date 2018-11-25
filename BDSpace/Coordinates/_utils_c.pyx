cpdef double[:] check_points_array(double[:] xyz):
    """
    Checks if input is a single 3d point or an array of 3d points
    :param xyz: single 3d point or Nx3 array of 3d points
    :return: Nx3 numpy array of 3d points
    """
    cdef:
        Py_ssize_t[:] shape = xyz.shape
        Py_ssize_t n_axes = shape.size
    if n_axes == 2:
        if shape[1] == 3:
            return xyz
        elif shape[0] == 3:
            return xyz.T
        else:
            raise ValueError('Input must be a single point or an array of points coordinates with shape Nx3')
    elif n_axes == 1 and xyz.size == 3:
        return xyz.reshape(1, 3)
    else:
        raise ValueError('Input must be a single point or an array of points coordinates with shape Nx3')
