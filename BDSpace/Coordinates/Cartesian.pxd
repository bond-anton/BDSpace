from BDQuaternions cimport Rotation

cdef class Cartesian(object):
    cdef:
        Rotation __rotation
        str __name
        double[:] __origin
        list __labels

    cpdef rotate(self, Rotation rotation, double[:] rot_center=*)
    cpdef rotate_axis_angle(self, double[:] axis, double theta, double[:] rot_center=*)
    cpdef rotate_euler_angles(self, double[:] euler_angles, double[:] rot_center=*)
    cpdef double[:, :] to_parent(self, double[:, :] xyz)
    cpdef double[:, :] to_local(self, double[:, :] xyz)
