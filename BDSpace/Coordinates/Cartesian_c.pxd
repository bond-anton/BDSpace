from BDQuaternions cimport Rotation

cdef class Cartesian(object):
    cdef:
        Rotation __rotation
        str __title

        int[4] __euler_next_axis
        str __name
        double[:] __origin
        list __labels

    cpdef rotate(self, Rotation rotation, double[:] rot_center=*)
