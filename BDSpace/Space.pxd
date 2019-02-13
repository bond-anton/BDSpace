from BDSpace.Coordinates.Cartesian cimport Cartesian

cdef class Space(object):
    cdef:
        str __name
        Cartesian __coordinate_system

        Space __parent
        dict __elements

    cpdef bint add_element(self, Space element)
    cpdef bint remove_element(self, Space element)
    cpdef void detach_from_parent(self)
    cpdef void print_tree(self, int level=*)

    cpdef double[:] to_global_coordinate_system_vector(self, double[:] xyz)
    cpdef double[:, :] to_global_coordinate_system(self, double[:, :] xyz)
    cpdef Cartesian basis_in_global_coordinate_system(self)
    cpdef double[:] to_local_coordinate_system_vector(self, double[:] xyz)
    cpdef double[:, :] to_local_coordinate_system(self, double[:, :] xyz)
