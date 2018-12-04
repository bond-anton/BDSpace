from BDSpace.Coordinates.Cartesian cimport Cartesian

cdef class Space(object):
    cdef:
        str __name
        Cartesian __coordinate_system

        Space __parent
        dict __elements
