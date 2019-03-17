from cython import boundscheck, wraparound
from cython.parallel import prange

from cpython.array cimport array, clone

from BDSpace.Coordinates.transforms cimport cartesian_to_spherical_point, cartesian_to_spherical
from BDSpace.Coordinates.transforms cimport spherical_to_cartesian_point, spherical_to_cartesian
from .Field cimport Field


cdef class SuperposedField(Field):

    def __init__(self, str name, list fields):
        self.__fields = []
        self.type = None
        self.fields = fields
        super(SuperposedField, self).__init__(name, self.type)

    @property
    def fields(self):
        return self.__fields

    @fields.setter
    def fields(self, list fields):
        self.__fields = []
        for field in fields:
            if not isinstance(field, Field):
                raise ValueError('Fields must be iterable of Field class instances')
            if self.type is None:
                self.type = field.type
            if self.type != field.type:
                raise ValueError('All fields must be iterable of Field class instances')
            self.__fields.append(field)

    @boundscheck(False)
    @wraparound(False)
    cpdef double scalar_field_point(self, double[:] xyz):
        """
        Calculates scalar field value at point xyz
        :param xyz: array of cartesian coordinates of the point
        :return: scalar field value
        """
        cdef:
            int i, n_fields = len(self.__fields)
            double[:] global_xyz = self.to_global_coordinate_system_vector(xyz)
            double total_field = 0.0
        for i in range(n_fields):
            total_field += self.__fields[i].scalar_field_point(
                self.__fields[i].to_local_coordinate_system_vector(global_xyz)
            )
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field(self, double[:, :] xyz):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
        cdef:
            int i, j, s = xyz.shape[0], n_fields = len(self.__fields)
            double[:, :] global_xyz = self.to_global_coordinate_system(xyz)
            double[:] total_field = self.__points_scalar(global_xyz, 0.0)
            double[:] field_contribution
        for j in range(n_fields):
            field_contribution = self.__fields[j].scalar_field(self.__fields[j].to_local_coordinate_system(global_xyz))
            with nogil:
                for i in prange(s):
                    total_field[i] += field_contribution[i]
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double scalar_field_polar_point(self, double[:] rtp):
        """
        Calculates scalar field value at point rtp in polar coordinates
        :param rtp: array of polar coordinates of the point
        :return: scalar field value
        """
        cdef:
            int i, n_fields = len(self.__fields)
            double[:] global_xyz = self.to_global_coordinate_system_vector(spherical_to_cartesian_point(rtp))
            double[:] local_rtp
            double total_field = 0.0
        for i in range(n_fields):
            local_rtp = cartesian_to_spherical_point(self.__fields[i].to_local_coordinate_system_vector(global_xyz))
            total_field += self.__fields[i].scalar_field_polar_point(local_rtp)
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] scalar_field_polar(self, double[:, :] rtp):
        """
        Calculates scalar field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: scalar values array
        """
        cdef:
            int i, j, s = rtp.shape[0], n_fields = len(self.__fields)
            double[:, :] global_xyz = self.to_global_coordinate_system(spherical_to_cartesian(rtp))
            double[:, :] local_rtp
            double[:] total_field = self.__points_scalar(global_xyz, 0.0)
            double[:] field_contribution
        for j in range(n_fields):
            local_rtp = cartesian_to_spherical(self.__fields[j].to_local_coordinate_system(global_xyz))
            field_contribution = self.__fields[j].scalar_field_polar(local_rtp)
            with nogil:
                for i in prange(s):
                    total_field[i] += field_contribution[i]
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] vector_field_point(self, double[:] xyz):
        """
        Calculates vector field value at point xyz
        :param xyz: array of cartesian coordinates of the point
        :return: vector field value
        """
        cdef:
            int i, n_fields = len(self.__fields)
            array[double] template = array('d')
            double[:] field_contribution
            double[:] total_field = clone(template, 3, zero=True)
            double[:] global_xyz = self.to_global_coordinate_system_vector(xyz)
        for i in range(n_fields):
            field_contribution = self.__fields[i].vector_field_point(
                self.__fields[i].to_local_coordinate_system_vector(global_xyz)
            )
            total_field[0] += field_contribution[0]
            total_field[1] += field_contribution[1]
            total_field[2] += field_contribution[2]
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field(self, double[:, :] xyz):
        """
        Calculates vector field value at points xyz
        :param xyz: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            int i, j, k, s = xyz.shape[0], n_fields = len(self.__fields)
            array[double] template = array('d')
            double[:, :] field_contribution
            double[:, :] total_field = self.__points_vector(xyz, clone(template, 3, zero=True))
            double[:, :] global_xyz = self.to_global_coordinate_system(xyz)
        for k in range(n_fields):
            field_contribution = self.__fields[k].vector_field(self.__fields[k].to_local_coordinate_system(global_xyz))
            with nogil:
                for i in prange(s):
                    for j in prange(3):
                        total_field[i][j] += field_contribution[i][j]
        return total_field

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] vector_field_polar_point(self, double[:] rtp):
        """
        Calculates vector field value at point rtp in polar coordinates
        :param rtp: array of polar coordinates of the point
        :return: vector field value
        """
        cdef:
            int i, n_fields = len(self.__fields)
            array[double] template = array('d')
            double[:] field_contribution
            double[:] total_field = clone(template, 3, zero=True)
            double[:] global_xyz = self.to_global_coordinate_system_vector(spherical_to_cartesian_point(rtp))
        for i in range(n_fields):
            field_contribution = self.__fields[i].vector_field_point(
                self.__fields[i].to_local_coordinate_system_vector(global_xyz)
            )
            total_field[0] += field_contribution[0]
            total_field[1] += field_contribution[1]
            total_field[2] += field_contribution[2]
        return cartesian_to_spherical_point(total_field)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:, :] vector_field_polar(self, double[:, :] rtp):
        """
        Calculates vector field value at points rtp in polar coordinates
        :param rtp: array of N points with shape (N, 3)
        :return: vector field values array
        """
        cdef:
            int i, j, k, s = rtp.shape[0], n_fields = len(self.__fields)
            array[double] template = array('d')
            double[:, :] field_contribution
            double[:, :] global_xyz = self.to_global_coordinate_system(spherical_to_cartesian(rtp))
            double[:, :] total_field = self.__points_vector(rtp, clone(template, 3, zero=True))
        for k in range(n_fields):
            field_contribution = self.__fields[k].vector_field(self.__fields[k].to_local_coordinate_system(global_xyz))
            with nogil:
                for i in prange(s):
                    for j in prange(3):
                        total_field[i][j] += field_contribution[i][j]
        return cartesian_to_spherical(total_field)
