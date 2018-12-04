from .Field cimport Field

cdef class SuperposedField(Field):
    cdef:
        list __fields
