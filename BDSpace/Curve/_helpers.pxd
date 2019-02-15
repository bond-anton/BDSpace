from BDMesh.Mesh1DUniform cimport Mesh1DUniform

cdef double trapz_1d(double[:] y, double[:] x) nogil
cdef int refinement_chunks(Mesh1DUniform mesh, double threshold)
cdef long[:, :] refinement_points(Mesh1DUniform mesh, double threshold)
