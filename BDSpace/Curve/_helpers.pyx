import numpy as np
from cython import boundscheck, wraparound

from BDMesh.Mesh1DUniform cimport Mesh1DUniform


@boundscheck(False)
@wraparound(False)
cdef double trapz_1d(double[:] y, double[:] x):
    cdef:
        int nx = x.shape[0], i
        double result = 0.0
    for i in range(nx - 1):
        result += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2
    return result

@boundscheck(False)
@wraparound(False)
cdef int refinement_chunks(Mesh1DUniform mesh, double threshold):
    cdef:
        int i, last = -2, n = mesh.num, result = 0
    for i in range(n):
        if mesh.residual[i] > threshold:
            if i - last > 1:
                result += 1
            last = i
    return result

@boundscheck(False)
@wraparound(False)
cdef long[:, :] refinement_points(Mesh1DUniform mesh, double threshold):
    cdef:
        int i, j = 0, n = mesh.num
        long[:, :] result = np.empty((refinement_chunks(mesh, threshold), 2))
    for i in range(n):
        if mesh.residual[i] > threshold:
            if i - last > 1:
                result[j, 0] = i
            last = i
        elif i - last == 1:
            result[j, 1] = last
            j += 1
    return result
