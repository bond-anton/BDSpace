from __future__ import division
import numpy as np

from Space.Figure import Figure


class Parallelepiped(Figure):

    def __init__(self, name='Parallelepiped', coordinate_system=None,
                 vectors=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float)):
        self.vectors = vectors
        super(Parallelepiped, self).__init__(name, coordinate_system=coordinate_system)

    def volume(self):
        return np.dot(self.vectors[2], np.cross(self.vectors[0], self.vectors[1]))


class ParallelepipedTriclinic(Parallelepiped):

    def __init__(self, name='Parallelepiped', coordinate_system=None,
                 a=1.0, b=1.0, c=1.0,
                 alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        v = np.sqrt(abs(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 +
                        2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)))
        v *= a * b * c
        orthogonalization_matrix = np.array([[a, b * np.cos(gamma), c * np.cos(beta)],
                                             [0, b * np.sin(gamma),
                                              c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
                                             [0, 0, v / (a * b * np.sin(gamma))]])
        vectors_fractional = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float)
        vectors_cartesian = np.dot(vectors_fractional, orthogonalization_matrix.T)
        super(ParallelepipedTriclinic, self).__init__(name, coordinate_system=coordinate_system,
                                                      vectors=vectors_cartesian)


class Cuboid(ParallelepipedTriclinic):

    def __init__(self, name='Cuboid', coordinate_system=None, a=1.0, b=1.0, c=1.0):
        super(Cuboid, self).__init__(name, coordinate_system=coordinate_system,
                                     a=a, b=b, c=c,
                                     alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2)


class Cube(Cuboid):

    def __init__(self, name='Cube', coordinate_system=None, a=1.0):
        super(Cube, self).__init__(name, coordinate_system=coordinate_system, a=a, b=a, c=a)
