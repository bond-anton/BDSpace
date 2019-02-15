import numpy as np

from BDSpace.Figure import Figure


class Parallelepiped(Figure):

    def __init__(self, name='Parallelepiped', coordinate_system=None,
                 vectors=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)):
        self.__vectors = None
        self.vectors = vectors
        super(Parallelepiped, self).__init__(name, coordinate_system=coordinate_system)

    @property
    def vectors(self):
        return self.__vectors

    @vectors.setter
    def vectors(self, vectors):
        if isinstance(vectors, (tuple, list, np.ndarray)):
            v = np.array(vectors, dtype=np.float64)
            if v.shape == (3, 3):
                self.__vectors = v
            else:
                raise ValueError('Needed 3 vectors, received %s' % str(vectors))
        else:
            raise ValueError('Needed 3 vectors, received %s' % str(vectors))

    def inner_volume(self):
        return 0.0

    def external_volume(self):
        return np.dot(self.vectors[2], np.cross(self.vectors[0], self.vectors[1]))

    def inner_surface_area(self):
        return 0.0

    def external_surface_area(self):
        p1 = np.cross(self.vectors[0], self.vectors[1])
        p2 = np.cross(self.vectors[0], self.vectors[2])
        p3 = np.cross(self.vectors[1], self.vectors[2])
        return 2 * (np.sqrt(np.dot(p1, p1)) + np.sqrt(np.dot(p2, p2)) + np.sqrt(np.dot(p3, p3)))


class ParallelepipedTriclinic(Parallelepiped):

    def __init__(self, name='Parallelepiped', coordinate_system=None,
                 a=1.0, b=1.0, c=1.0,
                 alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2):
        self.__a = None
        self.a = a
        self.__b = None
        self.b = b
        self.__c = None
        self.c = c
        self.__alpha = None
        self.alpha = alpha
        self.__beta = None
        self.beta = beta
        self.__gamma = None
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

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, a):
        self.__a = np.float64(a)

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, b):
        self.__b = np.float64(b)

    @property
    def c(self):
        return self.__c

    @c.setter
    def c(self, c):
        self.__c = np.float64(c)

    @property
    def alpha(self):
        return self.__alpha

    @alpha.setter
    def alpha(self, alpha):
        self.__alpha = np.float64(alpha)

    @property
    def beta(self):
        return self.__beta

    @beta.setter
    def beta(self, beta):
        self.__beta = np.float64(beta)

    @property
    def gamma(self):
        return self.__gamma

    @gamma.setter
    def gamma(self, gamma):
        self.__gamma = np.float64(gamma)


class Cuboid(ParallelepipedTriclinic):

    def __init__(self, name='Cuboid', coordinate_system=None, a=1.0, b=1.0, c=1.0):
        super(Cuboid, self).__init__(name, coordinate_system=coordinate_system,
                                     a=a, b=b, c=c,
                                     alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2)


class Cube(Cuboid):

    def __init__(self, name='Cube', coordinate_system=None, a=1.0):
        super(Cube, self).__init__(name, coordinate_system=coordinate_system, a=a, b=a, c=a)
