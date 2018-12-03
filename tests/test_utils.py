from __future__ import division, print_function
import unittest
import numpy as np
import timeit

from BDSpace.Coordinates._utils import check_points_array as check_points
from BDSpace.Coordinates.Cartesian_c import check_points_array as check_points_c


class TestCoordinatesUtils(unittest.TestCase):

    def test_check_points_array(self):
        point = np.arange(3)
        np.testing.assert_allclose(check_points(point), point.reshape(1, 3))
        point = np.arange(36).reshape(12, 3)
        np.testing.assert_allclose(check_points(point), point.reshape(12, 3))
        point = np.arange(36, dtype=np.double).reshape(12, 3).T
        np.testing.assert_allclose(check_points(point), point.T)
        for point in ['a', 1, np.arange(2), np.arange(4), np.arange(36).reshape(6, 6)]:
            with self.assertRaises(ValueError):
                check_points(point)

    def test_check_points_array_c(self):
        point = np.arange(3)
        np.testing.assert_allclose(check_points_c(point), point.reshape(1, 3))
        point = np.arange(36).reshape(12, 3)
        np.testing.assert_allclose(check_points_c(point), point.reshape(12, 3))
        point = np.arange(36, dtype=np.double).reshape(12, 3).T
        np.testing.assert_allclose(check_points_c(point), point.T)
        for point in ['a', 1, np.arange(2), np.arange(4), np.arange(36).reshape(6, 6)]:
            with self.assertRaises(ValueError):
                check_points_c(point)

    def test_speed(self):
        print()
        s1 = timeit.timeit('check_points(np.arange(36, dtype=np.double).reshape(12, 3).T)',
                          setup='import numpy as np\nfrom BDSpace.Coordinates._utils import check_points_array as check_points',
                          number=100000)
        print('Check Array of Points Py:', s1)
        s2 = timeit.timeit('check_points(np.arange(36, dtype=np.double).reshape(12, 3).T)',
                          setup='import numpy as np\nfrom BDSpace.Coordinates.Cartesian_c import check_points_array as check_points',
                          number=100000)
        print('Check Array of Points Cy:', s2)
        print('Check Array of Points Cy speedup: %2.2f%%' % ((s1 - s2) / s1 * 100))
