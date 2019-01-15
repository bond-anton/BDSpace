from __future__ import division, print_function
import unittest
import numpy as np

from BDSpace.Coordinates._utils import check_points_array as check_points


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
