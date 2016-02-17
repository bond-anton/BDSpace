import unittest
import numpy as np
from Space.Coordinates.transforms import reduce_angle


class TestTransforms(unittest.TestCase):

    def test_reduce_angles_constant(self):
        angle = 0.1
        self.assertEqual(reduce_angle(angle), angle)

    def test_reduce_angles_array(self):
        angle = np.array([0.1, ])
        self.assertEqual(reduce_angle(angle), angle)

    def test_reduce_angles_list(self):
        angle = [0.1, ]
        self.assertEqual(reduce_angle(angle), angle)

    def test_reduce_angles_tuple(self):
        angle = (0.1, )
        self.assertEqual(reduce_angle(angle), angle)

    def test_reduce_angles_negative(self):
        angle = -0.1
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle + 2 * np.pi)
