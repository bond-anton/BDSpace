import unittest
import numpy as np
from Space.Coordinates.transforms import reduce_angle


class TestTransforms(unittest.TestCase):

    def test_reduce_angle_constant_positive(self):
        angle = 0.1
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle)

    def test_reduce_angle_constant_negative(self):
        angle = -0.1
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle + 2 * np.pi)

    def test_reduce_angle_constant_zero(self):
        angle = 0.0
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle)

    def test_reduce_angle_constant_2pi(self):
        angle = -2 * np.pi
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle + 2 * np.pi)
        angle = 2 * np.pi
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle)

    def test_reduce_angle_array(self):
        angle = np.array([0.1, ])
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=True), angle)
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=False), angle)
        angle = np.array([0.1, -0.2])
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=True), angle)
        angle = np.array([0.1, 0.2])
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=True), angle)
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=False), angle)

    def test_reduce_angle_list(self):
        angle = [0.1, 0.2]
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=True), angle)
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=False), angle)

    def test_reduce_angle_tuple(self):
        angle = (0.1, 0.2)
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=True), angle)
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=False), angle)


    def test_reduce_angle_random_array_positive(self):
        size = 1000
        angle = np.random.random(size) * 2 * np.pi
        full_turns = np.random.random_integers(0, 100, size=size) * 2 * np.pi * 1
        np.testing.assert_allclose(reduce_angle(angle + full_turns, keep_sign=False),
                                   reduce_angle(angle, keep_sign=False))
        np.testing.assert_allclose(reduce_angle(angle + full_turns, keep_sign=True),
                                   reduce_angle(angle, keep_sign=True))
        np.testing.assert_allclose(reduce_angle(angle + full_turns, keep_sign=True), angle)
        np.testing.assert_allclose(reduce_angle(angle + full_turns, keep_sign=False), angle)

    def test_reduce_angle_random_array_loops(self):
        size = 1000
        angle = (np.random.random(size) - 0.5) * 4 * np.pi
        full_turns = np.random.random_integers(-100, 100, size=size) * 2 * np.pi * 1
        np.testing.assert_allclose(reduce_angle(angle + full_turns, keep_sign=False),
                                   reduce_angle(angle, keep_sign=False))

    def test_reduce_angle_random_array(self):
        size = 1000
        angle = (np.random.random(size) - 0.5) * 4 * np.pi
        np.testing.assert_allclose(reduce_angle(angle, keep_sign=True), angle)

