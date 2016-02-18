import unittest
import numpy as np
from Space.Coordinates import Cartesian


class TestCoordinates(unittest.TestCase):

    def setUp(self):
        self.CS = Cartesian()

    def test_equality(self):
        CS1 = Cartesian()
        self.assertEqual(self.CS, CS1)
        CS1.set_origin([1, 0, 0])
        self.assertNotEqual(self.CS, CS1)

    def test_rotation_axis_angle(self):
        CS1 = Cartesian()
        order = 3
        axis = [1, 1, 2]
        steps = 10**order # per turn
        step = 2 * np.pi / steps
        for k in range(steps):
            CS1.rotate_axis_angle(axis, step)
        self.assertEqual(self.CS, CS1)
        np.testing.assert_allclose(self.CS.basis,
                                   CS1.basis, atol=np.finfo(float).eps*steps)
        axis = [1, 0, 0]
        self.CS.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.CS.basis,
                                   np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), atol=np.finfo(float).eps)
        self.CS.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.CS.basis,
                                   np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), atol=2*np.finfo(float).eps)

