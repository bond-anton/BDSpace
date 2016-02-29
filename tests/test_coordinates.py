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

    def test_euler_angles(self):
        axis = [1, 0, 0]
        self.CS.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.CS.euler_angles,
                                   np.array([0, np.pi, 0]), atol=np.finfo(float).eps)
        self.CS.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.CS.euler_angles,
                                   np.array([0, 0, 0]), atol=np.finfo(float).eps)
        self.CS.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.CS.euler_angles,
                                   np.array([0, np.pi * 0.5, 0]), atol=np.finfo(float).eps)
        self.CS.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.CS.euler_angles,
                                   np.array([0, np.pi, 0]), atol=np.finfo(float).eps)
        self.CS.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.CS.euler_angles,
                                   np.array([np.pi, np.pi * 0.5, np.pi]), atol=np.finfo(float).eps)
        self.CS.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.CS.euler_angles,
                                   np.array([0, 0, 0]), atol=np.finfo(float).eps)

    def test_to_parent_to_local(self):
        origin = (np.random.random(3) - 0.5) * 100
        CS1 = Cartesian(origin=origin)
        axis = (np.random.random(3) - 0.5) * 100
        angle = (np.random.random() - 0.5) * 100
        CS1.rotate_axis_angle(axis, angle)
        point_global = ((np.random.random(3) - 0.5) * 100).reshape(1, 3)
        point_local = CS1.to_local(point_global)
        np.testing.assert_allclose(CS1.to_parent(point_local),
                                   point_global, atol=np.finfo(float).eps)
        point_local = ((np.random.random(3) - 0.5) * 100).reshape(1, 3)
        point_global = CS1.to_parent(point_local)
        np.testing.assert_allclose(CS1.to_local(point_global),
                                   point_local, atol=np.finfo(float).eps)
