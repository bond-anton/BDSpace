import unittest
import numpy as np
from Space.Coordinates import Cartesian


class TestCoordinates(unittest.TestCase):

    def setUp(self):
        self.coordinate_system = Cartesian()

    def test_equality(self):
        other_coordinate_system = Cartesian()
        self.assertEqual(self.coordinate_system, other_coordinate_system)
        other_coordinate_system.set_origin([1, 0, 0])
        self.assertNotEqual(self.coordinate_system, other_coordinate_system)

    def test_rotation_axis_angle(self):
        other_coordinate_system = Cartesian()
        order = 3
        axis = [1, 1, 2]
        steps = 10**order  # per turn
        step = 2 * np.pi / steps
        for k in range(steps):
            other_coordinate_system.rotate_axis_angle(axis, step)
        self.assertEqual(self.coordinate_system, other_coordinate_system)
        np.testing.assert_allclose(self.coordinate_system.basis,
                                   other_coordinate_system.basis, atol=np.finfo(float).eps*steps)
        axis = [1, 0, 0]
        self.coordinate_system.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.coordinate_system.basis,
                                   np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), atol=np.finfo(float).eps)
        self.coordinate_system.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.coordinate_system.basis,
                                   np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), atol=2*np.finfo(float).eps)

    def test_euler_angles(self):
        axis = [1, 0, 0]
        self.coordinate_system.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.coordinate_system.euler_angles,
                                   np.array([0, np.pi, 0]), atol=np.finfo(float).eps)
        self.coordinate_system.rotate_axis_angle(axis, np.pi)
        np.testing.assert_allclose(self.coordinate_system.euler_angles,
                                   np.array([0, 0, 0]), atol=np.finfo(float).eps)
        self.coordinate_system.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.coordinate_system.euler_angles,
                                   np.array([0, np.pi * 0.5, 0]), atol=np.finfo(float).eps)
        self.coordinate_system.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.coordinate_system.euler_angles,
                                   np.array([0, np.pi, 0]), atol=np.finfo(float).eps)
        self.coordinate_system.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.coordinate_system.euler_angles,
                                   np.array([np.pi, np.pi * 0.5, np.pi]), atol=np.finfo(float).eps)
        self.coordinate_system.rotate_axis_angle(axis, np.pi * 0.5)
        np.testing.assert_allclose(self.coordinate_system.euler_angles,
                                   np.array([0, 0, 0]), atol=np.finfo(float).eps)

    def test_to_parent_to_local(self):
        origin = (np.random.random(3) - 0.5) * 100
        other_coordinate_system = Cartesian(origin=origin)
        axis = (np.random.random(3) - 0.5) * 100
        angle = (np.random.random() - 0.5) * 100
        other_coordinate_system.rotate_axis_angle(axis, angle)
        point_global = (np.random.random(3) - 0.5) * 100
        point_local = other_coordinate_system.to_local(point_global)
        np.testing.assert_allclose(other_coordinate_system.to_parent(point_local),
                                   point_global, atol=np.finfo(float).eps)
        point_local = (np.random.random(3) - 0.5) * 100
        point_global = other_coordinate_system.to_parent(point_local)
        np.testing.assert_allclose(other_coordinate_system.to_local(point_global),
                                   point_local, atol=np.finfo(float).eps)
