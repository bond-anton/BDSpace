import unittest
import numpy as np
from Space.Coordinates.transforms import reduce_angle, unit_vector
from Space.Coordinates.transforms import rotation_matrix, rotation_matrix_euler_angles, rotate_vector
from Space.Coordinates.transforms import to_polar, to_cartesian


class TestTransforms(unittest.TestCase):

    def test_reduce_angle_constant_positive(self):
        angle = 0.1
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle)
        angle = 1
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle)

    def test_reduce_angle_constant_negative(self):
        angle = -0.1
        self.assertEqual(reduce_angle(angle, keep_sign=True), angle)
        self.assertEqual(reduce_angle(angle, keep_sign=False), angle + 2 * np.pi)
        angle = -1
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
        angle = np.array([1, ])
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

    def test_unit_vector_1d(self):
        v = 0.1
        self.assertEqual(unit_vector(v), 1.0)
        v = -0.1
        self.assertEqual(unit_vector(v), -1.0)
        v = [0.1]
        self.assertEqual(unit_vector(v), [1.0])
        v = [-0.1]
        self.assertEqual(unit_vector(v), [-1.0])
        v = np.array([0.1])
        self.assertEqual(unit_vector(v), [1.0])

    def test_unit_vector_null_vector(self):
        self.assertRaises(ValueError, unit_vector, 0.0)
        max_dimensions = 100
        for dim in range(max_dimensions):
            v = np.zeros(dim)
            self.assertRaises(ValueError, unit_vector, v)

    def test_unit_vector_random_vector(self):
        max_dimensions = 100
        for dim in range(max_dimensions):
            v = np.random.random(dim+1) * 100
            np.testing.assert_allclose(unit_vector(v), v / np.sqrt(np.dot(v, v)))

    def test_rotation_matrix_zero_angle(self):
        axis = (np.random.random(3) - 0.5) * 100
        angle = 0
        np.testing.assert_allclose(rotation_matrix(axis, angle), np.eye(3, 3))

    def test_rotation_matrix_two_pi(self):
        axis = (np.random.random(3) - 0.5) * 100
        angle = 2 * np.pi
        np.testing.assert_allclose(rotation_matrix(axis, angle), np.eye(3, 3), atol=2*np.finfo(float).eps)

    def test_rotation_matrix_half_pi(self):
        rx = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
        ry = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
        rz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        np.testing.assert_allclose(rotation_matrix([1, 0, 0], -np.pi/2), rx, atol=np.finfo(float).eps)
        np.testing.assert_allclose(rotation_matrix([0, 1, 0], -np.pi/2), ry, atol=np.finfo(float).eps)
        np.testing.assert_allclose(rotation_matrix([0, 0, 1], -np.pi/2), rz, atol=np.finfo(float).eps)

    def test_rotation_matrix_null_vector(self):
        axis = np.zeros(3)
        angle = np.random.random_sample()
        self.assertRaises(ValueError, rotation_matrix, axis, angle)

    def test_rotattion_matrix_euler_angles_zero(self):
        np.testing.assert_allclose(rotation_matrix_euler_angles([0, 0, 0]), np.eye(3, 3), atol=np.finfo(float).eps)

    def test_rotattion_matrix_euler_angles_two_pi(self):
        np.testing.assert_allclose(rotation_matrix_euler_angles(2 * np.pi * np.ones(3)),
                                   np.eye(3, 3), atol=3*np.finfo(float).eps)

    def test_rotattion_matrix_euler_angles_half_pi(self):
        rx = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        rz = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
        np.testing.assert_allclose(rotation_matrix_euler_angles([0, np.pi/2, 0]), rx, atol=np.finfo(float).eps)
        np.testing.assert_allclose(rotation_matrix_euler_angles([np.pi/2, 0, 0]), rz, atol=np.finfo(float).eps)

    def test_rotate_vecor_null_vector(self):
        v = [0, 0, 0]
        axis = np.random.random(3) * 100
        angle = np.random.random_sample() * 2 * np.pi
        np.testing.assert_allclose(rotate_vector(v, axis, angle), v)

    def test_rotate_vecor_zero_angle(self):
        v = np.random.random(3) * 100
        axis = np.random.random(3) * 100
        angle = 0
        np.testing.assert_allclose(rotate_vector(v, axis, angle), v)

    def test_rotate_vecor_two_pi(self):
        v = np.random.random(3) * 100
        axis = np.random.random(3) * 100
        angle = 2 * np.pi
        np.testing.assert_allclose(rotate_vector(v, axis, angle), v)

    def test_rotate_vector_pi(self):
        angle = np.pi
        v = np.array([1, 0, 0])
        axis = np.array([0, 1, 0])
        np.testing.assert_allclose(rotate_vector(v, axis, angle), -v, atol=np.finfo(float).eps)
        v = np.array([1, 0, 0])
        axis = np.array([0, 0, 1])
        np.testing.assert_allclose(rotate_vector(v, axis, angle), -v, atol=np.finfo(float).eps)
        v = np.array([0, 1, 0])
        axis = np.array([0, 0, 1])
        np.testing.assert_allclose(rotate_vector(v, axis, angle), -v, atol=np.finfo(float).eps)

    def test_to_polar_single_point(self):
        xyz = [1, 0, 0]
        np.testing.assert_allclose(to_polar(xyz), [1, np.pi / 2, 0], atol=np.finfo(float).eps)
        xyz = [0, 1, 0]
        np.testing.assert_allclose(to_polar(xyz), [1, np.pi / 2, np.pi / 2], atol=np.finfo(float).eps)
        xyz = [0, 0, 0]
        np.testing.assert_allclose(to_polar(xyz), [0, 0, 0], atol=np.finfo(float).eps)

    def test_to_polar_wrong_arguments(self):
        xyz = 0
        self.assertRaises(ValueError, to_polar, xyz)
        xyz = [0, 0]
        self.assertRaises(ValueError, to_polar, xyz)
        xyz = [0, 0, 0, 0]
        self.assertRaises(ValueError, to_polar, xyz)
        xyz = [0, 0, 0, 0, 0, 0]
        self.assertRaises(ValueError, to_polar, xyz)

    def test_to_polar_and_back_many_points(self):
        points_num = 100
        xyz = ((np.random.random(points_num * 3) - 0.5) * 200).reshape((points_num, 3))
        rtp = to_polar(xyz)
        np.testing.assert_allclose(to_cartesian(rtp), xyz)

    def test_to_cartesian_single_point(self):
        rtp = [0, 0, 0]
        np.testing.assert_allclose(to_cartesian(rtp), [0, 0, 0], atol=np.finfo(float).eps)
        rtp = [0, 1, 0]
        np.testing.assert_allclose(to_cartesian(rtp), [0, 0, 0], atol=np.finfo(float).eps)
        rtp = [1, 0, 0]
        np.testing.assert_allclose(to_cartesian(rtp), [0, 0, 1], atol=np.finfo(float).eps)
        rtp = [1, np.pi/2, 0]
        np.testing.assert_allclose(to_cartesian(rtp), [1, 0, 0], atol=np.finfo(float).eps)
        rtp = [1, np.pi, 0]
        np.testing.assert_allclose(to_cartesian(rtp), [0, 0, -1], atol=np.finfo(float).eps)

    def test_to_cartesian_wrong_arguments(self):
        rtp = 0
        self.assertRaises(ValueError, to_cartesian, rtp)
        rtp = [0, 0]
        self.assertRaises(ValueError, to_cartesian, rtp)
        rtp = [0, 0, 0, 0]
        self.assertRaises(ValueError, to_cartesian, rtp)
        rtp = [0, 0, 0, 0, 0, 0]
        self.assertRaises(ValueError, to_cartesian, rtp)
        rtp = [-1, 0, 0]
        self.assertRaises(ValueError, to_cartesian, rtp)
        rtp = [1, 2*np.pi, 0]
        self.assertRaises(ValueError, to_cartesian, rtp)
