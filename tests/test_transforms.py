import unittest
import numpy as np
from Space.Coordinates.transforms import adjust_rotation_angles

class TestTransforms(unittest.TestCase):

    def test_adjust_rotation_angles(self):
        angle = [0.1]
        self.assertEqual(adjust_rotation_angles(angle), angle)
