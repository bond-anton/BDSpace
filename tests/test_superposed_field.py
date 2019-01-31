from __future__ import division
import unittest
import numpy as np
from BDSpace.Coordinates import Cartesian
from BDSpace.Field import Field, SuperposedField
from BDSpace.Coordinates._utils import check_points_array
from BDQuaternions import Conventions


class TestField(unittest.TestCase):

    def setUp(self):
        self.Field = Field('My Field', 'The Type of My Field')
        self.SField = SuperposedField('My super Field', [])

    def test_name(self):
        self.assertEqual(self.SField.name, 'My super Field')
        self.SField.name = 'Another name for My super Field'
        self.assertEqual(self.SField.name, 'Another name for My super Field')

    def test_type(self):
        self.assertIsNone(self.SField.type)
        self.SField.type = 'Another type for My super Field'
        self.assertEqual(self.SField.type, 'Another type for My super Field')
        with self.assertRaises(TypeError):
            self.SField.type = 3.14

    def test_elements(self):
        with self.assertRaises(NotImplementedError):
            self.SField.add_element(Field('My another Field', 'The Type of My Field'))
        with self.assertRaises(NotImplementedError):
            self.SField.remove_element(Field('My another Field', 'The Type of My Field'))

    def test_scalar_field(self):
        xyz = np.ones((100, 3), dtype=np.double)
        result = self.SField.scalar_field(xyz)
        np.testing.assert_allclose(result, np.zeros(100, dtype=np.double))

    def test_vector_field(self):
        xyz = np.ones((100, 3), dtype=np.double)
        result = self.SField.vector_field(xyz)
        np.testing.assert_allclose(result, np.zeros((100, 3), dtype=np.double))
