from __future__ import division
import unittest
import numpy as np
from BDSpace.Field import Field, ConstantScalarConservativeField, ConstantVectorConservativeField


class TestField(unittest.TestCase):

    def setUp(self):
        self.Field = Field('My Field', 'The Type of My Field')

    def test_name(self):
        self.assertEqual(self.Field.name, 'My Field')
        self.Field.name = 'Another name for My Field'
        self.assertEqual(self.Field.name, 'Another name for My Field')

    def test_type(self):
        self.assertEqual(self.Field.type, 'The Type of My Field')
        self.Field.type = 'Another type for My Field'
        self.assertEqual(self.Field.type, 'Another type for My Field')
        with self.assertRaises(TypeError):
            self.Field.type = 3.14

    def test_elements(self):
        self.assertFalse(self.Field.add_element(Field('My another Field', 'The Type of My Field')))
        self.assertFalse(self.Field.remove_element(Field('My another Field', 'The Type of My Field')))

    def test_scalar_field(self):
        xyz = np.ones((100, 3), dtype=np.double)
        result = self.Field.scalar_field(xyz)
        np.testing.assert_allclose(result, np.zeros(100, dtype=np.double))

    def test_vector_field(self):
        xyz = np.ones((100, 3), dtype=np.double)
        result = self.Field.vector_field(xyz)
        np.testing.assert_allclose(result, np.zeros((100, 3), dtype=np.double))

    def test_constant_scalar_field(self):
        CField = ConstantScalarConservativeField('My field', 'My type', np.pi)
        xyz = np.random.random((100, 3))
        result = CField.scalar_field(xyz)
        np.testing.assert_allclose(result, np.ones(100, dtype=np.double) * np.pi)
        result = CField.vector_field(xyz)
        np.testing.assert_allclose(result, np.zeros((100, 3), dtype=np.double))

    def test_constant_vector_field(self):
        VField = ConstantVectorConservativeField('My field', 'My type', np.array([1.0, 0.0, 0.0], dtype=np.double))
        xyz = np.random.random((100, 3))
        result = VField.scalar_field(xyz)
        np.testing.assert_allclose(result, xyz[:, 0])
        result = VField.vector_field(xyz)
        check = np.zeros((100, 3), dtype=np.double)
        check[:, 0] += np.ones(100, dtype=np.double)
        np.testing.assert_allclose(result, check)
