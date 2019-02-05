from __future__ import division
import unittest
import numpy as np
from BDSpace.Coordinates import Cartesian
from BDSpace.Field import Field, SuperposedField
from BDSpace.Field import ConstantScalarConservativeField, ConstantVectorConservativeField
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

    def test_constant_scalar_fields(self):
        CField1 = ConstantScalarConservativeField('My field 1', 'My type', np.pi)
        CField2 = ConstantScalarConservativeField('My field 2', 'My type', -np.pi)
        CField3 = ConstantScalarConservativeField('My field 3', 'My type', np.pi)
        SField = SuperposedField('My super Field', [CField1, CField2])
        xyz = np.ones((100, 3), dtype=np.double)
        result = SField.scalar_field(xyz)
        np.testing.assert_allclose(result, np.zeros(100, dtype=np.double))
        SField.fields += [CField3]
        result = SField.scalar_field(xyz)
        np.testing.assert_allclose(result, np.ones(100, dtype=np.double) * np.pi)
        xyz = np.ones((100, 3), dtype=np.double)
        result = SField.vector_field(xyz)
        check = np.zeros((100, 3), dtype=np.double)
        np.testing.assert_allclose(result, check)

    def test_constant_vector_fields(self):
        VField1 = ConstantVectorConservativeField('My field', 'My type', np.array([1.0, 0.0, 0.0], dtype=np.double))
        VField2 = ConstantVectorConservativeField('My field', 'My type', np.array([0.0, 1.0, 0.0], dtype=np.double))
        VField3 = ConstantVectorConservativeField('My field', 'My type', np.array([0.0, 0.0, 1.0], dtype=np.double))
        SField = SuperposedField('My super Field', [VField1, VField2])
        xyz = np.ones((100, 3), dtype=np.double)
        result = SField.vector_field(xyz)
        check = np.zeros((100, 3), dtype=np.double)
        check[:, 0] += np.ones(100, dtype=np.double)
        check[:, 1] += np.ones(100, dtype=np.double)
        np.testing.assert_allclose(result, check)
        SField.fields += [VField3]
        result = SField.vector_field(xyz)
        check[:, 2] += np.ones(100, dtype=np.double)
        np.testing.assert_allclose(result, check)

    def test_constant_vector_fields(self):
        VField1 = ConstantVectorConservativeField('My field', 'My type', np.array([1.0, 0.0, 0.0], dtype=np.double))
        CField1 = ConstantScalarConservativeField('My field 1', 'My type', np.pi)
        SField = SuperposedField('My super Field', [VField1, CField1])
        xyz = np.ones((100, 3), dtype=np.double)
        result = SField.vector_field(xyz)
        check = np.zeros((100, 3), dtype=np.double)
        check[:, 0] += np.ones(100, dtype=np.double)
        np.testing.assert_allclose(result, check)
        result = SField.scalar_field(xyz)
        np.testing.assert_allclose(result, np.pi + np.asarray(VField1.scalar_field(xyz)))
