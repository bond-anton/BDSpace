import unittest
from BDSpace import Space


class TestSpace(unittest.TestCase):

    def setUp(self):
        self.solar_system = Space('Solar System')
        mercury = Space('Mercury')
        venus = Space('Venus')
        earth = Space('Earth')
        mars = Space('Mars')
        self.solar_system.add_element(mercury)
        self.solar_system.add_element(venus)
        self.solar_system.add_element(earth)
        self.solar_system.add_element(mars)

    def space_constructor(self):
        self.assertRaises(ValueError, Space, 'xxx', '0')

    def test_add_remove_subspaces(self):
        self.assertFalse(self.solar_system.add_element(self.solar_system))
        self.assertRaises(TypeError, self.solar_system.add_element, 'Any object except BDSpace')
        earth = self.solar_system.elements['Earth']
        moon = Space('Moon')
        earth.add_element(moon)
        lunohod = Space('Lunohod')
        moon.add_element(lunohod)
        mars = self.solar_system.elements['Mars']
        phobos = Space('Phobos')
        deimos = Space('Deimos')
        mars.add_element(phobos)
        mars.add_element(deimos)
        self.assertFalse(self.solar_system.add_element(lunohod))
        self.assertRaises(TypeError, self.solar_system.remove_element, 'Any object, not BDSpace')

    def test_add_same_name_subspaces(self):
        mars = self.solar_system.elements['Mars']
        count = 105
        for i in range(count):
            phobos = Space('Phobos')
            mars.add_element(phobos)
        self.assertEqual(phobos.name, 'Phobos %d' % (count-1))

    def test_basis_in_global_coordinates(self):
        print('Basis in GCS:', self.solar_system.basis_in_global_coordinate_system())
