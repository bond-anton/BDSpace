from BDSpace import Space


class Figure(Space):

    def __init__(self, name, coordinate_system=None):
        super(Figure, self).__init__(name, coordinate_system=coordinate_system)

    def __str__(self):
        description = 'Figure: %s\n' % self.name
        description += str(self.coordinate_system)
        return description

    def inner_volume(self):
        """
        calculates volume of the outer shell of the Figure.
        For figures without voids and inner cavities returns the same value as volume function.
        :return: Volume float
        """
        return 0.0

    def external_volume(self):
        """
        calculates volume of the outer shell of the Figure.
        For figures without voids and inner cavities returns the same value as volume function.
        :return: Volume float
        """
        return 0.0

    def volume(self):
        """
        Calculates volume of the Figure.
        :return: Volume (float)
        """
        return self.external_volume() - self.inner_volume()

    def inner_surface_area(self):
        """
        Calculates area of inner surface (zero if no such surface)
        :return: Area of inner surface if any exist.
        """
        return 0.0

    def external_surface_area(self):
        """
        Calculates external surface area of the Figure.
        :return: Area of external surface.
        """
        return 0.0

    def surface_area(self):
        return self.inner_surface_area() + self.external_surface_area()
