from Space import Space


class Figure(Space):

    def __init__(self, name, coordinate_system=None):
        super(Figure, self).__init__(name, coordinate_system=coordinate_system)

    def __str__(self):
        return 'Figure: ' + self.name

    def volume(self):
        """
        Calculates volume of the Figure.
        :return: Volume (float)
        """
        return 0.0

    def external_volume(self):
        """
        calculates volume of the outer shell of the Figure.
        For figures without voids and inner cavities returns the same value as volume function.
        :return: Volume float
        """
        return self.volume()
