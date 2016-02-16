from distutils.core import setup

setup(
    name='Space',
    version='0.0.1',
    packages=['Space', 'Space.Coordinates',
              'Space.Curve', 'Space.Figures',
              'Space.Pathfinder'],
    url='https://github.com/bond-anton/Space',
    license='Apache Software License',
    author='Anton Bondarenko',
    author_email='bond.anton@gmail.com',
    description='3D space positioning and motion'
)
