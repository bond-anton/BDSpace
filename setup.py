from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

from codecs import open
from os import path
import sys
import re


here = path.abspath(path.dirname(__file__))
package_name = 'BDSpace'
version_file = path.join(here, package_name, '_version.py')
with open(version_file, 'rt') as f:
    version_file_line = f.read()
version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(version_re, version_file_line, re.M)
if mo:
    version_string = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in %s.' % (version_file,))

readme_file = path.join(here, 'README.md')
with open(readme_file, encoding='utf-8') as f:
    long_description = f.read()

extensions = [
    Extension(
        'BDSpace.Coordinates.Cartesian',
        ['BDSpace/Coordinates/Cartesian.pyx'],
        depends=['BDSpace/Coordinates/Cartesian.pxd'],
    ),
    Extension(
        'BDSpace.Coordinates._utils',
        ['BDSpace/Coordinates/_utils.pyx'],
        depends=['BDSpace/Coordinates/_utils.pxd'],
    ),
    Extension(
        'BDSpace.Coordinates.transforms',
        ['BDSpace/Coordinates/transforms.pyx'],
        depends=['BDSpace/Coordinates/transforms.pxd'],
    ),
]

setup(
    name=package_name,
    version=version_string,

    description='3D space positioning and motion',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://github.com/bond-anton/BDSpace',

    author='Anton Bondarenko',
    author_email='bond.anton@gmail.com',

    license='Apache Software License',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='3D coordinate Space paths trajectory',

    packages=find_packages(exclude=['demo', 'tests', 'docs', 'contrib', 'venv']),
    ext_modules=cythonize(extensions, compiler_directives={'language_level': sys.version_info[0]}),
    package_data={
        'Coordinates.Cartesian_c': ['Cartesian.pxd', '_utils.pxd'],
        'Coordinates.transforms_c': ['transforms.pxd'],
        # 'Schottky.Metal': ['Schottky/Metal/__init__.pxd'],
    },
    install_requires=['numpy', 'Cython',
                      'BDQuaternions>=0.1.2'],
    test_suite='nose.collector',
    tests_require=['nose']
)
