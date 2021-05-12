import sys
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize

from codecs import open
from os import path, remove
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
        'BDSpace.Space',
        ['BDSpace/Space.pyx'],
        depends=['BDSpace/Space.pxd'],
    ),
    Extension(
        'BDSpace.Coordinates.Cartesian',
        ['BDSpace/Coordinates/Cartesian.pyx'],
        depends=['BDSpace/Coordinates/Cartesian.pxd'],
    ),
    Extension(
        'BDSpace.Coordinates.transforms',
        ['BDSpace/Coordinates/transforms.pyx'],
        depends=['BDSpace/Coordinates/transforms.pxd'],
    ),
    Extension(
        'BDSpace.Field.Field',
        ['BDSpace/Field/Field.pyx'],
        depends=['BDSpace/Field/Field.pxd'],
    ),
    Extension(
        'BDSpace.Field.SphericallySymmetric',
        ['BDSpace/Field/SphericallySymmetric.pyx'],
        depends=['BDSpace/Field/SphericallySymmetric.pxd'],
    ),
    Extension(
        'BDSpace.Field.SuperposedField',
        ['BDSpace/Field/SuperposedField.pyx'],
        depends=['BDSpace/Field/SuperposedField.pxd'],
    ),
    Extension(
        'BDSpace.Curve._helpers',
        ['BDSpace/Curve/_helpers.pyx'],
        depends=['BDSpace/Curve/_helpers.pxd'],
    ),
    Extension(
        'BDSpace.Curve.Parametric',
        ['BDSpace/Curve/Parametric.pyx'],
        depends=['BDSpace/Field/Field.pxd'],
    ),
    Extension(
        'BDSpace.Field.CurveField',
        ['BDSpace/Field/CurveField.pyx'],
        depends=['BDSpace/Field/CurveField.pxd'],
    ),
]

copt = {'msvc': [],
        'mingw32': [],
        'unix': []}
lopt = {'mingw32': [],
        'unix': []}


# check whether compiler supports a flag
def has_flag(compiler, flagname):
    import tempfile
    from distutils.errors import CompileError
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            res = compiler.compile([f.name], extra_postargs=[flagname])
            for item in res:
                remove(item)
        except CompileError:
            return False
    return True


# filter flags, returns list of accepted flags
def flag_filter(compiler, flags):
    result = []
    for flag in flags:
        if has_flag(compiler, flag):
            result.append(flag)
    return result


class CustomBuildExt(build_ext):
    def build_extensions(self):
        c = self.compiler.compiler_type
        print('Compiler:', c)
        opts = flag_filter(self.compiler, copt.get(c, []))
        lopts = flag_filter(self.compiler, lopt.get(c, []))
        for e in self.extensions:
            e.extra_compile_args = opts
            e.extra_link_args = lopts
        build_ext.build_extensions(self)


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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],

    keywords='3D coordinate Space paths trajectory',

    packages=find_packages(exclude=['demo', 'tests', 'docs', 'contrib', 'venv']),
    ext_modules=cythonize(extensions, compiler_directives={'language_level': sys.version_info[0]}),
    package_data={
        'BDSpace': ['*.pxd'],
        'BDSpace.Coordinates': ['*.pxd'],
        'BDSpace.Field': ['*.pxd'],
        'BDSpace.Curve': ['*.pxd'],
    },
    install_requires=['numpy', 'scipy',
                      'BDQuaternions>=0.2.11',
                      'BDMesh>=0.2.11'],
    cmdclass={'build_ext': CustomBuildExt},
    zip_safe=False,
)
