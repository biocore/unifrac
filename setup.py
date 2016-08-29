# ----------------------------------------------------------------------------
# Copyright (c) 2013--, BP development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import sys
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_py import build_py

# TODO: actually resolve where bp/_bp.pxd is for BP headers
sys.path.insert(0, os.path.abspath("../improved-octo-waddle"))


import numpy as np

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Programming Language :: Python
    Programming Language :: Python :: 3.5
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = """An state-based unifrac"""

from Cython.Compiler.Options import directive_defaults

USE_CYTHON = os.environ.get('USE_CYTHON', True)
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("q2_su._su",
                        ["q2_su/_su" + ext],
                        include_dirs=['BitArray/'],
                        library_dirs=['BitArray/'],
                        libraries=['bitarr'],
                            #    extra_compile_args=['-fopenmp'],
                             #           extra_link_args=['-fopenmp'],
                        define_macros=[('CYTHON_TRACE', '0')]),
              ]

extensions.extend([Extension("q2_su._ba",
                            ["q2_su/_ba" + ext],
                            include_dirs=['BitArray/'],
                             library_dirs=['BitArray/'],
                             libraries=['bitarr'],
                        define_macros=[('CYTHON_TRACE', '0')])])

extensions.extend([Extension("q2_su.test_su",
                            ["q2_su/test_su" + ext],
                            include_dirs=['BitArray/'],
                             library_dirs=['BitArray/'],
                             libraries=['bitarr'])])

extensions.extend([Extension("q2_su.foo",
                            ["q2_su/foo" + ext],
                            include_dirs=['BitArray/'],
                             library_dirs=['BitArray/'],
                             libraries=['bitarr'])])

#from Cython.Compiler.Options import directive_defaults

#directive_defaults['linetrace'] = True
#directive_defaults['binding'] = True

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, annotate=True)


import subprocess
bitarr = os.path.join(os.path.abspath(__file__).rsplit('/', 1)[0], 'BitArray')


class BitArrayInstall(build_py):
    def run(self):
        subprocess.run(['make', '-C', bitarr, 'libbitarr.a'])
        build_py.run(self)


setup(
    name="q2-state-unifrac",
    version='0.0.1',
    packages=find_packages(),
    install_requires=['iow >= 0.1.0, < 0.2',
                      'qiime >= 2.0.2', 'q2-types'],
    author="Daniel McDonald",
    author_email="mcdonadt@colorado.edu",
    description="State UniFrac implementation",
    include_dirs=[np.get_include(), 'BitArray/'],
    license='BSD-3-Clause',
    ext_modules=extensions,
    url="https://github.com/wasade/q2-state-unifrac",
    entry_points={
        'qiime.plugins':
        ['q2-state-unifrac=q2_su.plugin_setup:plugin']},
    cmdclass={'build_py': BitArrayInstall})

