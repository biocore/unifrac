# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.extension import Extension
import numpy as np

import subprocess
import os


SUCPP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'sucpp/')


def compile_ssu():
    """Clean and compile the SSU binary"""
    # clean the target
    subprocess.call(['make', 'clean'], cwd=SUCPP)

    cmd = ['make', 'test']
    ret = subprocess.call(cmd, cwd=SUCPP)
    if ret != 0:
        raise Exception('Error compiling ssu!')

    cmd = ['make', 'main']
    ret = subprocess.call(cmd, cwd=SUCPP)
    if ret != 0:
        raise Exception('Error compiling ssu!')

    cmd = ['make', 'api']
    ret = subprocess.call(cmd, cwd=SUCPP)
    if ret != 0:
        raise Exception('Error compiling ssu!')


class PreBuildCommand(install):
    """Pre-installation for development mode."""
    def run(self):
        self.execute(compile_ssu, [], 'Compiling SSU')
        self.copy_file(os.path.join(SUCPP, 'ssu'),
                       os.path.join(self.install_libbase, 'q2_state_unifrac/'))
        install.run(self)


class PreDevelopCommand(develop):
    """Pre-installation for development mode (i.e. `pip install -e ...`)."""
    def run(self):
        self.execute(compile_ssu, [], 'Compiling SSU')
        self.copy_file(os.path.join(SUCPP, 'ssu'),
                       os.path.join(self.egg_path, 'q2_state_unifrac/'))
        develop.run(self)


import os
USE_CYTHON = os.environ.get('USE_CYTHON', True)
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [Extension("q2_state_unifrac._api",
                        sources=["q2_state_unifrac/_api" + ext, "sucpp/api.cpp"],
                        language="c++",
                        extra_compile_args=["-std=c++11"],
                        extra_link_args=["-std=c++11"],
                        include_dirs=[np.get_include()] + ['sucpp/'],
                        library_dirs=[os.getcwd() + '/sucpp/'],
                        libraries=['ssu'])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name="q2-state-unifrac",
    version="2017.2.0",
    packages=find_packages(),
    author="Daniel McDonald",
    license='BSD-3-Clause',
    author_email="wasade@gmail.com",
    url="https://github.com/wasade/q2-state-unifrac",
    description="An implementation of Strided State UniFrac",
    entry_points={
        "qiime2.plugins":
        ["q2-state-unifrac=q2_state_unifrac.plugin_setup:plugin"]
    },
    ext_modules=extensions,
    cmdclass={'install': PreBuildCommand, 'develop': PreDevelopCommand},
    package_data={
        'q2_state_unifrac.tests': ['data/*', ]}
)
