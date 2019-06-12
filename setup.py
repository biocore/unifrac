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
from setuptools.command.build_ext import build_ext as build_ext_orig
import numpy as np

import subprocess
import os
import sys


SUCPP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'sucpp/')


PREFIX = os.environ.get('PREFIX', "")


# https://stackoverflow.com/a/33308902/379593
if sys.platform == 'darwin':
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.12'


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


class build_ext(build_ext_orig):

    def run(self):
        self.run_compile_ssu()
        super().run()
    
    def run_compile_ssu(self):
        self.execute(compile_ssu, [], 'Compiling SSU')
        if PREFIX:
            self.copy_file(os.path.join(SUCPP, 'libssu.so'),
                           os.path.join(PREFIX, 'lib/'))


# class PreBuildCommand(install):
#     """Pre-installation for development mode."""
#     def run(self):
#         self.execute(compile_ssu, [], 'Compiling SSU')
#         if PREFIX:
#             self.copy_file(os.path.join(SUCPP, 'libssu.so'),
#                            os.path.join(PREFIX, 'lib/'))
#         install.run(self)
# 
# 
# class PreDevelopCommand(develop):
#     """Pre-installation for development mode (i.e. `pip install -e ...`)."""
#     def run(self):
#         self.execute(compile_ssu, [], 'Compiling SSU')
#         if PREFIX:
#             self.copy_file(os.path.join(SUCPP, 'libssu.so'),
#                            os.path.join(PREFIX, 'lib/'))
#         develop.run(self)


if sys.platform == "darwin":
    LINK_ARGS = ['-Wl,sucpp/libssu.so']
else:
    LINK_ARGS = []

USE_CYTHON = os.environ.get('USE_CYTHON', True)
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [Extension("unifrac._api",
                        sources=["unifrac/_api" + ext,
                                 "sucpp/api.cpp"],
                        language="c++",
                        extra_compile_args=["-std=c++11"],
                        extra_link_args=["-std=c++11"] + LINK_ARGS,
                        library_dirs=[os.path.join(PREFIX, 'lib/')],
                        include_dirs=[np.get_include()] + ['sucpp/',
                            os.path.join(PREFIX, 'bin/')],
                        libraries=['ssu'])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name="unifrac",
    version="0.9.2",
    packages=find_packages(),
    author="Daniel McDonald",
    license='BSD-3-Clause',
    author_email="wasade@gmail.com",
    url="https://github.com/biocore/unifrac",
    description="High performance UniFrac",
    ext_modules=extensions,
#     cmdclass={'install': PreBuildCommand, 
#               'develop': PreDevelopCommand,
#               'build_ext': build_ext},
    cmdclass={'build_ext': build_ext},
    package_data={
        'unifrac.tests': ['data/*', ]}
)
