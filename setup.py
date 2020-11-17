# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
import numpy as np

import subprocess
import os
import sys


SUCPP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'sucpp/')


PREFIX = os.environ.get('PREFIX', "")

base = ["cython >= 0.26", "biom-format", "numpy", "h5py >= 2.7.0",
        "scikit-bio >= 0.5.1"]

test = ["nose", "flake8"]

all_deps = base + test

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
    """Pre-installation for any time an Extension is built"""

    def run(self):
        self.run_compile_ssu()
        super().run()

    def run_compile_ssu(self):
        self.execute(compile_ssu, [], 'Compiling SSU')
        if PREFIX:
            self.copy_file(os.path.join(SUCPP, 'libssu.so'),
                           os.path.join(PREFIX, 'lib/'))


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
                        include_dirs=[np.get_include()] + ['sucpp/'],
                        libraries=['ssu'])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

with open('README.md') as f:
    long_description = f.read()

setup(
    name="unifrac",
    version="0.20.1",
    packages=find_packages(),
    author="Daniel McDonald",
    license='BSD-3-Clause',
    author_email="wasade@gmail.com",
    url="https://github.com/biocore/unifrac",
    description="High performance phylogenetic diversity calculations",
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=extensions,
    install_requires=base,
    extras_require={'test': test, 'all': all_deps},
    cmdclass={'build_ext': build_ext},
    package_data={
        'unifrac.tests': ['data/*', ]}
)
