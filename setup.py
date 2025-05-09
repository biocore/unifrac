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


PREFIX = os.environ.get('PREFIX', "")

base = ["cython >= 0.26", "biom-format", "numpy", "h5py >= 3.3.0",
        "scikit-bio >= 0.6.0", "iow"]

test = ["nose", "flake8"]

all_deps = base + test

# https://stackoverflow.com/a/33308902/379593
if sys.platform == 'darwin':
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.12'


def compile_ssu():
    """Clean and compile the SSU binary"""
    to_link = ["unifrac/task_parameters.hpp",
               "unifrac/api.hpp",
               "unifrac/status_enum.hpp"]

    # clean the target
    cmd = ["rm", "-f"] + to_link
    ret = subprocess.call(cmd)
    if ret != 0:
        raise Exception('Error removing temp unifrac files!')

    for f in to_link:
        # link to files from conda
        cmd = ["ln", "-s", os.environ.get('CONDA_PREFIX') + '/include/' + f,
               "unifrac/"]
        ret = subprocess.call(cmd)
        if ret != 0:
            raise Exception('Error removing linking unifrac files!')


class build_ext(build_ext_orig):
    """Pre-installation for any time an Extension is built"""

    def run(self):
        self.run_compile_ssu()
        super().run()

    def run_compile_ssu(self):
        self.execute(compile_ssu, [], 'Compiling SSU')


if sys.platform == "darwin":
    LINK_ARGS = ['-Wl,' + os.environ.get('CONDA_PREFIX') +
                 '/lib/libssu.so']
else:
    LINK_ARGS = []
COMPILE_ARGS = []

if 'CONDA_PREFIX' in os.environ:
    CONDA_INCLUDES = [os.environ.get('CONDA_PREFIX') + '/include']
else:
    CONDA_INCLUDES = []

USE_CYTHON = os.environ.get('USE_CYTHON', True)
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [Extension("unifrac._api",
                        sources=["unifrac/_api" + ext],
                        extra_link_args=LINK_ARGS,
                        extra_compile_args=COMPILE_ARGS,
                        include_dirs=([np.get_include()] +
                                      CONDA_INCLUDES),
                        libraries=['ssu'])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

with open('README.md') as f:
    long_description = f.read()

setup(
    name="unifrac",
    version="1.5",
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
